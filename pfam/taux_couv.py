#!/usr/bin/env python
# -*- coding: utf-8 -*-


''' Pour chaque proteome:
    calcul le taux de couverture en residues
    et le taux de couverture en domaine.
    renvoit un histogramme du taux de coubverture

Parametres:

-i dir (.fasta) dir (.pfam) -o [histogramme]


usage
==============
python taux_couv.py
'''
import os, sys,re
import matplotlib
import matplotlib.pyplot as plt
import numpy as np


def read_fasta(path_proteome):
    ''''Description de la fonction : lit un fichier .fasta et renvoit un dict avec key = proteome_name et val = dict('protein_name': len_sequence)
    '''
    with open(path_proteome) as f:
        protein_dict = {}
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                name_protein = line[1:].split()[0]

                if name_protein in protein_dict:
                    print("Error, protein {} already exists".format(name_protein))
                    sys.exit(1)

                protein_dict[name_protein] = 0
            else:
                len_seq = len(line.strip())
                protein_dict[name_protein] += len_seq
    return protein_dict


def read_pfam_outp (path_pfam_ouf):
    '''Description de la fonction : renvoit un dict pour :key = proteome_name, value = dict('protein_name':[(start, stop)]) du domaine trouve sur pfam
    :NB: un protein_id est present en fonction du nombre de domaines qu'il porte '''
    with open(path_pfam_ouf) as pf:
        pfam_dict = {}
        for line in pf:
            if line.strip() and line[0] != '#':
                line = line.rstrip()
                tmp = line.split()
                seq_id = tmp[0]
                start = int(tmp[1]) - 1
                stop  = int(tmp[2])
                pos = start, stop
                if seq_id not in pfam_dict:
                    pfam_dict[seq_id] = []
                pfam_dict[seq_id].append(pos)
    return pfam_dict

def long_couvert(pfam_dict):
    '''Description de la fonction: renvoit la longeur totale couvert par les domaines sur la proteine
    :NB: permet de passer de la redondance des protein_id de la fonction precedente à un seul id par protein et comme value : la couverture totale des domaine/protein
    '''
    len_domains = {}
    for seq_id in pfam_dict:
        values = pfam_dict[seq_id]
        for start, stop in values:
            ld = stop - start
            len_domains[seq_id] = len_domains.get(seq_id, 0) + ld
    return len_domains

def taux_couverture(directory, liste_fasta, pfam_out, list_pfam_out):
    '''Description de la fonction:
    prend comme entrees la localisation et la liste des fichiers fasta et pfam.
    execute les deux fonctions (read_pfam_outp) et (long_couvert)
    renvoit:
    Taux de couverture pour les proteomes en residus (len_domain_concat/len_totale_proteome)*100
    Et le taux de couverture de chaque proteome en domain (ayant au moins un domaine present) = (nmb_domain_present/nombre_total_protein)*100
    '''

    all_prot, all_dom = dict(), dict()

    for filename in liste_fasta:
        if filename.endswith('.fasta'):
            path_proteome = os.path.join(directory, filename)
            path_pfam_ouf = os.path.join(pfam_out, filename+".pfam")
            proteome_name = filename.split('.')[0]

            fasta = read_fasta(path_proteome)
            table_pfam_dict = read_pfam_outp(path_pfam_ouf)
            len_domainsTot = long_couvert(table_pfam_dict)


            nmb_domain1 = 0
            nmbr_prot = 0
            len_proteome_domain = 0
            len_proteome_fasta  = 0
            #Calcul de la longueur totale pour chaque  proteome
            for seq_id in fasta:
                nmbr_prot =nmbr_prot+ 1
                len_fasta = fasta[seq_id]
                len_proteome_fasta = len_proteome_fasta + int(len_fasta)
            for seq_id in len_domainsTot:
                nmb_domain1 = nmb_domain1+1
                ld = len_domainsTot[seq_id]
                len_proteome_domain = len_proteome_domain + int(ld)


            taux_prot = (len_proteome_domain/len_proteome_fasta)*100
            taux_dom1 = (nmb_domain1/nmbr_prot)*100

            all_prot[proteome_name] = taux_prot
            all_dom[proteome_name] = taux_dom1

    return all_prot, all_dom

def plotting(all_prot, all_dom):
    '''renvoit un histogramme pour le taux de couverture en residue
    et pour le taux de couverture en domaines de chaque proteome
    '''

    all_couv = []
    t_prot_dom_name=()
    for proteome_name in all_prot:


        taux_prt = all_prot[proteome_name]
        taux_dom = all_dom[proteome_name]
        t_prot_dom_name = taux_prt, taux_dom, proteome_name

        all_couv.append(t_prot_dom_name)

    all_couv.sort()

    list_de_taux = []
    list_domain1 = []
    list_prot_name = []
    for data_protein, data_dom, proteome_lab in all_couv:

        list_de_taux.append(data_protein)
        list_domain1.append(data_dom)
        list_prot_name.append(proteome_lab)
    # pour creer des partitions
    # for i in range(0, len(list_prot_name), 30):
    #     data_protein = list_de_taux[i : i +30]
    #     data_domain = list_domain1[i : i +30]
    #     name_proteome = list_prot_name[i : i +30]

        #histogramme
    fig, ax = plt.subplots()
    N = len(list_prot_name)
    ind = np.arange(N)
    width = 0.35

    prot = ax.bar(ind, list_de_taux, width=0.2, alpha=0.5, color='b', label= 'Taux de couverture en residues (%)')
    dom = ax.bar(ind+width, list_domain1, width=0.2, alpha=0.5, color='g', label= 'Taux de couverture en domaines (%)')
    #
    ax.set_xlim(-width,len(ind)+width)
    ax.set_xticks(ind+width)
    ax.set_xticklabels (list_prot_name, rotation='vertical', fontsize=10)
    for name in ax.set_xticklabels(list_prot_name, rotation= 'vertical', fontsize=10):
        if re.search('Synechococcus_sp_PCC_6312', str(name)) :
            name.set_color('green')
        elif re.search('Synechococcus_calcipolaris' , str(name)):
            name.set_color('green')
        elif re.search('Thermosynechococcus_elongatus_BP1' , str(name)):
            name.set_color('green')


        elif re.search('Gloeomargarita_lithophora' , str(name)):
            name.set_color('red')
        elif  re.search ('Cyanothece_sp_PCC_7425', str(name)):
            name.set_color('red')
        elif  re.search ('Chroococcidiopsis_thermalis_PCC_7203', str(name)):
            name.set_color('red')
    plt.grid(True)
    plt.ylabel('Taux de couverture', fontsize=15, color='g', alpha=0.8)
    plt.xlabel('Proteomes', fontsize=15, color='b', alpha=0.8)
    ax.set_ylim(0,100)
    plt.title(u"Taux de couverture des proteomes par Pfam", fontsize=17, fontdict={'family': 'monospace'})
    ax.legend(loc='upper left', fontsize=10)
    # # plt.savefig(filename+u'boxplotProteome.png')
    plt.show()


if __name__ == '__main__':

    directory = '/home/issa/Documents/stage/init_data/proteomes'
    liste_fasta = os.listdir(directory)

    pfam_out = '/home/issa/pfam/pfamScan/results_pfamscan'
    list_pfam_out = os.listdir(pfam_out)

    proteome_couvert, domain1_couvert = taux_couverture(directory, liste_fasta, pfam_out, list_pfam_out)
    plotting(proteome_couvert, domain1_couvert)
