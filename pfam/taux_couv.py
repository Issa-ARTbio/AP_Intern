#!/usr/bin/env python
# -*- coding: utf-8 -*-


''' Pour chaque proteome:
    calcul le taux de couverture des domaines
    et la frequence des proteines avec au moins 1 domaine present

Parametres:

-i dir (.fasta) dir (.pfam) -o []


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
    :NB: permet de passer de la redondance des protein_id de la fonction precedente à un seul id par protein et comme value : la couverture total des domaine/protein
    '''
    len_domains = {}
    for seq_id in pfam_dict:
        values = pfam_dict[seq_id]
        for start, stop in values:
            ld = stop - start
            len_domains[seq_id] = len_domains.get(seq_id, 0) + ld
    return len_domains

def taux_couverture(directory, liste_fasta, pfam_out, list_pfam_out):
    '''Description de la fonction: prend comme entrees deux dico avec les mm keys = proteome name
    dans proteome_dict : values = dict {'protein_name': seq_len}
    dans all_domains: values = dict {'protein_name': len_domain}
    renvoit:
    Taux de couvertures pour les proteomes
    Et la frequence de protein(ayant au moins un domaine) sur le nombre total de protein
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
    list_de_taux = []
    list_domain1 = []
    list_prot_name = []

    for proteome_name in all_prot:

        list_prot_name.append(proteome_name)

        taux_prt = all_prot[proteome_name]
        list_de_taux.append(taux_prt)

        taux_dom = all_dom[proteome_name]
        list_domain1.append(taux_dom)

    fig, ax = plt.subplots()
    N = 111
    ind = np.arange(N)
    width = 0.35


    zipper = list(zip(list_de_taux, list_prot_name))
    zipper_sort = sorted (zipper)
    data_prot, proteome_lab = zip(*zipper_sort)
    prot = ax.bar(ind, data_prot, width=0.2, alpha=0.5, color='b', label= 'Taux de couverture en residues (%)')


    zipper = list(zip(list_domain1, list_prot_name))
    zipper_sort = sorted (zipper)
    data_dom, proteome_nam = zip(*zipper_sort)
    dom = ax.bar(ind+width, data_dom, width=0.2, alpha=0.5, color='g', label= 'Taux de couverture en domaines (%)')

    ax.set_xlim(-width,len(ind)+width)

    ax.set_xticks(ind+width)
    ax.set_xticklabels (proteome_lab, rotation='vertical', fontsize=8)
    for i in ax.set_xticklabels(proteome_lab, rotation= 'vertical'):
        if re.search('Synechococcus_sp_PCC_6312', str(i)) :
            i.set_color('green')
        elif re.search('Synechococcus_calcipolaris' , str(i)):
            i.set_color('green')
        elif re.search('Thermosynechococcus_elongatus_BP1' , str(i)):
            i.set_color('green')


        elif re.search('Gloeomargarita_lithophora' , str(i)):
            i.set_color('red')
        elif  re.search ('Cyanothece_sp_PCC_7425', str(i)):
            i.set_color('red')
        elif  re.search ('Chroococcidiopsis_thermalis_PCC_7203', str(i)):
            i.set_color('red')
    plt.grid(True)
    plt.ylabel('Taux de couverture', fontsize=15, color='g', alpha=0.8)
    plt.xlabel('Proteomes', fontsize=15, color='b', alpha=0.8)
    plt.title(u"Taux de couverture des proteomes par Pfam", fontsize=17, fontdict={'family': 'monospace'})
    ax.legend(loc='upper left')
    # plt.savefig(filename+u'boxplotProteome.png')
    plt.show()


if __name__ == '__main__':

    directory = '/home/issa/Documents/STAGE/Init_data/proteomes'
    liste_fasta = os.listdir(directory)

    pfam_out = '/home/issa/pfam/pfamScan/results_pfamscan'
    list_pfam_out = os.listdir(pfam_out)

    proteome_couvert, domain1_couvert = taux_couverture(directory, liste_fasta, pfam_out, list_pfam_out)
    plotting(proteome_couvert, domain1_couvert)
