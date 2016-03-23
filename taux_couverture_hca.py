#!/usr/bin/env python
# -*- coding: utf-8 -*-

''' Pour chaque proteome:
    calcul le taux de couverture en residues
    et le taux de couverture en domaine pyHCA.
    renvoit un histogramme du taux de coubverture

Parametres:

-i dir (.fasta) dir (.pfam) -o [histogramme]


usage
==============
python taux_couverture_hca.py
'''

import os, sys,re
from collections import defaultdict
import matplotlib
import matplotlib.pyplot as plt
import numpy as np


def read_fasta(fasta_in):
    ''''Description de la fonction : lit un fichier .fasta et renvoit un dict avec key = proteome_name et val = dict('protein_name': len_sequence)
    '''
    with open(fasta_in) as f:
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
        # for x in protein_dict:
        #     print(x, protein_dict[x])
    return protein_dict


def read_pyHCA_outF (hca_in):
    """return domains match in sequences by pyHCA program
    defaultdict(list) : key = protein_id , value = tpl(domain_start, domain_end, domain_lenght)
    """
    lenght_sequences = defaultdict(list)
    domain_hca = defaultdict(list)

    with open(hca_in, 'r') as hca:
        for line in hca:
            line = line.strip()
            domain = 0
            if line.startswith('>'):
                header = line[1:].split()[0]
                seq_lenght = line[1:].split()[1]
                lenght_sequences[header] = (seq_lenght)
            if line.startswith('domain'):
                domain_start  = line.strip().split()[1]
                domain_end    = line.strip().split()[2]
                domain_lenght = int(domain_end) - int(domain_start)
                domain   = domain + domain_lenght
                domain_hca[header].append(domain)

    domain_couv, proteome_lenght = {}, {}
    for protein in domain_hca:
        values = domain_hca[protein]
        dom_lenght = sum(values)
        if protein == 'Pse6802_5039':
            print(protein, domain_hca[protein])
        domain_couv[protein] = dom_lenght

    for protein in proteome_lenght:
        seq_lenght = lenght_sequences[protein]
        lenght_total = sum(seq_lenght)
        proteome_lenght[protein] = lenght_total

    return domain_couv, lenght_sequences


def read_cdd_outF (cdd_in):
    """return domain match in sequences by CDD scan
    defaultdict(list) : key = protein_id , value = tpl(domain_start, domain_end, domain_lenght)
    """
    domain_cdd = defaultdict(set)
    with open(cdd_in, 'r') as cdd:
        for line in cdd:
            line = line.strip()
            domain = 0
            if line.startswith('Q#'):
                element = line.split()
                header = element[2][1:]
                domain_start = element[6]
                domain_end   = element[7]
                domain_lenght= int(domain_end) - int(domain_start)
                print(header, domain_start, domain_end)
                domain = domain + domain_lenght
                domain_cdd[header].add(domain)
    all_domain_cdd = {}
    for protein in domain_cdd:
        values = domain_cdd[protein]
        dom_lenght = sum(values)
        all_domain_cdd[protein] = dom_lenght
    return all_domain_cdd



def taux_couverture (directory, liste_fasta, hca_out, cdd_out):
    print('running ...')
    # For each file in the directory, execute the previous fonctions
    all_prot, all_dom = dict(), dict()
    for filename in liste_fasta:
        if filename.endswith('.fasta'):
            fasta_in = os.path.join(directory, filename)
            hca_in   = os.path.join(hca_out, filename+".hca")
            cdd_in   = os.path.join(cdd_out, filename+".cdd")
            proteome_name = filename.split('.')[0]

            fasta = read_fasta(fasta_in)
            hca, seq_len = read_pyHCA_outF (hca_in)
            cdd = read_cdd_outF(cdd_in)
            # print (filename, len(hca), len(seq_len))
            nmb_domain = 0
            nmbr_prot = 0
            len_proteome_domain = 0
            len_proteome_fasta  = 0
            #Calcul de la longueur totale pour chaque  proteome
            for seq_id in fasta:
                nmbr_prot+=1
                len_fasta = fasta[seq_id]
                len_proteome_fasta = len_proteome_fasta + int(len_fasta)
            for seq_id in hca:
                nmb_domain+=1
                ld = hca[seq_id]
                len_proteome_domain = len_proteome_domain + int(ld)

            # print(hca)
            taux_prot = (len_proteome_domain/len_proteome_fasta)*100
            domaine_count = (nmb_domain/nmbr_prot)*100
            # print(filename, '==', domaine_count, '==', taux_prot)
            all_prot[proteome_name] = taux_prot
            all_dom[proteome_name] = domaine_count
    return all_prot, all_dom

def plotting (all_prot, all_dom):

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
        # print (proteome_lab, '==', data_protein, '==', data_dom)
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
    ax.set_ylim(0,110)
    plt.title(u"Taux de couverture en domaines des proteomes par HCA", fontsize=17, fontdict={'family': 'monospace'})
    ax.legend(loc='upper left', fontsize=10)
    # # plt.savefig(filename+u'boxplotProteome.png')
    plt.show()

if __name__ == '__main__':



    directory = '/home/issa/Documents/stage/initial_data/proteomes'
    liste_fasta = os.listdir(directory)

    hca_out = '/home/issa/Documents/stage/CDD_pyHCA/pyHCA/analyse/result_hca/'
    cdd_out = '/home/issa/Documents/stage/CDD_pyHCA/CDD/analyse/results_CDD/'



    protein, domain = taux_couverture (directory, liste_fasta, hca_out, cdd_out)
    plotting(protein, domain)
    print('succefully done !!!')
