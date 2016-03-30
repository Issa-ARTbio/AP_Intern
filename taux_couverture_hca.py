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
    return protein_dict


def read_pyHCA_outF (hca_in):
    """return domains match in sequences by pyHCA program
    defaultdict(list) : key = protein_id , value = tpl(domain_start, domain_end, domain_lenght)
    """
    domain_hca = defaultdict(set)

    with open(hca_in, 'r') as hca:
        for line in hca:
            line = line.strip()
            if line.startswith('>'):
                header = line[1:].split()[0]
                seq_lenght = line[1:].split()[1]
            if line.startswith('domain'):
                domain_start  = int(line.strip().split()[1])
                domain_end    = int(line.strip().split()[2])
                domain_lenght = domain_end - domain_start
                if domain_lenght >= 30:
                    domain_hca[header].add(domain_lenght)

    domain_couv = {}
    for protein in domain_hca:
        values = domain_hca[protein]
        dom_lenght = sum(values)
        domain_couv[protein] = dom_lenght


    return domain_couv


def read_cdd_outF (cdd_in, protein_dict):
    """return domain match in sequences by CDD scan
    defaultdict(list) : key = protein_id , value = tpl(domain_start, domain_end, domain_lenght)
    """
    domain_cdd = defaultdict(set)
    with open(cdd_in, 'r') as cdd:
        for line in cdd:
            line = line.strip()
            if line.startswith('Q#'):
                domains_list = []
                element = line.split('\t')
                header = element[0].split()[2][1:]
                if header.find('['):
                    header= header.split('[')[0]
                domain_start = element[3]
                domain_start = int(domain_start) - 1
                domain_end   = int(element[4])
                domain_lenght= domain_end - domain_start
                domain = (domain_start, domain_end, domain_lenght)
                domain_cdd[header].add(domain)

    dico_bin = {}
    for protein in protein_dict:
        lenght_seq = protein_dict[protein]
        if protein in domain_cdd:
            liste_of_position = [0]*lenght_seq
            cdd_dom = domain_cdd[protein]
            for start, stop, lenght in cdd_dom:
                for i in range (int(start), int(stop)):
                    liste_of_position[i] = 1
                    total_pos = sum(liste_of_position)
            dico_bin[protein] = total_pos

    return dico_bin

def taux_couverture (directory, liste_fasta, hca_out, cdd_out):
    print('running ...')
    # For each file in the directory, execute the previous fonctions
    couverture_residues_hca, couverture_domaines_hca, couverture_domaines_cdd, couverture_residues_cdd = dict(), dict(), dict(), dict()
    for filename in liste_fasta:
        if filename.endswith('.fasta'):
            fasta_in = os.path.join(directory, filename)
            hca_in   = os.path.join(hca_out, filename+".hca")
            cdd_in   = os.path.join(cdd_out, filename+".cdd")
            proteome_name = filename.split('.')[0]

            fasta = read_fasta(fasta_in)
            hca = read_pyHCA_outF (hca_in)
            cdd = read_cdd_outF(cdd_in, fasta)

            # Count the number of domain in hca and in cdd and the number of proteins in the fasta
            nmb_domain = 0
            nmbr_prot  = 0
            nb_dom_cdd = 0

            # Count the number of residues
            len_proteome_hca    = 0
            len_proteome_fasta  = 0
            len_proteome_cdd    = 0

            #Compute the total leght for each proteome
            for seq_id in fasta:
                nmbr_prot+=1
                len_fasta = int(fasta[seq_id])
                len_proteome_fasta += len_fasta
            for seq_id in hca:
                nmb_domain+=1
                ld = int(hca[seq_id])
                len_proteome_hca =  len_proteome_hca + ld
            for seq_id in cdd:
                nb_dom_cdd+=1
                lenght_cdd = int(cdd[seq_id])
                len_proteome_cdd = len_proteome_cdd + lenght_cdd


            taux_couverture_residues = (len_proteome_hca/len_proteome_fasta)*100
            taux_couverture_domains  = (nmb_domain/nmbr_prot)*100

            tx_residues_cdd = (len_proteome_cdd/len_proteome_fasta)*100
            tx_domains_cdd  = (nb_dom_cdd/nmbr_prot)*100

            couverture_residues_hca[proteome_name] = taux_couverture_residues
            couverture_domaines_hca[proteome_name] = taux_couverture_domains

            couverture_residues_cdd[proteome_name] = tx_residues_cdd
            couverture_domaines_cdd[proteome_name] = tx_domains_cdd
            print(proteome_name, taux_couverture_residues)
    return couverture_residues_hca, couverture_domaines_hca, couverture_residues_cdd, couverture_domaines_cdd

def plotting (couverture_residues_hca, couverture_domaines_hca, couverture_residues_cdd, couverture_domaines_cdd):

    all_couverture = []
    my_set =()
    for proteome_name in couverture_domaines_hca:


        residues_hca = couverture_residues_hca[proteome_name]
        domain_hca = couverture_domaines_hca[proteome_name]

        residues_cdd = couverture_residues_cdd[proteome_name]
        domain_cdd = couverture_domaines_cdd[proteome_name]

        my_set = residues_hca, domain_hca, residues_cdd, domain_cdd, proteome_name

        all_couverture.append(my_set)

    all_couverture.sort()

    res_hca, dom_hca, res_cdd, dom_cdd, list_proteome_name  = [], [], [], [], []

    for r_hca, d_hca, r_cdd, d_cdd, proteome in all_couverture:
        res_hca.append(r_hca)
        dom_hca.append(d_hca)
        res_cdd.append(r_cdd)
        dom_cdd.append(d_cdd)
        list_proteome_name.append(proteome)

    # pour creer des partitions
    # for i in range(0, len(list_prot_name), 30):
    #     data_protein = list_de_taux[i : i +30]
    #     data_domain = list_domain1[i : i +30]
    #     name_proteome = list_prot_name[i : i +30]

        #histogramme
    fig, ax = plt.subplots()
    N = len(list_proteome_name)
    ind = np.arange(N)
    width = 0.35

    hca_residues_plot = ax.bar(ind, res_hca, width=0.2, alpha=0.5, color='k', label= 'Taux de couverture en residues HCA (%)')
    # hca_domains_plot = ax.bar(ind, dom_hca, width=0.2, alpha=0.5, color='k', label= 'Taux de couverture en domaines HCA (%)')

    cdd_res_plot = ax.bar(ind+width, res_cdd, width=0.2, alpha=0.5, color='b', label= 'Taux de couverture en residues CDD (%)')
    # cdd_dom_plot = ax.bar(ind+width, dom_cdd, width=0.2, alpha=0.5, color='b', label= 'Taux de couverture en domaines CDD(%)')

    ax.set_xlim(-width,len(ind)+width)
    ax.set_xticks(ind+width)
    ax.set_xticklabels (list_proteome_name, rotation='vertical', fontsize=10)
    for name in ax.set_xticklabels(list_proteome_name, rotation= 'vertical', fontsize=10):
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
    plt.title(u"Taux de couverture en residues des proteomes par HCA et CDD", fontsize=17, fontdict={'family': 'monospace'})
    ax.legend(loc='upper left', fontsize=10)
    plt.show()

if __name__ == '__main__':



    directory = '/home/issa/Documents/stage/CDD_pyHCA/CDD/analyse/test_couverture/'
    liste_fasta = os.listdir(directory)

    hca_out = '/home/issa/Documents/stage/CDD_pyHCA/CDD/analyse/test_couverture/'
    cdd_out = '/home/issa/Documents/stage/CDD_pyHCA/CDD/analyse/test_couverture/'



    resid_HCA, domain_HCA, resid_CDD, domain_CDD = taux_couverture (directory, liste_fasta, hca_out, cdd_out)
    plotting(resid_HCA, domain_HCA, resid_CDD, domain_CDD)
    print('succefully done !!!')
