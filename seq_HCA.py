#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This programm is for matching the orphan domains (domains that found in HCA and not in CDD)
"""
import os, re
from collections import defaultdict
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
def read_fasta (fasta_in):
    """return a dict which key = protein_id and value = list of amino acid in the sequence
    """
    number_protein = 0
    with open(fasta_in, 'r') as fasta:
        dico_fasta = {}
        for line in fasta:
            line = line.strip()
            if line.startswith('>'):
                number_protein+=1
                sequence = ''
                header = line[1:].split()[0]
                dico_fasta[header] = []
            else:
                sequence += line.strip()
                for aa in sequence:
                    dico_fasta[header].append(aa)
    return dico_fasta, number_protein


def read_pyHCA_outF (hca_in):
    """return domains match in sequences by pyHCA program with lenght >= 30 aa
    dict : key = protein_id , value (list of tuples) = tpl(domain_start, domain_end, domain_lenght)
    """

    domain_hca = defaultdict(set)
    with open(hca_in, 'r') as hca:
        for line in hca:
            line = line.strip()
            domain = ()
            if line.startswith('>'):
                pos = line[1:].split()
                header = pos[0]
            if line.startswith('domain'):
                domain_start  = line.strip().split()[1]
                domain_start  = int(domain_start) - 1
                domain_end    = int(line.strip().split()[2])
                domain_lenght = domain_end - domain_start
                if domain_lenght >= 30:
                    domain   = (domain_start, domain_end, domain_lenght)
                    domain_hca[header].add(domain)

    return domain_hca

def read_cdd_outF (proteome_name, cdd_in):
    """return domain match in sequences by CDD scan
    dict : key = protein_id , value(list of tuples) = tpl(domain_start, domain_end, domain_lenght)
    """
    domain_cdd = defaultdict(set)
    with open(cdd_in, 'r') as cdd:
        for line in cdd:
            line = line.strip()
            domain = ()
            if line.startswith('Q#'):
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

    return domain_cdd

def binary_position(dico_fasta, domain_cdd):

    dico_bin = {}
    for protein in dico_fasta:
        fasta_seq = dico_fasta[protein]
        lenght_seq = len(fasta_seq)

        if protein in domain_cdd:
            liste_of_position = [0]*lenght_seq
            cdd_dom = domain_cdd[protein]
            for start, stop, lenght in cdd_dom:
                for i in range (int(start), int(stop)):
                    liste_of_position[i] = 1
                    tot = sum(liste_of_position)
            dico_bin[protein] = liste_of_position
    return dico_bin

def match_domain(proteome_name, dico_fasta, domain_hca, domain_cdd, dico_bin, number_protein):

    """find the orphan domains between HCA and CDD domain_list with coverage <= 20%
    """
    orphan_domains = defaultdict(list)
    prot_has_both_dom = set(domain_hca.keys()).intersection(set(domain_cdd.keys()))
    for protein in domain_hca:
        if protein in domain_cdd:
            hca_pos_domain = domain_hca[protein]
            list_position_bin = dico_bin[protein]
            for start, stop, lenght in hca_pos_domain:
                couverture_hca = []
                dom_hca = (start, stop, lenght)
                for i in range(int(start), int(stop)):
                    couverture_hca.append(list_position_bin[i])
                couverture = sum(couverture_hca)/lenght
                if couverture <= 0.2:
                    orphelin = (start, stop,couverture)
                    orphan_domains[protein].append(orphelin)
        else:
            protein_orph = domain_hca[protein]
            for start, stop, lenght in protein_orph:
                false_lengh = 0.00000
                orphelin = (start, stop, false_lengh)
                orphan_domains[protein].append(orphelin)

    return orphan_domains

def write_outFile (orphan_domains, outF, dico_fasta): #dico_fasta

    # Write positions of orphan domains in out files
    # with open (outF, 'w') as outfile:
    #     for protein in orphan_domains:
    #         pos = orphan_domains[protein]
    #         for start, stop, lenght in pos:
    #             outfile.write('>'+protein+'\n'+'orphan domain'+'\t'+str(start)+'\t'+str(stop)+'\t'+str(lenght)+'\n')

    # Write orphan domains in fasta format
    with open (outF, 'w') as outfile:
        for protein in dico_fasta:
            sequence = dico_fasta[protein]
            if protein in orphan_domains:
                pos = orphan_domains[protein]
                for start, stop, lenght in pos:
                    start = int(start)
                    stop = int (stop)
                    domain = (sequence[start:stop])
                    sequence_domain = ''.join(str(aa) for aa in domain)
                    print(sequence_domain)
                    outfile.write('>'+protein+'\n'+str(sequence_domain)+'\n')
def main():

    list_of_count = []

    for filename in liste_fasta:
        count = ()
        if filename.endswith('.fasta'):
            fasta_in = os.path.join(directory, filename)
            hca_in = os.path.join(directory, filename+".hca")
            cdd_in = os.path.join(directory, filename+".cdd")
            proteome_name = filename.split('.')[0]
            outF = os.path.join(dir_out, filename+".orp")

            fasta, protein_number = read_fasta(fasta_in)
            hca = read_pyHCA_outF(hca_in)
            cdd = read_cdd_outF(proteome_name, cdd_in)
            positions_dom = binary_position(fasta, cdd)
            orph = match_domain (proteome_name, fasta, hca, cdd, positions_dom, protein_number)
            out_file = write_outFile(orph, outF, fasta)


            count = (len(orph), proteome_name, protein_number)
            list_of_count.append(count)

    list_of_count.sort(reverse=True)


    orphan_dom_count, proteome_labels, protein_count  = [], [], []

    protein_pour, orphelin_count_all = [], []
    for orphelin, proteome, protein_num in list_of_count:

        protein = (protein_num/protein_num)*100
        domain_pourcentage = (orphelin/protein_num)*100
        protein_pour.append(protein)
        orphelin_count_all.append(domain_pourcentage)

        orphan_dom_count.append(orphelin)
        proteome_labels.append(proteome)
        protein_count.append(protein_num)

    fig, ax = plt.subplots()
    N = len(proteome_labels)
    ind = np.arange(N)
    width = 0.20

    plot_line = ax.bar(ind, protein_pour, width=0.2, alpha=1, color='w', label= 'Nombre de proteines total dans proteome')
    plot_lin = ax.bar(ind, orphelin_count_all, width=0.2, alpha=0.5, color='r', label= 'Nombre de proteines avec des domaines orphelins')

    for a,b in zip(ind, protein_count):
        plt.text(a, protein_pour[a], str(b), va = 'bottom', fontsize=10, fontdict={'family': 'serif', 'color':  'k', 'weight': 'normal','size': 12})
    for a,b in zip(ind, orphan_dom_count):
        plt.text(a+0.1, orphelin_count_all[a], str(b), va = 'baseline', rotation = 'vertical', fontsize=12, fontdict={'family': 'serif', 'color':  'k', 'weight': 'normal','size': 12})


    ax.set_xlim(-width,len(ind)+width)
    ax.set_xticks(ind+width)
    ax.set_xticklabels (proteome_labels, rotation='vertical', fontsize=10)
    for name in ax.set_xticklabels(proteome_labels, rotation= 'vertical', fontsize=10):
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
    plt.title(u"Taux de couverture en domaines orphelins\n {Domaines presents sur HCA et absents sur CDD}", fontsize=17, fontdict={'family': 'monospace'})
    ax.legend(loc='upper left', fontsize=10)
    plt.show()
if __name__ == '__main__':
    print('running...')
    # directory = '/home/issa/Documents/stage/CDD_pyHCA/data/'
    directory = '/home/issa/Documents/stage/CDD_pyHCA/CDD/analyse/test_couverture/'
    # dir_out   = '/home/issa/Documents/stage/CDD_pyHCA/data/orphans_domaines_pos/'
    dir_out   = '/home/issa/Documents/stage/CDD_pyHCA/CDD/analyse/orphans_domaines_faa/'
    liste_fasta = os.listdir(directory)
    main()
