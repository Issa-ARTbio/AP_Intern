#!/usr/bin/env python
# -*- coding: utf-8 -*-

#Suite du program find_OD_in_protein

'''find the OD positions in protein sequences
'''

import os, sys
from collections import defaultdict

def read_protein_id (list_id_OD):

    dict_protein_id = defaultdict(list)
    with open(list_id_OD) as f:
        for line in f:
            line = line.strip()
            element = line.split(',')
            protein = element[0][2:-1]
            long_protein = element[2]
            cluster_OD = element[-1][2:-2]
            tpl = (protein, long_protein)
            dict_protein_id[cluster_OD].append(tpl)

    return dict_protein_id

def read_fasta(directory, list_fasta, dict_protein_id):


    nb=0
    for filename in list_fasta:
        if filename.endswith('.fasta'):
            fasta_in = os.path.join(directory, filename)
            cluster_id = filename.split('.')[0]
            with open(fasta_in, 'r') as fasta:
                dico_fasta = {}
                for line in fasta:
                    line = line.strip()
                    if line.startswith('>'):
                        header = line[1:].split()[0]
                        protein_id = header.split('|')[1]
                        dico_fasta[protein_id] = 0
                    else:
                        dico_fasta[protein_id] += len(line)

                for cluster in dict_protein_id:
                    list_pr = dict_protein_id[cluster]
                    for protein, aa in list_pr:
                        if protein in dico_fasta:
                            nb+=1
                            print(cluster, protein, aa, cluster_id)
                            print('***************************')
    print(nb)
    return dico_fasta
if __name__ == '__main__':

    list_id_OD = '/home/issa/Documents/stage/cd-hit/biominerales_clusters/strout_OD_matches.txt'
    id_OD = read_protein_id (list_id_OD)
    directory = '/home/issa/Documents/stage/cd-hit/biominerales_clusters/cluster_orthoMCL_found/'
    list_fasta = os.listdir(directory)
    fasta = read_fasta(directory, list_fasta, id_OD)
