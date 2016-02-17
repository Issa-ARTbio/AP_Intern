#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''Extraire les sequences fasta correspondant à partir des noms de proteins

Parametres

-in file clusters file fasta all -o fasta file for each cluster

usage
======
python sequence_extraction.py
'''

import sys, os

def read_All_fasta(directory, liste_fasta):
    '''extract fasta sequence corresponding to the names in core file
    '''
    dict_proteome = {}
    list_protein_name = []
    for filename in liste_fasta:
        if filename.endswith('.fasta'):
            path_proteome = os.path.join(directory, filename)
            proteome_name = filename.split('.')[0]
            for line in open(path_proteome):
                line = line.strip()
                if line.startswith('>'):
                    name_protein = line[1:].split()[0]
                    # print(name_protein)
                    if name_protein in dict_proteome:
                        print("Error, protein {} already exists".format(name_protein))
                        sys.exit(1)
                    dict_proteome[name_protein] = ""
                else:
                    dict_proteome[name_protein] += line.strip()

    return dict_proteome

def extract_sequence(dict_proteome,core, cluster):

    nb = 0
    for line in open(core):
        nb+=1
        line = line.strip()
        tmp  = line.split()
        with open(cluster+str(nb)+'.fasta' , 'w' ) as outfile:
            for indice in range(len(tmp)):
                protein = tmp[indice]
                for name_protein in dict_proteome:
                    if name_protein==protein:
                        print(name_protein)
                        protein_name = '>'+str(protein)
                        outfile.write(protein_name+'\n')
                        outfile.write(dict_proteome[name_protein]+'\n')


if __name__ == '__main__':
    core= '/home/issa/Documents/stage/orthomcl/orthomcl_results/new_groups.txt'
    # all_fasta='/home/issa/Documents/stage/init_data/CyanobacteriaProteins.fas'
    cluster= '/home/issa/Documents/stage/orthomcl/orthomcl_results/core_cluster/cluster_'
    directory = '/home/issa/Documents/stage/orthomcl/proteomes_format/'
    liste_fasta = os.listdir(directory)
    all_fasta = read_All_fasta(directory, liste_fasta)
    extract_sequence (all_fasta, core, cluster)
