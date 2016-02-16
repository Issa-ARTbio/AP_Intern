#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''Extraire les sequences fasta correspondant à partir des noms de proteins

Parametres

-in file clusters file fasta all -o fasta file for each cluster

usage
======
python sequence_extraction.py
'''
def extract_seq(core, all_fasta, cluster):
    '''extract fasta sequence corresponding to the names in core file
    '''
    dict_proteome = {}
    list_protein_name = []
    with open(all_fasta) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                name_protein = line[1:].split()[0]
                # list_protein_name.append(name_protein)
                if name_protein in dict_proteome:
                    print("Error, protein {} already exists".format(name_protein))
                    sys.exit(1)
                dict_proteome[name_protein] = ""
            else:
                dict_proteome[name_protein] += line.strip()

    nb = 0
    for line in open(core, 'r'):
        nb+=1
        line = line.rstrip()
        tmp  = line.split(',')
        with open(cluster_+str(nb)+'.fasta' , 'w' ) as outfile:
            for indice in range(len(tmp)):
                protein = tmp[indice]
                for name_protein in dict_proteome:
                if name_protein==protein:
                    protein_name = '>'+str(protein)
                    outfile.write(protein_name+'\n'+dict_proteome[name_protein])


if __name__ == '__main__':
    core= 'test2.txt'
    all_fasta='test1.txt'
    cluster= 'cluster_'
    extract_seq(core, all_fasta, cluster_)
