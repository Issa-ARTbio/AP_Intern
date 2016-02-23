#!/usr/bin/env python
# -*- coding: utf-8 -*-


'''
'''

import os
from Bio import SeqIO


def read_fasta_aln (directory, liste_cluster):
    '''extract fasta sequences from dir, convert sequence in one string
    Return a dico with proteome_id as key and all corresponding sequence in files as value
    '''
    dico_all = {}
    for cluster in liste_cluster:
        if cluster.endswith('.aln'):
            path_cluster = os.path.join(directory, cluster)
            cluster = cluster.split('.')[0]
            with open(path_cluster) as f:
                for line in f:
                    if line.startswith('>'):
                        line = line.strip()
                        name = line[1:]
                        tmp  = name.split('|')
                        name_protein       = tmp[1]
                        full_proteome_name = tmp[0]
                    else:
                        sequence = line.strip()

                        if full_proteome_name not in dico_all:
                            dico_all[full_proteome_name] = sequence
                        else:
                            dico_all[full_proteome_name] +=sequence

    print(len(dico_all))
    print (len(dico_all[full_proteome_name]))
    return dico_all

def write_seq (dico_all, core_genome):
    '''write in a output file the fusionned sequence
    '''

    with open(core_genome, 'w') as outf:
        for proteome_id in dico_all:
            outf.write('>'+proteome_id+'\n'+dico_all[proteome_id]+'\n')



rep = '/home/issa/Documents/stage/muscle/final_alignement/'
liste_cl = os.listdir(rep)
core_genome = 'clusters_fusionned.fasta'
dico = read_fasta_aln(rep, liste_cl)
write_seq (dico, core_genome)
