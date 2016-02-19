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
        if cluster.endswith('.fasta'):
            path_cluster = os.path.join(directory, cluster)
            cluster = cluster.split('.')[0]
            f = SeqIO.parse(open(path_cluster),'fasta')
            for line in f:
                name, sequence = line.id, str(line)
                tmp  = name.split('|')
                name_protein       = tmp[1]
                full_proteome_name = tmp[0]
                # line = ('>',name, '\n' , sequence, '\n')
                if full_proteome_name not in dico_all:
                    dico_all[full_proteome_name] = sequence
                else:
                    dico_all[full_proteome_name] +=sequence
    print(len(dico_all))
    return dico_all

def write_seq (dico_all, core_genome):


    with open(core_genome, 'w') as outf:
        for proteome_id in dico_all:
            outf.write('>'+proteome_id+'\n'+dico_all[proteome_id]+'\n')



rep = '/home/issa/Documents/stage/muscle/core_to_align_fullname/'
liste_cl = os.listdir(rep)
core_genome = 'clusters_fusionned.fasta'
dico = read_fasta_aln(rep, liste_cl)
write_seq (dico, core_genome)
