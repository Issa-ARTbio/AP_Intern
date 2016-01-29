#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''Description des données
Pour chaque fichier fasta: renvoie le nom et la longeur de chaque séquence
Se_description.py -i .fasta -o .txt Protein:__ Taille:__'''

import sys
from pylab import *
import matplotlib
import matplotlib.pyplot as plt
import os

def decrip_prot (proteome_directory, out_file) :
    ''
    ''
    nbseq = 0
    seqLen = 0
    list_seq_lenght= []
    list_seq_name = []
    # aa = "'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'O', 'U', 'S', 'T', 'W', 'Y', 'V'"
    f = open (proteome_directory, 'r')
    # s = open (out_file, 'w')
    for line in f:
        line = line[:-1]
        if line[0] == '>':
            nbseq = nbseq + 1
            seq_name = line [1:]
            list_seq_name.append(seq_name)
        if line [0] != '>':
            seqLen = len(line)
            list_seq_lenght.append (seqLen)
            # s.write ( 'Potein:'  + (seq_name) + '\t' + 'Taille:'  + str(seqLen) + '\n' )
    print (list_seq_lenght)

    # s.close()
    f.close()
    # Construire l'histogramme
    plt.hist(list_seq_lenght,  50, facecolor='green')
    #labels
    plt.xlabel(u'Nombre de Proteines')
    plt.ylabel(u'Taille de Sequence')
    plt.title(u'Distribution de la Taille des séquences')
    # plt.axis ([0 , 2000])
    plt.grid(True)
    plt.savefig('/home/issa/Documents/STAGE/Data/Données_Initiales/Proteomes_test/'+proteome_directory+'.pdf')
    # plt.show()

# proteome_directory = '/home/issa/Documents/STAGE/Analyses_Preliminaires/Script_Analyses/example.fasta'
proteome_directory = '/home/issa/Documents/STAGE/Data/Données_Initiales/Proteomes_test/'
out = '/home/issa/Documents/STAGE/Data/Données_Initiales/Proteomes_test/'
# print ( 'Proteome\t'+"\t".join(list_aa) )
for filename in os.listdir(proteome_directory):
    if filename.endswith(".fasta"):
        path_proteome = os.path.join(proteome_directory, filename)
        decrip_prot(path_proteome, out)
