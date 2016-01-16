#!/usr/bin/env python
# -*- coding: utf-8 -*-
''' Distribution de la Taille des séquences dans les Protéomes
Pour chaque proteome, construire un plot des tailles de sequences 
'''

import sys
from pylab import *
import matplotlib
import matplotlib.pyplot as plt
import os


def plot_P20 (proteome_dir):
    list_taille = []
    for filename in os.listdir(proteome_dir):
        if filename.endswith(".fasta"):
            path_proteome = os.path.join(proteome_dir, filename)
            seqLen = 0
            list_seq_lenght= []
            list_seq_name = []
            Taille_proteome = 0
            f = open (path_proteome, 'r')
            for line in f:
                if line != '':
                    line = line[:-1]
                    if line [0] != '>':
                        Taille_proteome= Taille_proteome+1
                        seqLen = len(line)
                        list_seq_lenght.append(seqLen)
            list_taille.append(list_seq_lenght)
            f.close()
            
    #Construire un plot de liste des tailles de sequences pour chaque proteomes
    plt.figure()
    plt.boxplot(list_taille,  0, 'gD')
    plt.grid(True)
    plt.ylabel('Taille des Séquences')
    plt.xlabel('Protéomes')
    plt.title(u"Distribution de la Taille des séquences dans les Protéomes")
    plt.savefig(proteome_dir+u'boxplotProteome.png')
    plt.show()

# proteome_dir = '/home/issa/Documents/STAGE/Data/Init_data/Proteomes_test/ex/'
proteome_dir = '/home/issa/Documents/STAGE/Data/Init_data/Proteomes_test/'
plot_P20(proteome_dir)
