#!/usr/bin/env python
# -*- coding: utf-8 -*-

''' Renvoie la Distribution d'un acide aminé donné pour chaque protéome
-i directory(.fasta) -o plot (.png)
'''
import sys
from pylab import *
import matplotlib
import matplotlib.pyplot as plt
import os

def Freq_AA (proteome_dir, list_aa):
    #Creation d'un dico vide pour stocker la fréquence de chaque aa
    dict_aa = {}
    list_filename = []
    for filename in os.listdir(proteome_dir):
        if filename.endswith(".fasta"):
            list_filename.append(filename)
            path_proteome = os.path.join(proteome_dir, filename)
            Taille_proteome = 0 # Initialise la taille de chaque protéome
            with open(path_proteome) as f:
                Prot = {}
                for line in f:
                    line = line[:-1]
                    if line [0] != '>':
                        for aa in line:
                            Taille_proteome= Taille_proteome+1
                            Prot[aa] = Prot.get(aa, 0) +1
                for aa in list_aa:
                    freq = Prot.get(aa, 0) / Taille_proteome
                    if aa not in dict_aa:
                        dict_aa[aa] =  []
                    dict_aa[aa].append(freq)
    # print (dict_aa)
    N = 10
    ind = np.arange(N)
    width = 0.20
    for aa in list_aa:

        fig, ax = plt.subplots()
        ax.plot(ind, dict_aa[aa], color='g', label='Issa',linewidth=3)
        # bar_Num = ax.bar(dict_aa[aa], color='g')
        ax.set_xticks(ind+width)
        ax.set_xticklabels (list_filename, rotation=70)
        plt.grid(True)
        ax.set_ylabel(u'Taille des Séquences')
        ax.set_xlabel(u'Protéomes')
        ax.set_title(u"Distribution de la Fréquences de "+aa+" dans les Protéomes")
        save_file = plt.savefig(proteome_dir+u'_Freq_aa.png')
        # plt.show()
        plt.close()
    return (save_file)
# proteome_dir = '/home/issa/Documents/STAGE/Data/Init_data/Proteomes_test/ex/'
proteome_dir = '/home/issa/Documents/STAGE/Data/Init_data/Proteomes_test/'
list_aa = ['A', 'R', 'N', 'D', 'C', 'Q', 'E','G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y','V']
print ('Loading Plots...')
Freq_AA(proteome_dir, list_aa)
print ('Finished')
