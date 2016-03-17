#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''plot de la frequence de chaque aa contenue dans les proteomes


Parameters
=============
-i(file.dat) -o (plot.pdf pour chaque acide amine)


Usage:
=====
$ python3 frequence_AA.py -i(.dat)
'''

import os, sys, argparse, re
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec


def Freq_AA (proteom_composition_aa):
    """ 'Description de la fonction' : A partir du fichier Taille_Proteome.dat contenant la composition en acide amine pour chaque proteome
    renvoit un plot de la frequence de chaque acide amine sur les 110 proteomes

    Parameters:
    ===========
    proteom_composition_aa: file
        un fichier avec les 20 acides amines en colonnes separes par '\ŧ' en premiere colonne = nom des proteomes
        le nombre et le pourcentage de chaque aa en ligne separe par un space

    outf: file for each aa
        plot de la frequence de chaque aa dans les differents proteomes """

    dico_freq = dict()
    with open(proteom_composition_aa) as f:
        line = f.readline()
        line = line.strip()
        liste_aa = line.split('\t')[1:]
        for line in f:
            line = line.strip()
            tmp = line.split('\t')
            proteome_name = tmp[0]
            proteome_name = proteome_name.split('.')[0]
            values = tmp[1:]
            for i, val in enumerate (values):
                dico_freq.setdefault (liste_aa[i], dict())
                dico_freq[liste_aa[i]][proteome_name] = (float(val.split()[1]))/100

    N = 111
    ind = np.arange(N)
    width = 0.35

    for aa in dico_freq.keys():
        data = dico_freq[aa]
        values_aa = []
        for proteome in data.keys():
            data_aa= data[proteome]
            values_aa.append(data_aa)

        zipper = list(zip( values_aa, data.keys()))
        zipper_sort = sorted(zipper)
        val_trie, key_trie = zip(*zipper_sort)
        fig, ax = plt.subplots()
        ax.plot(ind, val_trie, color='black',linewidth=3)
        ax.set_xticks(ind)
        ax.set_xticklabels (key_trie, rotation= 'vertical', fontsize=8)
        for i in ax.set_xticklabels(key_trie, rotation= 'vertical', fontsize=8):
            if re.search('Synechococcus_sp_PCC_6312', str(i)) :
                i.set_color('green')
            elif re.search('Synechococcus_calcipolaris' , str(i)):
                i.set_color('green')
            elif re.search('Thermosynechococcus_elongatus_BP1' , str(i)):
                i.set_color('green')


            elif re.search('Gloeomargarita_lithophora' , str(i)):
                i.set_color('red')
            elif  re.search ('Cyanothece_sp_PCC_7425', str(i)):
                i.set_color('red')
            elif  re.search ('Chroococcidiopsis_thermalis_PCC_7203', str(i)):
                i.set_color('red')
        plt.grid(True)
        ax.set_ylabel(u'Frequence', fontsize=15, color='blue')
        ax.set_ylim (0, 0.15)
        ax.set_xlabel(u'Proteomes', fontsize=15, color='blue')
        ax.set_title(u"Distribution de la Frequence de "+aa+" dans les Proteomes", fontsize=17, fontdict={'family': 'monospace'})

        plt.show()
        # plt.savefig(aa+u'_Freqence.pdf')


proteom_composition_aa = '/home/issa/Documents/STAGE/Firsts_Analyses/Resultats/proteome_composition_AA.csv'
frequence_aa = Freq_AA (proteom_composition_aa)
