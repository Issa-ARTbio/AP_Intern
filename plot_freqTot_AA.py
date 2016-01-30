#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''plot sur une meme figure la frequence de chaque aa contenue dans les proteomes
-i(file.dat)


usage:
=====
$ python3 plot_freqTot_AA.py -i(.dat)
'''

import os, sys, argparse, re
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import itertools

from collections import OrderedDict

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
            # print(proteome_name)
            values = tmp[1:]
            for i, val in enumerate (values):
                dico_freq.setdefault (liste_aa[i], dict())
                dico_freq[liste_aa[i]][proteome_name] = (float(val.split()[1]))/100

    N = 111
    n = 20
    ind = np.arange(N)
    indp= np.arange(n)
    datas= []
    d=[]
    # fig = plt.figure()
    width = 0.35
    labels_aa= []
    labels_prt = []
    my_handle = []
    for aa in dico_freq.keys():
        labels_aa.append(aa)
        data = dico_freq[aa]
        datas.append(data)
        values_aa = []
        # linestyles = {'1': 'solid', '2': 'dashed', '3': 'dashdot', '4': 'dotted'}
        # markers = {'d': 'thin_diamond', '*': 'star', '^': 'triangle_up', '+': 'plus', 's': 'square'
        # color = {}
        for proteome in data:
            # mk = markers[m]
            data_aa=(data[proteome])
            values_aa.append(data_aa)
            # print(ind, values_aa)
            d.append(values_aa)
            values_aa = sorted(values_aa)
            labels_prt.append(proteome)
        zipper = list(zip(labels_prt, values_aa))
        zipper_sort = sorted (zipper)
        key_trie, val_trie = zip(*zipper_sort)

        petit = ['S', 'T','C','P']
        apolaire = ['A', 'G' , 'V', 'L', 'M', 'I', 'F', 'Y', 'W']
        polaire = ['D', 'E', 'K', 'R', 'H', 'N', 'Q']
        # filled_markers = ('o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd')

        if aa in apolaire:
            my_handle.append(aa)
            if aa == 'A': plt.plot(ind, values_aa, linewidth = 2 , color = 'purple',marker= '>', label= 'apolaire (A)')#, zorder=1)
            elif aa == 'G': plt.plot(ind, values_aa, linewidth = 2 , color = 'purple',marker= '*', label= 'apolaire (G)')#, zorder=1)
            elif aa == 'V': plt.plot(ind, values_aa, linewidth = 2 , color = 'purple',marker= 'v', label= 'Aliphatique (V)')#, zorder=1)
            elif aa == 'L': plt.plot(ind, values_aa, linewidth = 2 , color = 'purple',marker= 'd', label= 'Aliphatique (L)')#, zorder=1)
            elif aa == 'M': plt.plot(ind, values_aa, linewidth = 2 , color = 'purple',marker= 'p', label= 'Hydrophobe fort(M)')#, zorder=1)
            elif aa == 'I': plt.plot(ind, values_aa, linewidth = 2 , color = 'purple',marker= '1', label= 'Aliphatique (I)')#, zorder=1)
            elif aa == 'F': plt.plot(ind, values_aa, linewidth = 2 , color = 'purple',marker= '^', label= 'aromatique (F)')#, zorder=1)
            elif aa == 'Y': plt.plot(ind, values_aa, linewidth = 2 , color = 'purple',marker= 's', label= 'aromatique (Y)')#, zorder=1)
            elif aa == 'W': plt.plot(ind, values_aa, linewidth = 2 , color = 'purple',marker= 'h', label= 'aromatique (W)')#, zorder=1)
        elif aa in polaire:
            my_handle.append(aa)
            if aa == 'K': plt.plot(ind, values_aa, linewidth = 2 , color = 'green',marker= 'o', label= 'Polaire charge positif (K)')#, zorder=2)
            elif aa == 'R': plt.plot(ind, values_aa, linewidth = 2 , color = 'green',marker= '*', label= 'Polaire charge positif (R)')#, zorder=2)
            elif aa == 'H': plt.plot(ind, values_aa, linewidth = 2 , color = 'green',marker= 'v', label= 'Polaire charge positif (H)')#, zorder=2)
            elif aa == 'D': plt.plot(ind, values_aa, linewidth = 2 , color = 'green',marker= '>', label= 'Polaire charge negatif (D)')#, zorder=2)
            elif aa == 'E': plt.plot(ind, values_aa, linewidth = 2 , color = 'green',marker= '<', label= 'Polaire charge negatif (E)')#, zorder=2)
            elif aa == 'N': plt.plot(ind, values_aa, linewidth = 2 , color = 'green',marker= 's', label= 'polaire non charge (N)')#, zorder=2)
            elif aa == 'Q': plt.plot(ind, values_aa, linewidth = 2 , color = 'green',marker= 'd', label= 'polaire non charge (Q)')#, zorder=2)
        elif aa in petit:
            my_handle.append(aa)
            if aa == 'S': plt.plot(ind, values_aa, linewidth = 2 , color = 'r',marker= 'o', label= 'polaire non charge (S)', zorder=3)
            elif aa == 'T': plt.plot(ind, values_aa, linewidth = 2 , color = 'r',marker= '*', label= 'polaire non charge (T)', zorder=2)
            elif aa == 'C': plt.plot(ind, values_aa, linewidth = 2 , color = 'r',marker= 'v', label= 'polaire non charge (C)', zorder=4)
            elif aa == 'P': plt.plot(ind, values_aa, linewidth = 2 , color = 'r',marker= '|', label= 'apolaire (P)', zorder=1)
            # plt.legend (ncol = 1)
    # print (my_handle)

    f = (ind+width)
    f = sorted (f)
    plt.grid(True)
    plt.xticks(f , fontsize=2)
    plt.xlabel('Proteomes', fontsize=15, color='blue')
    plt.ylabel('Frequence', fontsize=15, color='red')
    plt.title('Distribution de la Frequence des aa dans les Proteomes', fontdict={'family': 'monospace'})
    ncol = 3
    for aa in liste_aa:
        if aa in apolaire:
            plt.legend(ncol = [i for i in range(ncol)])
        if aa in apolaire:
            plt.legend(ncol = [i+1 for i in range(ncol)])
        if aa in apolaire:
            plt.legend(ncol = [i+2 for i in range(ncol)])
    # plt.legend(ncol = 1, loc='upper left', fontsize=8, title = 'apolaire')
    # l = plt.legend(ncol=3, loc='upper left', fontsize=8)
    # itertools.chain(*[aa[i::ncol] for i in range(ncol)])
    plt.show()
    # plt.savefig(aa+u'_Freqence.pdf')


proteom_composition_aa = './old_/proteome_composition_AA.csv'
frequence_aa = Freq_AA (proteom_composition_aa)
