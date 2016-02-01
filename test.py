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
        petit = ['S', 'T','C','P']
        filled_markers = itertools.cycle(('o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd'))

        my_label_dico = {A:'apolaire (A)', G:'apolaire (G)', V:'Aliphatique (V)' , L:'Aliphatique (L)', M:'Hydrophobe fort(M)', I:'Aliphatique (I)', F:'aromatique (F)', Y:'aromatique (Y)', W:'aromatique (W)',
        K:'Polaire charge positif (K)', R:'Polaire charge positif (R)', H:'Polaire charge positif (H)', D:'Polaire charge negatif (D)', E:'Polaire charge negatif (E)', H:'polaire non charge (N)', Q:'polaire non charge (Q)',
        S:'polaire non charge (S)', T:'polaire non charge (T)', C:'polaire non charge (C)', P:'apolaire (P)'}
        for aa in my_label_dico:
            if aa in polaire:
                print(aa)
                plt.plot(ind, values_aa, linewidth = 2 , color = 'purple',marker= filled_markers.next(), label= my_label_dico[aa])


        f = (ind+width)
        f = sorted (f)
        plt.grid(True)
        plt.xticks(f , fontsize=2)
        plt.xlabel('Proteomes', fontsize=15, color='blue')
        plt.ylabel('Frequence', fontsize=15, color='red')
        plt.title('Distribution de la Frequence des aa dans les Proteomes', fontdict={'family': 'monospace'})
        plt.legend(ncol = 3, loc='upper left', fontsize=8, title = 'apolaire')
        plt.show()
        # plt.savefig(aa+u'_Freqence.pdf')


proteom_composition_aa = '/home/issa/Documents/STAGE/Results/proteome_composition_AA.csv'
frequence_aa = Freq_AA (proteom_composition_aa)
