#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Renvoit deux polar_plots pour les frequences borders et les frequence bio
"""

import os, sys, argparse, re
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import random

def read_freq (proteom_composition_aa):
    """ 'Description de la fonction' : A partir du fichier Taille_Proteome.dat contenant la composition en acide amine pour chaque proteome
    renvoit un dict{ aa, {data_aa}}

    """

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

    return (dico_freq)

def get_all_freq(proteome_name, dico_freq, list_aa):
    """ return the list of AA frequencies for a particular proteome
    """
    list_freq = []

    for aa in list_aa:
        fr_aa = dico_freq[aa][proteome_name]
        list_freq.append(fr_aa)

    pass



def get_name_min_med_max (dico_freq, amino, list_aa):
    """renvoit le nom du proteome avec la frequence minimale de chaque aa dans amino
    """
    freq_All_min, freq_All_med, freq_All_max = [], [], []

    my_dico_min = {}
    my_dico_med = {}
    my_dico_max = {}

    for aa in amino:

        my_dico_min[aa] = []
        my_dico_med[aa] = []
        my_dico_max[aa] = []
        # name_min_proteome, name_med_proteome, name_max_proteome = None, None, None
        data_aa = dico_freq[aa]
        proteome_name, values_aa = zip(*data_aa.items())
        zipper = list(zip(values_aa, proteome_name))
        zipper_sort = sorted(zipper)
        val_trie, key_trie = zip(*zipper_sort)
        # print (key_trie , val_trie)
        medium = len(val_trie)/2 + 0.5
        name_min_proteome = key_trie [0]
        my_dico_min[aa].append(name_min_proteome)

        name_med_proteome = key_trie[int(medium)]
        my_dico_med[aa].append(name_med_proteome)

        name_max_proteome = key_trie [-1]
        my_dico_max[aa].append(name_max_proteome)

        for aa2 in list_aa:

            if name_min_proteome in data_aa:
                fr = dico_freq[aa2][name_min_proteome]
                my_dico_min[aa].append(fr)

            if name_med_proteome in data_aa:
                fr = dico_freq[aa2][name_med_proteome]
                my_dico_med[aa].append(fr)

            if name_med_proteome in data_aa:
                fr = dico_freq[aa2][name_max_proteome]
                my_dico_max[aa].append(fr)
    # print ('my_dico_min', my_dico_min, 'my_dico_med',my_dico_med,  'my_dico_med', my_dico_med)

    return my_dico_min, my_dico_med,  my_dico_max


def polar_plotting (my_dico_min, my_dico_med,  my_dico_max, list_aa):

    # fig = super(PlotWindPowerDensity, self).get_figure()
    for aa in my_dico_min:
        # if aa == 'L':
        data_min,  data_med, data_max = None, None, None
        liste_min, liste_med, liste_max = [], [], []

        data_min = my_dico_min[aa]
        liste_min.append(data_min)
        name_min_proteome = liste_min[0]


        data_med = my_dico_med[aa]
        liste_med.append(data_med)
        name_med_proteome = liste_med[0]

        data_max = my_dico_max[aa]
        liste_max.append(data_max)
        name_max_proteome = liste_max[0]

        ax = plt.subplot(111, polar=True)



        # print (data_min_s, list_aa)
        r = np.arange(len(list_aa))
        rad = (r*2*np.pi/len(r))
        rad = list(rad)
        rad.append(rad[0])
        # ax.plot(data,  color='r',linewidth=1,alpha= 0.5, label = data_min[0], marker= '*')
        # N = 20
        # theta = np.linspace(0.0, 2 * np.pi, N, endpoint=False)
        # fig = ax.bar(r, data_min_s,  color='r',linewidth=5,alpha= 0.5, label = data_min[0], )
        # i = float([0])
        data_min_s = list(data_min[1:])
        data_min_s.append(data_min_s[0])

        data_med_s = list(data_med[1:])
        data_med_s.append(data_med_s[0])

        data_max_s = list(data_max[1:])
        data_max_s.append(data_max_s[0])

        ax.plot(rad, data_min_s,  color='r',linewidth=5,alpha= 0.3, label = data_min[0]+'\n Minimale', marker= '*')
        ax.fill (rad, data_min_s, color= 'r', alpha= 0.5)

        ax.plot(rad, data_med_s,  color='b',linewidth=5,alpha= 0.3, label = data_med[0]+'\n Mediane', marker= 's')
        ax.fill (rad, data_med_s, color= 'b', alpha= 0.5)

        ax.plot(rad, data_max_s,  color='g',linewidth=5,alpha= 0.7, label = data_max[0]+'\n Maximale', marker= '>')
        ax.fill (rad, data_max_s,color= 'g', alpha= 0.5)

        plt.grid(True)
        rad = rad[:-1]
        ax.set_xticks(rad)
        ax.set_rmax(0.15)
        # ax.set_xticklabels (list_aa) #, color= 'b') #for i in ['A', 'K', 'N', 'L'] )
        hydrophobe = (['V', 'I', 'L', 'M', 'F', 'W', 'Y'])
        polaire = (['K', 'R', 'H', 'D', 'E', 'N', 'S', 'Q',  'T'])
        petit_apolaire = (['P', 'G', 'A', 'C'])
        list_lab = (hydrophobe, polaire, petit_apolaire)
        colors = [ 'b', 'm', 'grey']
        for xtick, i in zip (ax.get_xticklabels(), list_aa):

            if i in hydrophobe:
                xtick.set_color(colors[0])
            if i in polaire:
                xtick.set_color(colors[1])
            if i in petit_apolaire:
                xtick.set_color(colors[2])
                ax.set_xticklabels(list_aa)
        ax.set_title('Frequence des differents aa dans les proteomes\n avec une grande variation de la frequence de '+aa, va ='bottom', color = 'k', fontdict={'family': 'monospace'})
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.10),
  fancybox=True, shadow=True, ncol=3)
        # , label= 'b' if i == rad [0] else '')
        # fig.add_axes(ax)
        plt.show()



if __name__ == "__main__":
    proteom_composition_aa = 'proteome_composition_AA.csv'
    dico_aa = read_freq (proteom_composition_aa)
    amino = ['A', 'K', 'N', 'L']
    list_aa = ['V', 'I', 'L', 'M', 'F', 'W', 'Y', 'K', 'R', 'H', 'D', 'E', 'N', 'S', 'Q',  'T', 'P', 'G', 'A',
     'C']
    dico_min, dico_med, dico_max = get_name_min_med_max(dico_aa,amino, list_aa)
    my_plot = polar_plotting(dico_min, dico_med, dico_max, list_aa)
    # polar_plotting (dico_aa, names)
