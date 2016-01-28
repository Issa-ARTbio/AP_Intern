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
                fr = dico_freq[aa2][name_med_proteome]
                my_dico_max[aa].append(fr)
    # print ('my_dico_min', my_dico_min, 'my_dico_med',my_dico_med,  'my_dico_med', my_dico_med)

    return my_dico_min, my_dico_med,  my_dico_max


def polar_plotting (my_dico_min, my_dico_med,  my_dico_max, list_aa):

    # fig = super(PlotWindPowerDensity, self).get_figure()
    for aa in my_dico_min:
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


        r = np.arange(20)
        ax = plt.subplot(111, polar=True)
        data_min = sorted (data_min[1:])
        ax.fill(data_min[1:],  facecolor='r',linewidth=3,alpha= 0.5, label = data_min[0])
        ax.fill(data_med[1:],  facecolor='m',linewidth=3,alpha= 0.5, label = data_min[0])
        ax.fill(data_max[1:],  facecolor='g',linewidth=3,alpha= 0.5, label = data_min[0])
#test : http://stackoverflow.com/questions/33072450/matplotlib-polar-plot-axis-tick-label-location
        # ax.plot ( color='black',linewidth=3)
        plt.grid(True)
        ax.set_xticks(r)
        ax.set_xticklabels (list_aa)
        ax.set_title('Plot de la frequence de '+ aa, va ='bottom', color = 'k')

        # fig.add_axes(ax)
        plt.show()



if __name__ == "__main__":
    proteom_composition_aa = 'proteome_composition_AA.csv'
    dico_aa = read_freq (proteom_composition_aa)
    amino = ['A', 'K', 'N', 'L']
    list_aa = ['A', 'R', 'N', 'D', 'C', 'Q', 'E','G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    dico_min, dico_med, dico_max = get_name_min_med_max(dico_aa,amino, list_aa)
    my_plot = polar_plotting(dico_min, dico_med, dico_max, list_aa)
    # polar_plotting (dico_aa, names)
