#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Renvoit un polar_plot pour les frequences des acides amines dans la liste amino
la liste amino correspond aux 3 aa avec les plus grandes variations de la frequence sur l'ensemble des 111 proteomes
la Leucine (L) a ete rajoute car elle presente une freq plus importante dans tous les proteomes
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


def get_name_min_med_max (dico_freq, amino, list_aa):
    """renvoit le nom du proteome avec la frequence minimale, mediane et maximale de chaque aa present dans amino
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
        print (aa, name_min_proteome, val_trie[0],  name_med_proteome, val_trie[int(medium)],  name_max_proteome ,  val_trie[-1])

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
    '''renvoit un plot de la Frequence des differents aa dans les proteomes
    '''

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


        # A Prochlorococcus_marinus_MIT9515 0.052763944272729624 Fischerella_muscicola_PCC_7414 0.08219473507646308 Synechococcus_sp_WH5701 0.11434684337252429

        #  K Synechococcus_sp_WH5701 0.0216 67929809500178 Nodularia_spumigena_CCY9414 0.0484 7890822986619 Prochlorococcus_marinus_MIT9215 0.0871 7695568848365

        # N Synechococcus_sp_WH5701 0.021664576332951036 Synechocystis_sp_PCC_7509 0.04383814914250497 Prochlorococcus_marinus_MIT9515 0.06562362370740592

        # L Rivularia_sp_PCC_7116 0.1046 7297187323174 Microchaete_sp_PCC_7126 0.1100 5250461551017 Synechococcus_sp_WH5701 0.1319 4253482585394



        ax.plot(rad, data_min_s,  color='r',linewidth=5,alpha= 0.3, label = data_min[0]+'\n Frequence Minimale = 0.1046', marker= '*')
        ax.fill (rad, data_min_s, color= 'r', alpha= 0.5)

        ax.plot(rad, data_med_s,  color='b',linewidth=5,alpha= 0.3, label = data_med[0]+'\n Frequence Mediane = 0.11', marker= 's')
        ax.fill (rad, data_med_s, color= 'b', alpha= 0.5)

        ax.plot(rad, data_max_s,  color='g',linewidth=5,alpha= 0.7, label = data_max[0]+'\n Frequence Maximale = 0.1319', marker= '>')
        ax.fill (rad, data_max_s,color= 'g', alpha= 0.5)

        plt.grid(True)
        rad = rad[:-1]
        ax.set_xticks(rad)
        ax.set_rmax(0.16)

        hydrophobe = (['V', 'I', 'L', 'M', 'F', 'W', 'Y'])
        polaire = (['K', 'R', 'H', 'D', 'E', 'N', 'S', 'Q',  'T'])
        petit_apolaire = (['P', 'G', 'A', 'C'])
        list_lab = (hydrophobe, polaire, petit_apolaire)
        colors = [ 'b', 'r', 'k']

        for xtick, i in zip (ax.get_xticklabels(), list_aa):
            if i in petit_apolaire:
                xtick.set_color(colors[2])
            if i in polaire:
                xtick.set_color(colors[1])
            if i in hydrophobe:
                xtick.set_color(colors[0])
                ax.set_xticklabels(list_aa)

                #Rendre gras un caractere du label
                # if i == 'L':
                    # xtick.set_color('purple')
                    # xtick.set_weight(1000)
                    # xtick.set_size('x-large')

        ax.set_title('Frequence des differents aa dans les proteomes\n avec une grande variation de la frequence de '+aa, va ='bottom', color = 'k', fontdict={'family': 'monospace'})
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.10),
  fancybox=True, shadow=True, ncol=3)
        # fig.add_axes(ax)
        plt.show()



if __name__ == "__main__":
    proteom_composition_aa = '/home/issa/Documents/STAGE/Results/proteome_composition_AA.csv'
    dico_aa = read_freq (proteom_composition_aa)
    amino = ['A', 'K', 'N', 'L']
    list_aa = ['V', 'I', 'L', 'M', 'F', 'W', 'Y', 'K', 'R', 'H', 'D', 'E', 'N', 'S', 'Q',  'T', 'P', 'G', 'A',
     'C']
    dico_min, dico_med, dico_max = get_name_min_med_max(dico_aa,amino, list_aa)
    my_plot = polar_plotting(dico_min, dico_med, dico_max, list_aa)
    # polar_plotting (dico_aa, names)
