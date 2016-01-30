#!/usr/bin/env python
# -*- coding: utf-8 -*-




import os, sys, argparse, re
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import random
import itertools
from polar_plot import read_freq



def plot_biomin (liste_proteome_biom, dico, list_aa):
    """ return a polar plor for a list of proteome names
    """

    # markers = ['o', 'v', '^', '<', '>', 's', 'p']
    # colors = ['b', 'g', 'r', 'cyan', 'm', 'y']
    # col = itertools.cycle(colors)
    # couleur = next(col)
    # print (couleur)
    for proteome_name in liste_proteome_biom:

        list_freq = []
        for aa in list_aa:
            data_aa = dico[aa]
            if proteome_name in data_aa:
                fr_aa = data_aa[proteome_name]
                list_freq.append(fr_aa)
        # print (list_freq)


        ax = plt.subplot(111, polar=True)
        r = np.arange(len(list_aa))
        rad = (r*2*np.pi/len(r))
        rad = list(rad)
        rad.append(rad[0])



        list_freq = list(list_freq)
        list_freq.append(list_freq[0])

        # for col, mark in zip (np.arange(len(liste_proteome_biom)), markers):
        #     f=0
        ax.plot(rad, list_freq, linewidth=5, alpha= 1, label = proteome_name)#, marker=[markers[i] for i in range,  (len(markers))]), marker= mark)
        ax.fill (rad, list_freq, alpha= 0.2)

        plt.grid(True)
        rad = rad[:-1]
        ax.set_xticks(rad)
        ax.set_rmax(0.13)

        hydrophobe = (['V', 'I', 'L', 'M', 'F', 'W', 'Y'])
        polaire = (['K', 'R', 'H', 'D', 'E', 'N', 'S', 'Q',  'T'])
        petit_apolaire = (['P', 'G', 'A', 'C'])
        list_lab = (hydrophobe, polaire, petit_apolaire)
        colors = [ 'b', 'r', 'k']
        for xtick, i in zip (ax.get_xticklabels(), list_aa):

            if i in hydrophobe:
                xtick.set_color(colors[0])
            if i in polaire:
                xtick.set_color(colors[1])
            if i in petit_apolaire:
                xtick.set_color(colors[2])
                ax.set_xticklabels(list_aa)
    ax.set_title(u'Frequence des differents aa dans les proteomes des cyanobacteries biomineralisantes'+u'\n avec des inclusions de CaCO3 localisees aux poles', va ='bottom', color = 'g', fontdict={'family': 'monospace'})
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.10),
    fancybox=True, shadow=True, ncol=3)
    plt.show()


if __name__ == "__main__":

    list_aa = ['V', 'I', 'L', 'M', 'F', 'W', 'Y', 'K', 'R', 'H', 'D', 'E', 'N', 'S', 'Q',  'T', 'P', 'G', 'A',
     'C']
    liste_proteome_biom = ['Synechococcus_sp_PCC_6312', 'Synechococcus_calcipolaris', 'Thermosynechococcus_elongatus_BP1']
    # proteome_name = 'Synechococcus_sp_PCC_6312' 'Synechococcus_calcipolaris', 'Thermosynechococcus_elongatus_BP1','Gloeomargarita_lithophora', 'Cyanothece_sp_PCC_7425', 'Chroococcidiopsis_thermalis_PCC_7203'
    proteom_composition_aa = './old_/proteome_composition_AA.csv'
    dico = read_freq (proteom_composition_aa)
    # for proteome_name in liste_proteome_biom:
    freq_prot = plot_biomin(liste_proteome_biom, dico, list_aa)
