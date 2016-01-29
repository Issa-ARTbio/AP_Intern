




import os, sys, argparse, re
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import random
from polar_plot import read_freq



def plot_biomin (liste_proteome_biom, dico, list_aa):
    """ return a polar plor for a list of proteome names
    """

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
        ax.plot(rad, list_freq,  color='b',linewidth=5,alpha= 0.3, label = proteome_name, marker= 's')
        ax.fill (rad, list_freq, color= 'b', alpha= 0.5)

        plt.grid(True)
        rad = rad[:-1]
        ax.set_xticks(rad)
        ax.set_rmax(0.15)

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
    ax.set_title('Frequence des differents aa dans les proteomes des cyanobacteries biomineralisantes', va ='bottom', color = 'k', fontdict={'family': 'monospace'})
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.10),
    fancybox=True, shadow=True, ncol=3)
    plt.show()


if __name__ == "__main__":

    list_aa = ['V', 'I', 'L', 'M', 'F', 'W', 'Y', 'K', 'R', 'H', 'D', 'E', 'N', 'S', 'Q',  'T', 'P', 'G', 'A',
     'C']
    liste_proteome_biom = [ 'Synechococcus_sp_PCC_6312', 'Synechococcus_calcipolaris', 'Thermosynechococcus_elongatus_BP1', 'Gloeomargarita_lithophora', 'Cyanothece_sp_PCC_7425', 'Chroococcidiopsis_thermalis_PCC_7203']
    # proteome_name = 'Synechococcus_sp_PCC_6312'
    proteom_composition_aa = 'proteome_composition_AA.csv'
    dico = read_freq (proteom_composition_aa)
    # for proteome_name in liste_proteome_biom:
    freq_prot = plot_biomin(liste_proteome_biom, dico, list_aa)
