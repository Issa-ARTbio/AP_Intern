#!/usr/bin/env python
# -*- coding: utf-8 -*-


'''Le pourcentage de gap dans les alignements
'''

import os, sys,re
import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter


def read_fasta (directory, liste_cluster):
    dict_cluster = {}
    for cluster in liste_cluster:
        if cluster.endswith('.fasta'):
            path_cluster = os.path.join(directory, cluster)
            cluster = cluster.split('.')[0]

            dict_proteome = {}
            for line in open(path_cluster):
                line = line.strip()
                if line.startswith('>'):
                    name_protein = line[1:].split()[0]

                    if name_protein in dict_proteome:
                        print("Error, protein {} already exists".format(name_protein))
                        sys.exit(1)
                    dict_proteome[name_protein] = ""

                else:
                    dict_proteome[name_protein] += line.strip()
            dict_cluster[cluster] = dict_proteome
    return dict_cluster

def gap_counter_plot (dict_cluster):
    dico_gap = {}
    for cluster in dict_cluster:
        dict_proteome = dict_cluster[cluster]
        gap_pourcen = []
        total = 0

        for protein in dict_proteome:
            sequence = dict_proteome[protein]
            gap = sequence.count('-')
            normaliz = gap/len(sequence)*100
            total += normaliz
        gap_coverage = total/len(dict_proteome)
        dico_gap[cluster] = float(gap_coverage)

    dico_gap = sorted(dico_gap.items(), key=itemgetter(1))
    list_gap_values, list_proteome_names = [], []
    for tupl in dico_gap:
        names, gap = tupl[0], tupl[1]
        list_proteome_names.append(names)
        list_gap_values.append(gap)


    for i in range(0, len(list_proteome_names), 150):
        cluster_lab = list_proteome_names[i : i +150]
        gap_pourc = list_gap_values[i : i +150]


        fig, ax = plt.subplots()
        N = len(cluster_lab)
        ind = np.arange(N)
        width = 0.35
        dom = ax.bar(ind, gap_pourc, width=0.6, alpha=0.4, color='r', label= 'Proteines conservees')
        ax.set_xlim(-width,len(ind)+width)
        ax.set_xticks(ind+width)
        ax.set_xticklabels (cluster_lab, rotation='vertical', fontsize=10)
        plt.grid(True)
        plt.ylabel('Pourcentage de Gaps', fontsize=15, color='r', alpha=0.8)
        plt.xlabel('Clusters', fontsize=15, color='b', alpha=0.8)
        ax.set_ylim(0,100)
        for a,b in zip(ind,gap_pourc):
            b = "%.2f"%b
            plt.text(a, b, (b)+'%', va = 'bottom',fontsize=10, rotation='vertical', fontdict={'family': 'serif', 'color':  'k', 'weight': 'normal','size': 16})
        plt.title(u"Pourcentage de Gap dans l'alignement des clusters \n {sur les blocks retenus par Gblock}", fontsize=17, fontdict={'family': 'monospace'})
        # ax.legend(loc='upper left', fontsize=10)
        plt.show()

if __name__ == '__main__':
    # path_proteome = 'test_interm_cluster.txt'
    directory = '/home/issa/Documents/stage/raxml/clusters_Gblock/alignment/fasta_cluster/'
    liste_cluster = os.listdir(directory)
    fasta = read_fasta(directory, liste_cluster)
    gap_counter_plot (fasta)
