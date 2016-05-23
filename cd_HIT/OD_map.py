#!/usr/bin/env python
# -*- coding: utf-8 -*-



"""Positionnement des domaines orphelins sur la protéine entière
"""

import sys, os
from collections import defaultdict
import matplotlib
import matplotlib.pyplot as plt
import numpy as np


def read_table(file_in):

    dico_data = defaultdict(list)
    with open(file_in) as f:
        for line in f:
            line = line.strip()
            element = line.split(',')
            protein_id = element[0]
            long_protein = element[1]
            sp = element[2]
            pos_aln = element[3]
            pos_fasta = element[4]
            len_DO = element[5]
            cluster_id = element[6]
            nb_prot_in_cluster = element[7]
            lenght_aln = element[8]
            tpl = (protein_id, pos_aln, len_DO, lenght_aln)
            dico_data[cluster_id].append(tpl)

    return dico_data

def plotting(dico_data):

    for cluster in dico_data:
        protein_lab = []
        data = dico_data[cluster]
        nb=1
        x1 = []
        x2 = []
        y = []
        for protein_id, pos_aln, len_DO, lenght_aln in data:
            nb+=1
            protein_lab.append(protein_id)
            x1.append(int(pos_aln))
            end_dom = int(pos_aln)+int(len_DO)
            x2.append((end_dom))
            y.append(nb)
            # print(protein_id,pos_aln, end_dom)

        #plt.plot((x1, x2), (y1, y2), 'k-') to draw a line from the point (x1, y1) to the point (x2, y2) in color k
        N = int(lenght_aln)
        ind = np.arange(N)
        fig = plt.figure()
        fig.subplots_adjust(bottom=0.2)
        ax = fig.add_subplot(111)
        plt.plot((x1, x2), (y,y), linewidth=3)
        l = int(lenght_aln)-1
        plt.plot((1, l), (0.5, 0.5), 'k-', linewidth=4.5, label='longueur alignement')
        ax.set_ylim(0,len(protein_lab)+2)
        # ax.set_xlim(0,N+5)
        plt.ylabel('Domaine Identifiant', fontsize=15, color='k', fontdict={'family': 'monospace'}, alpha=0.8)
        plt.xlabel('Domaines couverture', fontsize=15, color='k', fontdict={'family': 'monospace'}, alpha=0.8)
        ax.set_xticks(range(0,N,10))
        ax.set_yticks(range(2,len(x1)))
        ax.set_yticklabels(protein_lab)

        print(cluster, protein_lab)
        plt.title(u"Position des domaines orphelins sur l'alignement des proteines \n Pour le Cluster "+str(cluster), fontsize=17, fontdict={'family': 'monospace'})
        # plt.show()



if __name__ == '__main__':

    file_in = '/home/issa/Documents/stage/cd-hit/positions_DO.csv'
    data = read_table(file_in)
    plotting(data)
