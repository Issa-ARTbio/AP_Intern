#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np


def count_fasta_header(directory, liste_cluster):
    '''extract fasta sequence
    -in: -dir fasta files
    -out: dico: key: name_protein ,  value: sequence
    '''
    dict_cluster = {}
    list_protein_name = []
    for cluster in liste_cluster:
        nb= 0
        if cluster.endswith('fasta'):

            print(cluster)
            path_cluster = os.path.join(directory, cluster)
            cluster_name = cluster.split('.')[0]
            for line in open(path_cluster):
                line = line.strip()
                if line.startswith('>'):
                    nb+=1
                    name_protein = line[1:].split()[0]
        dict_cluster [cluster_name] = nb
    for cluster in dict_cluster:
        if cluster == 'cluster_118':
            print(dict_cluster[cluster])
    return dict_cluster


def plotting(dict_cluster):

    tpl = ()
    my_sort_list = []
    for cluster in dict_cluster:
        value = dict_cluster[cluster]
        tpl = value, cluster
        my_sort_list.append(tpl)


    my_sort_list.sort()

    cluster_label = []
    cluster_num = []
    for value, name in my_sort_list:
        cluster_num.append(int(value))
        cluster_label.append(name)



    for i in range(0, len(cluster_label), 150):
        cluster_lab = cluster_label[i : i +150]
        cluster_values = cluster_num[i : i +150]

        fig, ax = plt.subplots()
        N = len(cluster_lab)
        ind = np.arange(N)
        width = 0.35

        dom = ax.bar(ind, cluster_values, width=0.6, alpha=0.4, color='k')
        ax.set_xlim(-width,len(ind)+width)
        ax.set_xticks(ind+width)

        ax.set_xticklabels (cluster_lab, rotation='vertical', fontsize=10)

        plt.grid(True)
        plt.ylabel('Nombre de proteines', fontsize=15, color='r', alpha=0.8)
        plt.xlabel('Clusters', fontsize=15, color='b', alpha=0.8)
        for a,b in zip(ind, cluster_values):
            plt.text(a, cluster_values[a], str(b), va = 'bottom', fontsize=10, rotation='vertical', fontdict={'family': 'serif', 'color':  'red', 'weight': 'normal','size': 10})
        plt.title(u"Distribution du nombre de protein dans les clusters \n {Presence des proteomes de cyano biomineranlisantes dans tous les clusters}", fontsize=17, fontdict={'family': 'monospace'})
        plt.show()

if __name__ == '__main__':

    directory = '/home/issa/Documents/stage/orthomcl/Intermediaires_clusters/seq_full_names/'
    liste_cluster =  os.listdir(directory)
    my_cluster = count_fasta_header(directory, liste_cluster)
    plotting (my_cluster)
