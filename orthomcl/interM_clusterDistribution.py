#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Distribution des clusters intermediaires


usage
========

python3 interM_clusterDistribution.py
'''



import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

def count_fasta_header(directory, liste_cluster, biom):
    '''extract fasta sequence
    -in: -dir fasta files
    -out: dico: key: name_protein ,  value: sequence
    '''
    dict_cluster = {}
    list_protein_name = []
    for cluster in liste_cluster:
        liste_protein = []
        nb= 0
        if cluster.endswith('.fasta'):
            path_cluster = os.path.join(directory, cluster)
            cluster_name = cluster.split('.')[0]
            for line in open(path_cluster):
                line = line.strip()
                if line.startswith('>'):
                    nb+=1
                    name_protein = line[1:].split('|')[0]
                    liste_protein.append(name_protein)
        set_protein = set (liste_protein)
        dict_cluster [cluster_name] = len(set_protein)
    return dict_cluster


def plotting(dict_cluster):
    '''renvoit un barplot du nombre de clusters en fonction du nombre d'especes contenus dans les clusters
    '''

    data = {}
    my_list = defaultdict(list)
    for key, num in sorted(dict_cluster.items()):
        if num in my_list:
            my_list[num] += ( ' ' + key)
        else:
            my_list[num] = key

    list_values = []
    list_cluster_name = []

    for value in sorted(my_list):
        cluster_num=[]
        new_value = value - 6
        list_values.append(str(new_value))
        for number in my_list[value].split():
            tmp = number.split('_')
            num = tmp[1]
            cluster_num.append(num)
        list_cluster_name.append(int(len(cluster_num)))


    # for i in range(0, len(list_cluster_name), 150):
    #     cluster_lab = list_cluster_name[i : i +150]
    #     cluster_values = list_values[i : i +150]

    fig, ax = plt.subplots()
    N = len(list_cluster_name)
    ind = np.arange(N)
    width = 0.25

    dom = ax.bar(ind, list_cluster_name, width=0.25, alpha=0.4, color='b', label="+ : nombre de proteines (ou d'especes) contenu dans les clusters \n en plus des 6 cyano biominérales")
    ax.set_xlim(-width,len(ind)+width)
    ax.set_xticks(ind)
    labels = []
    for x_label in list_values:
        lab = '+'+x_label
        labels.append(lab)


    ax.set_xticklabels (labels, fontsize=10, color= 'b')

    plt.grid(True)
    plt.ylabel("Number of clusters", fontsize=15, color='k', alpha=0.8)
    plt.xlabel('6 cyano + X species', fontsize=15, color='b', alpha=0.8)
    for a,b in zip(ind, list_cluster_name):
        plt.text(a, list_cluster_name[a], str(b), va = 'bottom', fontsize=11, fontdict={'family': 'serif', 'color':  'k', 'weight': 'normal','size': 12})
    plt.title(u"Distribution du nombre de proteines dans les clusters de cyanobacteries biomineranlisantes", fontsize=17, fontdict={'family': 'monospace'})

    plt.show()

if __name__ == '__main__':

    directory = '/home/issa/Documents/stage/orthomcl/Intermediaires_clusters/clusters_Biominerales/'
    liste_cluster =  os.listdir(directory)
    biom = ['Synechococcus_sp_PCC_6312', 'Synechococcus_calcipolaris', 'Thermosynechococcus_elongatus_BP1', 'Gloeomargarita_lithophora', 'Cyanothece_sp_PCC_7425', 'Chroococcidiopsis_thermalis_PCC_7203']
    my_cluster = count_fasta_header(directory, liste_cluster, biom)
    plotting (my_cluster)
