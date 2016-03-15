#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter
import pandas as ps


def read_hgtFile(directory_in, liste_in): #numberOfhgt

    dico = {}
    for file_name in liste_in:
        if file_name.endswith('.outscr'):
            hgt_in = os.path.join(directory_in, file_name)
            file_name = file_name.split('.')[0]

            with open(hgt_in, 'r') as f:
                lines = f.read().splitlines()
                for i in range(len(lines)):
                    line = lines[i]
                    if line.startswith(my_search):
                        # numberOfhgt.write(file_name)
                        tmp = line.split()
                        hgt = tmp[-1]
                        dico[file_name] = str(hgt)


    return dico
def plotting (dico):



    dataF = {}
    dico_ght = sorted(dico.items(), key=itemgetter(1))
    list_hgt_values, list_cluster_names = [], []
    for tupl in dico_ght:
        cluster, hgt = tupl[0], tupl[1]
        list_hgt_values.append(int(hgt))
        list_cluster_names.append(cluster)
    # print(cluster, hgt)
    # for i in range(0, len(list_cluster_names), 150):
    #     cluster_lab = list_cluster_names[i : i +150]
    #     num_hgt = list_hgt_values[i : i +150]
    #     dataF[num_hgt]= cluster_lab
    #     # dataF['name']= cluster_lab

    list_hgt_group = ps.DataFrame(dico_ght, columns=['name', 'hgt'])


    # print(list_hgt_group)

    hgt_plot, names = [], []
    grouped_plot_data = list_hgt_group.groupby('hgt', sort=True)
    for x in grouped_plot_data:
        hgt_plot.append(int(x[0]))
        names.append(x[1])
    #     for y in names:
    #         if y.find('cluster'):
    #             print(y)
    #             # print(x[1])
    #             print('=============')
    # # groups = dict()
    # data = list_hgt_group.groupby(['hgt', 'name'])
    # n = data.name
    n = dict(list(grouped_plot_data))
    for y in n:
        print(y)
        print('=============')
        print(n[y])
        print('========================================')





    fig, ax = plt.subplots()
    N = len(hgt_plot)
    ind = np.arange(N)
    width = 0.35
    dom = ax.bar(ind, hgt_plot, width=0.2, alpha=0.4, color='r')
    # # fig = list_hgt_group.plot(kind='hexbin', x = 'name', y = 'ght', color='k', alpha=0.5, width=0.6, label= 'Tranfert de Genes Horizontaux')
    #
    #
    ax.set_xlim(-width,len(ind)+width)
    ax.set_xticks(ind)
    ax.set_xticklabels (grouped_plot_data['name'], rotation='vertical', fontsize=10)
    plt.grid(True)
    plt.ylabel('Nombre de HGT', fontsize=15, color='r', alpha=0.8)
    plt.xlabel('Clusters', fontsize=15, color='b', alpha=0.8)
    ax.set_ylim(0,100)
    for a,b in zip(ind,hgt_plot):
        # b = "%.2f"%b
        plt.text(a, b, str(b), va = 'bottom',fontsize=10, rotation='vertical', fontdict={'family': 'serif', 'color':  'k', 'weight': 'normal','size': 16})
    plt.title(u"Nombre de HGT dans l'alignement des clusters \n {sur les arbres construits à partir de Gblock}", fontsize=17, fontdict={'family': 'monospace'})
    # ax.legend(loc='upper left', fontsize=10)
    plt.show()

if __name__ == '__main__':

    directory_in = '/home/issa/Documents/stage/hgt/result_hgt/'
    liste_in = os.listdir(directory_in)
    # outf = 'sortie_model_prottest.txt'
    my_search = 'hgt : number of HGT(s) found'

    num_ght = read_hgtFile(directory_in, liste_in)
    plotting(num_ght)
