#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

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
    my_list = defaultdict(list)
    for key, num in sorted(dico.items()):
        if num in my_list:
            my_list[num] += ( ' ' + key)
        else:
            my_list[num] = key
    list_values = []
    list_cluster_name = []

    for value in sorted(my_list):
        cluster_num=[]
        # print(int(value))
        list_values.append(int(value))
        for number in my_list[value].split():
            tmp = number.split('_')
            num = tmp[1]
            cluster_num.append(num)
        # print(len(cluster_num))
        list_cluster_name.append((len(cluster_num)))
    print(sum(list_cluster_name))
    print(len(list_values))
    fig, ax = plt.subplots()
    N = len(list_cluster_name)
    ind = np.arange(N)
    width = 0.35
    dom = ax.bar(ind, list_cluster_name, width=0.2, alpha=0.4, color='r')
    # # fig = list_hgt_group.plot(kind='hexbin', x = 'name', y = 'ght', color='k', alpha=0.5, width=0.6, label= 'Tranfert de Genes Horizontaux')
    #
    #
    ax.set_xlim(-width,len(ind)+width)
    ax.set_xticks(ind)
    ax.set_xticklabels (list_values, fontsize=15)
    plt.grid(True)
    plt.ylabel('Nombre de Clusters ', fontsize=15, color='r', alpha=0.8)
    plt.xlabel('Nombre de HGTs', fontsize=15, color='b', alpha=0.8)
    # ax.set_ylim(0,100)
    for a,b in zip(ind,list_cluster_name):
        # b = "%.2f"%b
        plt.text(a, b, str( b), va = 'bottom',fontsize=13,fontdict={'family': 'serif', 'color':  'red', 'weight': 'normal','size': 19})
    plt.title(u"Nombre de HGT dans l'alignement des clusters \n {sur les arbres construits à partir de Gblock}", fontsize=17, fontdict={'family': 'monospace'})
    ax.legend(loc='upper left', fontsize=10)
    plt.show()

if __name__ == '__main__':

    directory_in = '/home/issa/Documents/stage/hgt/result_hgt/'
    liste_in = os.listdir(directory_in)
    # outf = 'sortie_model_prottest.txt'
    my_search = 'hgt : number of HGT(s) found'

    num_ght = read_hgtFile(directory_in, liste_in)
    plotting(num_ght)
