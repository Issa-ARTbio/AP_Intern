#!/usr/bin/env python
# -*- coding: utf-8 -*-


'''Represente la distribution de différence de position dans la sortie Gblocks
'''

import os, sys,re
import matplotlib
import matplotlib.pyplot as plt
import numpy as np


def read_file(directory_htm, list_dir_htm):

    cluster_labels = []
    initial= []
    gblock = []
    pourcentage_pos =[]
    for file_name in list_dir_htm:
        # if file_name.endswith('.htm'):
        if file_name in cluster_test:
            pos = os.path.join(directory_htm, file_name)
            file_name = file_name.split('.')[0]

            with open(pos, 'r') as f:
                lines = f.read().splitlines()
                for i in range(len(lines)):
                    line = lines[i]
                    if line.startswith(my_search):
                        name = line.split()[5]
                        gb = name.split('/')[-1]
                        cluster_numb = gb.split('.')[0]
                        cluster_labels.append(str(cluster_numb))
                        number_pos_Gblock = line.split()[7]
                        gblock.append(int(number_pos_Gblock))
                        pourc = line.split()[9]
                        pourc = pourc[1:-1]
                        pourcentage_pos.append(pourc)
                        num_initial = line.split() [13]
                        initial.append(int(num_initial))

    return cluster_labels, initial, gblock, pourcentage_pos

def read_fasta_aln (directory_aln, list_dir_aln):
    '''extract fasta sequences from dir, convert sequence in one string
    Return a dico with proteome_id as key and all corresponding sequence in files as value
    '''
    dico_all = {}
    for cluster in list_dir_aln:
        counter_gap=0
        lenght_tot = 0
        if cluster.endswith('.aln'):
            path_cluster = os.path.join(directory_aln, cluster)
            cluster = cluster.split('.')[0]
            with open(path_cluster) as f:
                for line in f:
                    if line.startswith('>'):
                        line = line[1:]
                        full_proteome_name = line.split('|')[0]
                    else:
                        line = line.strip()
                        seq_lenght = len(line)
                        gap = line.count('-')
                        lenght_tot += seq_lenght
                        counter_gap += gap
                        sequence = lenght_tot - counter_gap

    return lenght_tot, counter_gap, sequence


def plotting_gblock (cluster_labels, initial, gblock, pourcentage_pos):

    all_cluster = []
    t_cluster_name=()
    for i in range(len(cluster_labels)):

        name_cl = cluster_labels[i]
        init = initial[i]
        gbloc = gblock[i]
        pourcrnt = pourcentage_pos[i]
        t_cluster_name = init, name_cl, gbloc, pourcrnt,

        all_cluster.append(t_cluster_name)

    all_cluster.sort()

    list_name_cl = []
    list_init = []
    list_gbloc = []
    list_pourcrnt = []
    for initial_pos, name, pos_gb,pourct in all_cluster:

        list_name_cl.append(name)
        list_init.append(int(initial_pos))
        list_gbloc.append(int(pos_gb))
        list_pourcrnt.append(pourct)


    for i in range(0, len(list_name_cl), 40):
        pos_init = list_init[i : i +40]
        pos_gb = list_gbloc[i : i +40]
        cluster_number = list_name_cl[i : i +40]
        pourcentage_pos = list_pourcrnt[i : i +40]

        fig, ax = plt.subplots()
        N = len(cluster_number)
        ind = np.arange(N)
        width = 0.35


        init = ax.bar(ind, pos_init, width=0.2, alpha=0.5, color='grey', label= "Nombre de position initial dans l'alignement")
        gb = ax.bar(ind+width, pos_gb, width=0.2, alpha=0.5, color='r', label= 'Nombre de positions restant apres Gblock')

        ax.set_xlim(-width,len(ind)+width)
        ax.set_xticks(ind+width)
        ax.set_xticklabels (cluster_number, rotation='vertical', fontsize=10)

        plt.grid(True)
        plt.ylabel('Nombre de position', fontsize=15, color='g', alpha=0.8)
        plt.xlabel('Clusters', fontsize=15, color='b', alpha=0.8)
        # ax.set_ylim(0,100)
        for a,b in zip(ind, pourcentage_pos):
            plt.text(a+0.35, pos_gb[a], str(b)+'%', va = 'bottom', fontsize=20, rotation='vertical', fontdict={'family': 'serif', 'color':  'b', 'weight': 'normal','size': 10})
        # plt.title(u"Distribution du nombre de position sur le Gblock", fontsize=17, fontdict={'family': 'monospace'})
        plt.title(u"Distribution du nombre de position sur le Gblock \n {taille block 10, moitie gap}", fontsize=17, fontdict={'family': 'monospace'})
        ax.legend(loc='upper left', fontsize=10)

        plt.show()


if __name__ == '__main__':

    directory_htm = '/home/issa/Documents/stage/muscle/Gblock/Gblock_10_half_gap/htm/'
    list_dir_htm  = os.listdir(directory_htm)
    directory_aln = '/home/issa/Documents/stage/muscle/final_alignement/'
    list_dir_aln  = os.listdir(directory_aln)
    my_search = 'New number of positions in'
    cluster_test= ['cluster_333.fasta.aln-gb.htm', 'cluster_120.fasta.aln-gb.htm', 'cluster_243.fasta.aln-gb.htm', 'cluster_179.fasta.aln-gb.htm', 'cluster_264.fasta.aln-gb.htm', 'cluster_114.fasta.aln-gb.htm', 'cluster_176.fasta.aln-gb.htm', 'cluster_94.fasta.aln-gb.htm']

    name, inition, gb, pourcentage = read_file (directory_htm, list_dir_htm)
    gap_c = read_fasta_aln(directory_aln, list_dir_aln)
    plotting_gblock (name, inition, gb, pourcentage)
