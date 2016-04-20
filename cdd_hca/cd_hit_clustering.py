#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys
from collections import defaultdict
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def read_cdhit_out (cdhit30):

    dico_cluster = {}
    with open(cdhit30) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                nb_prot=[]
                cluster_id = line[1:]
                dico_cluster[cluster_id] = []
            else:
                nb_prot.append(1)
                protein =()
                # output.write(str(cluster_id)+'\t'+str(sum(nb_prot))+'\n')
                if line.split()[-1] == '*':
                    nb = line.split('\t')[0]
                    match = line.split('\t')[-1]
                    aa = match.split()[0][:-3]
                    protein_represent = match.split()[1][1:-3]
                    protein = (protein_represent, aa)
                    dico_cluster[cluster_id].append(protein)
                else:
                    nb = line.split('\t')[0]
                    match = line.split('\t')[-1]
                    aa = match.split()[0][:-3]
                    protein_simil = match.split()[1][1:-3]
                    # protein.append(protein_represent)
                    # protein.append(aa)
                    protein = (protein_simil, aa)
                    dico_cluster[cluster_id].append(protein)

        return dico_cluster

def write_output(dico_cluster, outF):
    """write in an output file the cluster_id, list and number of protein in cluster
    """
    tmp = []
    with open(outF , 'w') as output:
        for cluster in dico_cluster:
            list_protein = dico_cluster[cluster]
            output.write('>'+cluster+'\t'+str(len(list_protein))+'\n')
            for protein, aa in list_protein:
                output.write(str(protein)+'\t'+str(aa)+'aa'+'\n')

def plotting (dico_cluster):
    '''plot the number of cluster according the number of protein in each cluster
    '''

    dict_lenght = {}
    my_list = defaultdict(list)
    for cluster in dico_cluster:
        list_protein = dico_cluster[cluster]
        nb_prot_in_cluster = len(list_protein)
        my_list[str(nb_prot_in_cluster)].append(cluster)

    new_list = sorted(my_list) #= '1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '2', '20', '21', '22', '23', '24',...

    print(new_list)
    number_of_protein = number_of_cluster = []
    for val in my_list:
        number_of_protein.append(val)
        lenght = len(my_list[val])
        print(val,'proteins dans', lenght, 'clusters')
        number_of_cluster.append(lenght)

    N = len(number_of_protein)
    ind = np.arange(N)
    width = 0.1
    fig, ax = plt.subplots()

    ax.plot(ind, number_of_cluster)

    plt.grid(True)
    ax.set_xticks(ind)
    ax.set_xticklabels (number_of_protein, fontsize=8)
    ax.set_ylabel(u'Number of clusters', fontsize=15, color='blue')
    # ax.set_ylim (0, 0.15)
    ax.set_xlabel(u'Number of species', fontsize=15, color='blue')
    ax.set_title(u"Number of proteins in orphan domain clusters", fontsize=17, fontdict={'family': 'monospace'})
    # plt.show()

if __name__ == '__main__':
    cdhit30 = sys.argv[1]
    outF = sys.argv[2]
    cd_hit = read_cdhit_out(cdhit30)
    # write_output(cd_hit, outF)
    plotting(cd_hit)
