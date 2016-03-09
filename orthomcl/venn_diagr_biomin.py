#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib_venn  import venn3, venn3_circles

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
        # print (liste_protein)
        set_protein = set (liste_protein)
        liste_protein = str(liste_protein)
        print (cluster_name, ':' , liste_protein, len(set_protein))
        if liste_protein in dict_cluster:
            dict_cluster[liste_protein]+=1
        else:
            dict_cluster[liste_protein]=1
        # dict_cluster [cluster_name] = set_protein

    # for x in dict_cluster:
    #
    #     print (cluster_name, ':' , x, dict_cluster[x])
    return dict_cluster


# def venn_diagr (dict_cluster):


if __name__ == '__main__':

    directory = '/home/issa/Documents/stage/orthomcl/Intermediaires_clusters/clusterBiomOnly/'
    liste_cluster =  os.listdir(directory)
    biom = ['Synechococcus_sp_PCC_6312', 'Synechococcus_calcipolaris', 'Thermosynechococcus_elongatus_BP1', 'Gloeomargarita_lithophora', 'Cyanothece_sp_PCC_7425', 'Chroococcidiopsis_thermalis_PCC_7203']
    my_cluster = count_fasta_header(directory, liste_cluster, biom)
    # plotting (my_cluster)
