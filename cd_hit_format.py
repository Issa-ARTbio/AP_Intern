#!/usr/bin/env python
# -*- coding: utf-8 -*-

''' Orphan domains analyses'''


import sys
from collections import defaultdict
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import random


def read_cdhit_out (cdhit30):
    '''read the .clstr file of the cd-hit output.
        return a dict with the protein_id in the cluster
    '''
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
                    protein = (protein_simil, aa)
                    dico_cluster[cluster_id].append(protein)
    return dico_cluster


def protein_biomineral(dico_cluster, list_biominerale, list_non_biominerale):
    '''Return 2 dicos with all proteins in the biominerales and in the non - biominerales species
    '''
    dico_biominerale = {}
    nb=0
    for cluster in dico_cluster:
        list_protein = dico_cluster[cluster]
        l = len(list_protein)

        biom_count = 0
        no_biom_count = 0
        cluster_set = set()

        for protein, aa in list_protein:
            proteome = protein.split('|')[1]
            cluster_set.add(proteome)
        set_pro = set()
        for proteome in cluster_set:

            if proteome in list_biominerale:
                biom_count +=1
                set_pro.add(proteome)
            elif proteome in list_non_biominerale:
                no_biom_count +=1

        if biom_count >= 4 and no_biom_count <= 0:
            nb+=1
            # print(set_pro)
            # print(cluster, l)
            dico_biominerale[cluster] = list_protein

    return dico_biominerale
    # print('nombre de cluster =' , nb)

def find_ortho(dico_cluster, dico_biominerale):

    dico_search={}
    inter = set(dico_cluster) - set(dico_biominerale)
    for cluster in dico_cluster:
        if cluster in inter:
            dico_search[cluster]= dico_cluster[cluster]
    nb=0
    found_cluster = {}
    for list_prot_biom in dico_biominerale.values():
        for protein, aa in list_prot_biom:
            for cluster in dico_search:
                protein_set = set(dico_search[cluster])
                for p, a in protein_set:
                    if p == protein and len(protein_set) > 1:
                        nb+=1
                        found_cluster[cluster] = protein_set
                        # print(cluster)
                        # print('Protein =', protein)
                        # print('Lenght in biominerale cluster', aa)
                        # print('Lenght in found cluster        ',a)
                        # print('*********************************')
    # print('Total protein found = ', nb)

    return found_cluster

def write_fasta(dico_biominerale, found_cluster, fasta, output):

    dico_fasta={}
    with open(fasta) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                protein_id = line[1:]
                dico_fasta[protein_id]= ''
            else:
                seq = line
                dico_fasta[protein_id]=seq

    for cluster in dico_biominerale:
        list_protein = dico_biominerale[cluster]
        nb=0
        with open(output+str(cluster)+'.fasta', 'w') as out:
            for protein , aa in list_protein:
                if protein in dico_fasta:
                    seq = dico_fasta[protein]
                    if len(seq) == int(aa):
                        nb+=1
                        out.write('>'+str(protein)+'\n'+str(seq)+'\n')
                        print(protein, len(seq), aa)
                        # print(seq)
                        # print(cluster)
        print(cluster, 'contains', nb, 'proteins')
        print("++++++++++++++++++++++++++++++++++++")
if __name__ == '__main__':
    cdhit30 = '/home/issa/Documents/stage/cd-hit/30result.fasta.clstr'
    list_biominerale = ['Synechococcus_sp_PCC_6312', 'Synechococcus_calcipolaris', 'Thermosynechococcus_elongatus_BP1', 'Gloeomargarita_lithophora', 'Cyanothece_sp_PCC_7425', 'Chroococcidiopsis_thermalis_PCC_7203']
    list_non_biominerale = ['Gloeobacter_violaceus_PCC7421', 'Synechococcus_elongatus_PCC6301', 'Synechocystis_sp_PCC_7509', 'Synechocystis_sp_PCC6803', 'Gloeocapsa_sp_PCC_7428', 'Gloeobacter_kilaueensis_JS1']
    fasta = '/home/issa/Documents/stage/cd-hit/data2.fasta'
    output = '/home/issa/Documents/stage/cd-hit/tmp/'

    cd_hit = read_cdhit_out(cdhit30)
    cl_biom = protein_biomineral(cd_hit, list_biominerale, list_non_biominerale)
    cl_orth = find_ortho(cd_hit, cl_biom)
    outp = write_fasta(cl_biom, cl_orth, fasta, output)
