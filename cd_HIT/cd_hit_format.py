#!/usr/bin/env python
# -*- coding: utf-8 -*-

''' Orphan domains analyses'''


import sys, os
from collections import defaultdict
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import random
from interM_clusterExtraction import read_All_fasta

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

def find_orthoMCL(groups, dico_biominerale):

    nb=0
    found_cluster = set()
    with open(groups) as f:
        number_cluster = 0
        for line in f:
            number_cluster+=1
            nb_prot_biom = nb_prot_no_biom = nb_prot_total = 0
            line = line.strip()
            for name in line.split():
                nb_prot_total+=1
                proteome_SN = name.split('|')[0]
                protein = name.split('|')[1]


                for list_prot_biom in dico_biominerale.values():
                    for protein_id, aa in list_prot_biom:
                        tmp = protein_id.split('|')
                        p = tmp[0]
                        proteome = tmp[1]
                        if protein == p:
                            # z = name.replace(proteome_SN, proteome)
                            nb+=1
                            found_cluster.add(line)
    # print('Number of cluster OrthoMcl', nb)
    return found_cluster

def format_found_cluster(found_cluster, dico_biominerale, fasta, output):
    '''Give full name to the found proteins in the cluster
    And count the number of matched protein in orphan domain cluster
    '''
    nb=0
    bn=0
    tmp_set = set()
    for cluster in found_cluster:
        nb+=1
        # print(cluster)
        # with open(output+str(nb)+'.fasta', 'w') as outF:
        for protein_id in cluster.split():
            if protein_id in fasta:
                bn +=1
                tmp_set.add(protein_id)
                seq = fasta[protein_id]
                # outF.write('>'+str(protein_id)+'\n'+seq+'\n')
                # print(nb, '====>>>> ')
                # print(protein_id)
                # print(seq)
    print(len(tmp_set))
    print('number of orphan cluster =', nb)
    print('number of protein in fasta =', bn)

# def write_fasta(dico_biominerale, found_cluster, fasta, output):
#     '''write fasta sequences in output file
#     '''
#     dico_fasta={}
#     with open(fasta) as f:
#         for line in f:
#             line = line.strip()
#             if line.startswith('>'):
#                 protein_id = line[1:]
#                 dico_fasta[protein_id]= ''
#             else:
#                 seq = line
#                 dico_fasta[protein_id]+=seq
#
#     for cluster in dico_biominerale:
#         list_protein = dico_biominerale[cluster]
#         nb=bn=0
#         with open(output+str(cluster)+'.fasta', 'w') as out:
#             for protein , aa in list_protein:
#                 if protein in dico_fasta:
#                     seq = dico_fasta[protein]
#                     if len(seq) == int(aa):
#                         nb+=1
#                         out.write('>'+str(protein)+'\n'+str(seq)+'\n')
#                         # print(protein, len(seq), aa)
#                         # print(seq)
#                         # print(cluster)
#                     else: # protein ids are duplicated in the fasta file with differents segment lenghts...
#                         bn+=1
#                         print("Error:", cluster, protein, len(seq), aa)
#         print(cluster, 'contains', nb, 'proteins')
#         print(bn, 'proteins not write')
#         # print("++++++++++++++++++++++++++++++++++++")
if __name__ == '__main__':
    cdhit30 = '/home/issa/Documents/stage/cd-hit/30result.fasta.clstr'
    list_biominerale = ['Synechococcus_sp_PCC_6312', 'Synechococcus_calcipolaris', 'Thermosynechococcus_elongatus_BP1', 'Gloeomargarita_lithophora', 'Cyanothece_sp_PCC_7425', 'Chroococcidiopsis_thermalis_PCC_7203']
    list_non_biominerale = ['Gloeobacter_violaceus_PCC7421', 'Synechococcus_elongatus_PCC6301', 'Synechocystis_sp_PCC_7509', 'Synechocystis_sp_PCC6803', 'Gloeocapsa_sp_PCC_7428', 'Gloeobacter_kilaueensis_JS1']
    fasta = '/home/issa/Documents/stage/cd-hit/data2.fasta'
    output = '/home/issa/Documents/stage/cd-hit/biominerales_clusters/fasta_clusters_OrthoMCL/Cluster_'
    groups = '/home/issa/Documents/stage/orthomcl/Intermediaires_clusters/wo_para/Inter_Cluster_sansDupli.txt'
    directory = '/home/issa/Documents/stage/orthomcl/proteomes_format_shortname/'
    liste_fasta = os.listdir(directory)

    fasta = read_All_fasta(directory, liste_fasta)
    cd_hit = read_cdhit_out(cdhit30)
    # print(len(cd_hit["Cluster 677"]))
    cl_biom = protein_biomineral(cd_hit, list_biominerale, list_non_biominerale)
    cl_orth = find_orthoMCL(groups, cl_biom)
    # outp = write_fasta(cl_biom, cl_orth, fasta, output)
    format_found_cluster(cl_orth, cl_biom, fasta, output)
