#!/usr/bin/env python
# -*- coding: utf-8 -*-


'''Find the protein id which matches with an orphan domain in the output cd-hit clustering
'''
import os, sys

def read_fasta_protein(directory, liste_fasta):

    dict_cluster_ortMCL = {}
    list_protein_name = []
    for filename in liste_fasta:
        if filename.endswith('.fasta'):
            path_cluster = os.path.join(directory, filename)
            cluster_id = filename.split('.')[0]
            for line in open(path_cluster):
                line = line.strip()
                if line.startswith('>'):
                    name_protein = line[1:].split('|')[1]
                    # print(name_protein)
                    if name_protein in dict_cluster_ortMCL:
                        print("Duplication, protein {prot} already exists in {cl}".format(prot=name_protein, cl=cluster_id))
                    #     sys.exit(0)
                    dict_cluster_ortMCL[name_protein] = 0
                else:
                    dict_cluster_ortMCL[name_protein] += len(line.strip())
    return dict_cluster_ortMCL


def read_fasta_OD(directory_orphanD, liste_faa, dict_cluster_ortMCL):

    dict_cluster_OD = {}
    list_protein_name = []
    nb=bn=0
    for filename in liste_faa:
        if filename.endswith('.fasta'):
            path_cluster = os.path.join(directory_orphanD, filename)
            cluster_id = filename.split('.')[0]
            for line in open(path_cluster):
                line = line.strip()
                if line.startswith('>'):
                    name_protein = line[1:].split('|')[0]
                    # print(name_protein)
                    if name_protein in dict_cluster_ortMCL:
                        nb+=1
                        print(name_protein, 'with', dict_cluster_ortMCL[name_protein],'aa contains an OD in', cluster_id)
                    else:
                        bn+=1
                    dict_cluster_OD[name_protein] = 0
                else:
                    dict_cluster_OD[name_protein] += len(line.strip())
                    # else:
                    #     print(name_protein)
    #                     print("Error, protein {} already exists".format(name_protein))
    #                 #     sys.exit(0)
    #                 dict_cluster_OD[name_protein] = 0
    #             else:
    #                 dict_cluster_OD[name_protein] += len(line.strip())
    # return dict_cluster_OD
    print(nb, 'matches and',bn,'not found' )
# def find_OD (dict_cluster):



if __name__ == '__main__':

    directory = '/home/issa/Documents/stage/cd-hit/biominerales_clusters/cluster_orthoMCL_found/'
    liste_fasta = os.listdir(directory)
    list_biominerale = ['Synechococcus_sp_PCC_6312', 'Synechococcus_calcipolaris', 'Thermosynechococcus_elongatus_BP1', 'Gloeomargarita_lithophora', 'Cyanothece_sp_PCC_7425', 'Chroococcidiopsis_thermalis_PCC_7203']
    directory_orphanD = '/home/issa/Documents/stage/cd-hit/biominerales_clusters/fasta/'
    liste_faa = os.listdir(directory_orphanD)
    othomcl = read_fasta_protein(directory, liste_fasta)
    read_fasta_OD(directory_orphanD, liste_faa,othomcl)
