#!/usr/bin/env python
# -*- coding: utf-8 -*-


'''Extraire les noms des proteines dans les clusters intermediaires'''



import os, sys,re


def read_All_fasta(directory, liste_fasta):
    '''extract fasta sequence corresponding to the names in core file
    '''
    dict_proteome = {}
    list_protein_name = []
    for filename in liste_fasta:
        if filename.endswith('.fasta'):
            path_proteome = os.path.join(directory, filename)
            proteome_name = filename.split('.')[0]
            for line in open(path_proteome):
                line = line.strip()
                if line.startswith('>'):
                    name_protein = line[1:].split()[0]
                    # print(name_protein)
                    if name_protein in dict_proteome:
                        print("Error, protein {} already exists".format(name_protein))
                        sys.exit(1)
                    dict_proteome[name_protein] = ""
                else:
                    dict_proteome[name_protein] += line.strip()

    return dict_proteome

def read_mcl_output(groups, biom, proteomes_name):

    biom_cluster = []
    biomin = set(biom)
    # cluster_all, cluster_inter, cluster_spec, cluster_spec_single = [0]*len(proteomes_name), [0]*len(proteomes_name), [0]*len(proteomes_name), [0]*len(proteomes_name)
    with open (groups, 'r') as clusters:
        nb_cluster=0
        for line in clusters:
            line = line.rstrip()
            tmp = line.split()
            my_cluster = (tmp[1:])
            cluster_id = tmp[0]
            nb_cluster+=1
            list_cluter_inter = []
            tmp_lis = []
            for name in my_cluster:
                tmpr = name.split('|')
                proteome = tmpr[0]
                tmp_lis.append(proteome)
                list_cluter_inter.append(name)

            if 1 < len(list_cluter_inter) < len(proteomes_name):
                my_set = set(tmp_lis)
                if my_set >= biomin:
                    biom_cluster.append(list_cluter_inter)

    return biom_cluster

def write_cluster_file (biom_cluster, dict_proteome, outf):

    nb = 0
    for line in biom_cluster:
        nb+=1
        with open(outf+str(nb)+'.fasta' , 'w' ) as outfile:
            for indice in range(len(line)):
                protein = line[indice]
                for name_protein in dict_proteome:
                    if name_protein==protein:
                        print(name_protein)
                        protein_name = '>'+str(protein)
                        outfile.write(protein_name+'\n')
                        outfile.write(dict_proteome[name_protein]+'\n')



if __name__ == '__main__':

    outf = '/home/issa/Documents/stage/orthomcl/Intermediaires_clusters/sequences/cluster_'
    # groups = 'test_interm_cluster.txt'
    groups = '/home/issa/Documents/stage/orthomcl/orthomcl_results/groups.txt'
    # biom = ['Synechococcus_sp_PCC_6312', 'Synechococcus_calcipolaris', 'Thermosynechococcus_elongatus_BP1', 'Gloeomargarita_lithophora', 'Cyanothece_sp_PCC_7425', 'Chroococcidiopsis_thermalis_PCC_7203']
    biom = ['S_6312', 'S_calcipolaris', 'T_BP1', 'G_lithophora', 'C_7425', 'C_thermalis7203']
    proteome_dir = '/home/issa/Documents/stage/init_data/proteomes/'
    full_names = os.listdir(proteome_dir)
    directory = '/home/issa/Documents/stage/orthomcl/proteomes_format_shortname/'
    liste_fasta = os.listdir(directory)

    proteomes_name=[]

    for filename in liste_fasta:
        if filename.endswith('.fasta'):
            proteome_name = filename.split('.')[0]
            proteomes_name.append(proteome_name)

    fas = read_All_fasta(directory, liste_fasta)
    interm = read_mcl_output(groups, biom, proteomes_name)
    write_cluster_file (interm, fas, outf)
