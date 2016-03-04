#!/usr/bin/env python
# -*- coding: utf-8 -*-


'''Extraire les noms des proteines dans les clusters intermediaires'''



import os, sys,re


def read_All_fasta(directory, liste_fasta):
    '''extract fasta sequence
    -in: -dir fasta files
    -out: dico: key: name_protein ,  value: sequence
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
    '''renvoit la liste des clusters intermediaires dans la sortie du mcl
    '''
    biom_cluster = []
    biomin = set(biom)
    cluster_all, cluster_inter, cluster_spec = [0]*len(proteomes_name), [0]*len(proteomes_name), [0]*len(proteomes_name)
    with open (groups, 'r') as clusters:
        nb_cluster=0
        for line in clusters:
            line = line.rstrip()
            tmp = line.split()
            my_cluster = (tmp[1:])
            cluster_id = tmp[0]
            nb_cluster+=1

            #make a dict to regroup duplications
            dico_proteome_name = {}

            for names in my_cluster:
                tmpr = names.split('|')
                proteome = tmpr[0]
                protein  = tmpr[1]

                if proteome in dico_proteome_name:
                    dico_proteome_name[proteome]+=1
                else:
                    dico_proteome_name[proteome]=1

            #make a liste for intermediaire clusters
            list_cluter_inter = []
            tmp_liste = []
            for name in my_cluster:
                if 1 < len(dico_proteome_name) < len(proteomes_name)/2:
                    tmpr = name.split('|')
                    proteome = tmpr[0]
                    tmp_liste.append(proteome)
                    list_cluter_inter.append(name)
                    # Test if in intermediare clusters we have all biomineralizing
                    my_set = set(tmp_liste)
                    if my_set >= biomin and len(list_cluter_inter)< len(proteomes_name)/2:
                        biom_cluster.append(list_cluter_inter)

    # print(len(biom_cluster))
    return biom_cluster

def give_full_name(full_names, proteomes_name, biom_cluster):
    ''''Description de la fonction' : retransforme les noms des Proteomes (sous format lisible par orthomcl) aux noms complets initiaux
    '''

    new_biom_cluster = []
    proteomes_name= [thermalis.replace('C_thermalis7203','C_7203') for thermalis in proteomes_name]
    match_name={}
    label_names=[]
    for filename in full_names:
        if filename.endswith('.fasta'):
          filename = filename.split('.')[0]
          sub_id = filename.split('_')
          identifiant = sub_id[0][0]+'_'+sub_id[-1]
          match_name[identifiant]= filename
    for short_name in proteomes_name:
        full_name = match_name[short_name]
        label_names.append(full_name)


    for my_list in biom_cluster:
        tmp_dico = {}
        my_new_list = []
        for name in my_list:
            proteome = name.split('|')[0]
            protein  = name.split('|')[1]
            if proteome in match_name:
                tmp_dico[match_name[proteome]] = protein
                new_name = str(match_name[proteome])+'|'+str(protein)
                my_new_list.append(new_name)
            new_biom_cluster.append(my_new_list)
    return label_names, new_biom_cluster

def write_cluster_file (biom_cluster, dict_proteome, outf):

    nb = 0
    for line in biom_cluster:
        nb+=1
        with open(outf+str(nb)+'.fasta' , 'w' ) as outfile:
            for indice in range(len(line)):
                protein = line[indice]
                if protein in dict_proteome:
                    print(protein)
                    protein_name = '>'+str(protein)
                    outfile.write(protein_name+'\n')
                    outfile.write(dict_proteome[protein]+'\n')
                    print('OK')


if __name__ == '__main__':

    proteome_dir = '/home/issa/Documents/stage/init_data/proteomes/'
    full_names = os.listdir(proteome_dir)
    outf = '/home/issa/Documents/stage/orthomcl/Intermediaires_clusters/seqshort_names/cluster_'
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
    give_full_name(full_names, proteomes_name, interm)
    write_cluster_file (interm, fas, outf)