#!/usr/bin/env python
# -*- coding: utf-8 -*-


''' A partir de la sortie du mcl programme:
repartit les differents clusters en fonction du nombre de proteome presents
clusters avec tous les proteomes representes (conserves)
clusters avec au moins 2 proteomes et au plus n-1 avec n = nombre de proteomes (111 dans notre cas) (intermediaires)
cluster specifique avec un seul proteome represente


Parametres:

-i  proteome names (.list) file groups.txt -o [histogramme]


usage
==============
python cluster_repartition.py
'''

import os, sys,re
import matplotlib
import matplotlib.pyplot as plt
import numpy as np


def count_clusters(good_prot, proteomes_name, groups):
    ''''Description de la fonction' : renvoit 3 listes avec les nombres de clusters contenant des proteines conservees, intermediaires et specifiques
    '''
    with open (groups, 'r') as out_mcl, open (good_prot, 'r') as f:

        all_proteins = set()
        my_ortho_cluster = set()
        for ortho in out_mcl :

            ortho = ortho.rstrip()
            tmp = ortho.split()
            my_cluster = (tmp[1:])
            for names in my_cluster:
                my_ortho_cluster.add(names)

        for line in f:
            line = line.strip()
            if line.startswith('>'):
                name_protein = line[1:].split()[0]
                all_proteins.add(name_protein)

        singleton = all_proteins - my_ortho_cluster
        # Format des singleton idem que celui d'orthomcl pour fichier de sortie (loop optionnel)
        # specifiq = set()
        # singleton=list(singleton)
        # for seq_id in range(len(singleton)):
        #     my_format = 'cyano'+str(seq_id)+': '+str(singleton[seq_id])
        #     specifiq.add(my_format)


        singleton = list(singleton)
        prot_spec=[]

        cluster_all, cluster_inter, cluster_spec = [0]*len(proteomes_name), [0]*len(proteomes_name), [0]*len(proteomes_name)

        for names in singleton:
            tmpr_sing = names.split('|')
            proteome_single = tmpr_sing[0]
            protein  = tmpr_sing[1]
            pos_sing = proteomes_name.index(proteome_single)
            prot_spec.append(proteome_single)
            cluster_spec[pos_sing]+=1

        with open(groups, 'r') as clusters:

            nb_cluster=0
            for line in clusters:
                line = line.rstrip()
                tmp = line.split()
                my_cluster = (tmp[1:])
                cluster_id = tmp[0]
                nb_cluster+=1

                dico_proteome_name = {}

                for names in my_cluster:
                    tmpr = names.split('|')
                    proteome = tmpr[0]
                    protein  = tmpr[1]

                    if proteome in dico_proteome_name:
                        dico_proteome_name[proteome]+=1
                    else:
                        dico_proteome_name[proteome]=1
                for proteome in dico_proteome_name:
                    for pos, item in enumerate(proteomes_name):
                        if item == proteome:
                            if len(dico_proteome_name)<=1:
                                cluster_spec[pos]+= dico_proteome_name[proteome]
                            elif len(dico_proteome_name) >= len(proteomes_name):
                                cluster_all[pos]+=dico_proteome_name[proteome]
                            else:
                                cluster_inter[pos]+=dico_proteome_name[proteome]

        return cluster_all, cluster_inter, cluster_spec

def give_full_name(full_names, proteomes_name):
    ''''Description de la fonction' : retransforme les noms des Proteomes (sous format lisible par orthomcl) aux noms complets initiaux
    '''

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
    # print(len(label_names))
    return label_names



def plotting (label_names, cluster_all, cluster_inter, cluster_spec):

    ''''Description de la fonction' : renvoit un histogramme superpose aux nombres de clusters avec des proteines conservees intermediaires et specifiques
    '''
    fig, ax = plt.subplots()
    N = len(label_names)
    ind = np.arange(N)
    width = 0.35

    dom = ax.bar(ind, cluster_all, width=0.6, alpha=0.4, color='b', label= 'conserves')
    dom = ax.bar(ind, cluster_inter, bottom=cluster_all, width=0.7, alpha=0.6, color='grey', label= 'intermediaires')
    bottom_spec = []
    for x,y in zip(cluster_all,cluster_inter):
        z=x+y
        bottom_spec.append(z)
    dom = ax.bar(ind, cluster_spec, bottom=bottom_spec, width=0.7, alpha=0.7, color='red', label= 'specifiques')
    all_values=[]
    for x,y in zip(bottom_spec,cluster_spec):
        z=x+y
        all_values.append(z)
    # print((all_values))
    ax.set_xlim(-width,len(ind)+width)
    ax.set_xticks(ind+width)
    ax.set_xticklabels (label_names, rotation='vertical', fontsize=13)
    for name in ax.set_xticklabels(label_names, rotation= 'vertical', fontsize=13):
        if re.search('Synechococcus_sp_PCC_6312', str(name)) :
            name.set_color('green')
        elif re.search('Synechococcus_calcipolaris' , str(name)):
            name.set_color('green')
        elif re.search('Thermosynechococcus_elongatus_BP1' , str(name)):
            name.set_color('green')


        elif re.search('Gloeomargarita_lithophora', str(name)):
            name.set_color('red')
        elif  re.search ('Cyanothece_sp_PCC_7425', str(name)):
            name.set_color('red')
        elif  re.search ('Chroococcidiopsis_thermalis_PCC_7203', str(name)):
            name.set_color('red')
    plt.grid(True)
    plt.ylabel('Nombre de proteines', fontsize=30, color='r', alpha=0.8)
    plt.xlabel('Proteomes', fontsize=30, color='b', alpha=0.8)
    # ax.set_ylim(0,500)
    print(sum(cluster_all))
    print(sum(cluster_inter))
    print(sum(cluster_spec))
    for a,b in zip(ind,all_values):
        plt.text(a, b, str(b), va = 'bottom',fontsize=10, rotation='vertical', fontdict={'family': 'serif', 'color':  'k', 'weight': 'normal','size': 16})
    # plt.title(u"Le nombre total de proteines par proteome dans les clusters", fontsize=17, fontdict={'family': 'monospace'})
    ax.legend(loc='upper left', fontsize=32)
    plt.show()

if __name__ == '__main__':

    proteome_dir = '/home/issa/Documents/stage/initial_data/proteomes/'
    full_names = os.listdir(proteome_dir)
    good_prot = '/home/issa/Documents/stage/orthomcl/orthomcl_results/goodProteins.fasta'
    groups = '/home/issa/Documents/stage/orthomcl/orthomcl_results/groups.txt'

    directory = '/home/issa/Documents/stage/orthomcl/proteomes_format_shortname/'
    liste_fasta = os.listdir(directory)

    proteomes_name=[]

    for filename in liste_fasta:
        if filename.endswith('.fasta'):
            proteome_name = filename.split('.')[0]
            proteomes_name.append(proteome_name)

    core, inter, spec = count_clusters(good_prot,proteomes_name, groups)
    liste_names = give_full_name(full_names, proteomes_name)
    my_plot= plotting(liste_names, core, inter, spec)
