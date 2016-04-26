#!/usr/bin/env python
# -*- coding: utf-8 -*-

''' Orphan domains analyses'''


# import sys
# from collections import defaultdict
# import matplotlib
# import matplotlib.pyplot as plt
# import numpy as np
# import random


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

# def write_output(dico_cl_no_dupli, outF):
#     """write in an output file the cluster_id, list and number of protein in cluster
#     """
#     tmp = []
#     with open(outF , 'w') as output:
#         for cluster in dico_cl_no_dupli:
#             list_protein = dico_cl_no_dupli[cluster]
#             output.write('>'+cluster+'\t'+str(len(list_protein))+'\n')
#             for protein, aa in list_protein:
#                 output.write(str(protein)+'\t'+str(aa)+'aa'+'\n')
def delete_duplicate(dico_cluster):
    ''' Remove the duplicated protein ids in clusters
    '''
    dico_cl_no_dupli = {}
    for cluster in dico_cluster:
        dico_cl_no_dupli[cluster]= []
        liste_protein = dico_cluster[cluster]
        l = set()
        g = []
        for prot, a in liste_protein:
            protein = prot
            aa = a
            if prot in l:
                # print('{prot} has more than one OD'.format(prot=protein))
                tpl= (protein, aa)
                g.append(tpl)
            else:
                l.add(protein)
                tpl = (protein, aa)
                dico_cl_no_dupli[cluster].append(tpl)
    return dico_cl_no_dupli

def protein_biomineral(fasta, dico_cluster, list_biominerale, list_non_biominerale):
    '''Return 2 dico with all proteins in the biominerales and in the non - biominerales species
    '''
    dico_fasta = {}
    with open(fasta) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                protein_id = line[1:]
                dico_fasta[protein_id]= ''
            else:
                seq = line
                dico_fasta[protein_id]=seq

    dico_biominerale = {}
    nb=0
    for cluster in dico_cluster:
        biom_count = 0
        no_biom_count = 0
        cluster_set = set()
        # dico_biominerale[cluster] = []
        list_protein = dico_cluster[cluster]
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

        if biom_count >= 4 and no_biom_count <= 1:
            nb+=1
            # with open(output, 'w') as out:
    #         print(cluster, len(list_protein))
            # print(protein, aa)
    #
    # print('nombre de cluster =' , nb)



# def write_fasta (fasta, cluster_biom)
def plotting (dico_cluster, dico_cl_no_dupli):
    '''plot the number of cluster according the number of protein in each cluster
    '''

    dict_lenght = defaultdict(list)
    my_list = defaultdict(list)
    for cluster in dico_cl_no_dupli:
        list_protein = dico_cl_no_dupli[cluster]
        nb_prot_in_cluster = len(list_protein)
        my_list[str(nb_prot_in_cluster)].append(cluster)

    # new_list = sorted(my_list, key=int) #= '1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '2', '20', '21', '22', '23', '24',...


    pe = sorted(my_list, key=int)
    number_of_protein = number_of_cluster = []
    for val in my_list:
        number_of_protein.append(val)
        lenght = len(my_list[val])
        number_of_cluster.append(lenght)
        nums = [x for x in range(10)]
        if val in nums:
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
    list_biominerale = ['Synechococcus_sp_PCC_6312', 'Synechococcus_calcipolaris', 'Thermosynechococcus_elongatus_BP1', 'Gloeomargarita_lithophora', 'Cyanothece_sp_PCC_7425', 'Chroococcidiopsis_thermalis_PCC_7203']
    list_non_biominerale = ['Gloeobacter_violaceus_PCC7421', 'Synechococcus_elongatus_PCC6301', 'Synechocystis_sp_PCC_7509', 'Synechocystis_sp_PCC6803', 'Gloeocapsa_sp_PCC_7428', 'Gloeobacter_kilaueensis_JS1']
    fasta = '/home/issa/Documents/stage/cd-hit/30result.fasta'
    output = '/home/issa/Documents/stage/cd-hit/biominerales_clusters/'
    cdhit30 = '/home/issa/Documents/stage/cd-hit/30result.fasta.clstr'
    # outF = sys.argv[2]

    cd_hit = read_cdhit_out(cdhit30)
    # cd_hit_no_dupli = delete_duplicate(cd_hit)
    # write_output(cd_hit_no_dupli, outF)
    # plotting(cd_hit, cd_hit_no_dupli)
    protein_biomineral(cd_hit, list_biominerale, list_non_biominerale, fasta)
