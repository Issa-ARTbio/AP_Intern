#!/usr/bin/env python
# -*- coding: utf-8 -*-



""" Recuperer les clusters dont les NON - biominerales (avec phenotypes connus) sont absents et dont 4 biominerales sont au moins presents

"""
from interM_clusterExtraction import read_mcl_output
import os
#
    # """ lire chaque proteine dans chaque proteome de NoBiomes et des biominerales
    # trouver les clusters avec au moins 4 proteines biominerales =!
    # Verifier s'il y'a pas plus de 2 copies de proteines des NoBiomes
    # """
def read_fasta(directory_fasta, liste_fasta):
    dico_all = {}
    for filename in liste_fasta:
        if filename.endswith('.fasta'):
            path_proteome = os.path.join(directory_fasta, filename)
            proteome_name = filename.split('.')[0]
            dico_all.setdefault(proteome_name, [])
            for line in open(path_proteome):
                line = line.strip()
                if line.startswith('>'):
                    name_protein = line[1:].split()[0]
                    dico_all[proteome_name].append(name_protein)
    # for x in dico_all:
    #     print(x)
    return dico_all

def protein_biomineral(dico_all, list_biominerale, list_non_biominerale):
    '''Return 2 dico with all proteins in the biominerales and in the non - biominerales species
    '''
    dico_biominerale = dico_non_biominerale = {}
    for proteome in dico_all:
        if proteome in list_biominerale:
            dico_biominerale[proteome] = dico_all[proteome]
        elif proteome in list_non_biominerale:
            dico_non_biominerale[proteome] = dico_all[proteome]

    return dico_biominerale, dico_non_biominerale

def remove_cluster_nonBiom(group, dico_biominerale, dico_non_biominerale):
    """Count for each cluster(= a line in group file), the number of protein in biomineral and no biomineral sp
    return cluster with at least 4 biomineral proteins and 0 or one no biomineral
    """
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

                for proteome in dico_biominerale:
                    list_biomin = dico_biominerale[proteome]
                    if protein in list_biomin:
                        nb_prot_biom+=1

                for proteom in dico_non_biominerale:
                    list_non_biomin = dico_non_biominerale[proteom]
                    if protein in list_non_biomin:
                        nb_prot_no_biom+=1

            if nb_prot_biom > 4 and nb_prot_no_biom <= 1:
                # norm_biom = nb_prot_biom/nb_prot_total
                # norm_no_biom = nb_prot_no_biom/nb_prot_total
                # if norm_biom > norm_no_biom:

                print(line)

        #     print('biominerales=', nb_prot_biom)
        #     print('non biomin  =', nb_prot_no_biom)
        print('number_cluster', number_cluster)
if __name__ == '__main__':

    list_biominerale = ['Synechococcus_sp_PCC_6312', 'Synechococcus_calcipolaris', 'Thermosynechococcus_elongatus_BP1', 'Gloeomargarita_lithophora', 'Cyanothece_sp_PCC_7425', 'Chroococcidiopsis_thermalis_PCC_7203']
    list_non_biominerale = ['Synechocystis_sp_PCC6803', 'Nostoc_sp_PCC_7107', 'Gloeocapsa_sp_PCC_73106', 'Synechococcus_elongatus_PCC7942', 'Acaryochloris_marina_MBIC11017']
    groups = '/home/issa/Documents/stage/orthomcl/Intermediaires_clusters/wo_para/Inter_Cluster_sansDupli.txt'
    # proteomes_names =


    directory_fasta = '/home/issa/Documents/path2/selected_proteomes'
    liste_fasta  = os.listdir(directory_fasta)
    protein_name = read_fasta(directory_fasta, liste_fasta)
    biom, non_biom = protein_biomineral(protein_name, list_biominerale, list_non_biominerale)
    remove_cluster_nonBiom(groups, biom, non_biom)
    # dico_cluster = read_mcl_output(groups, proteomes_names)
