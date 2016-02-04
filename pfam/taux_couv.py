#!/usr/bin/env python
# -*- coding: utf-8 -*-


''' Pour chaque proteome: calcul le taux de couverture pour chaque proteome
'''
import os, sys, argparse,re


def read_fasta(path_proteome):
    ''''Description de la fonction : lit un fichier .fasta et renvoit un dict avec key = protein_name et val = sequence
    '''
    dict_proteome = {}
    with open(path_proteome) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                name_protein = line[1:].split()[0]

                if name_protein in dict_proteome:
                    print("Error, protein {} already exists".format(name_protein))
                    sys.exit(1)

                dict_proteome[name_protein] = 0
            else:
                len_seq = len(line.strip())
                # dict_proteome[name_protein] = dict_proteome.get(name_protein, 0) + len_seq
                dict_proteome[name_protein] += len_seq
    # print (dict_proteome)
    return dict_proteome


def read_pfam_outp (pfam_ouf):
    '''Description de la fonction : renvoit un dict pour :key = protein_id, value = [(start, stop)] du domaine trouve sur pfam
    :un protein id est present en fonction du nombre de domaines qu'il porte '''
    pfam_dict = {}
    with open(pfam_ouf) as pf:
        for line in pf:
            if line.strip() and line[0] != '#':
                line = line.rstrip()
                tmp = line.split()
                seq_id = tmp[0]
                start = int(tmp[1]) - 1
                stop  = int(tmp[2])
                pos = start, stop
                if seq_id not in pfam_dict:
                    pfam_dict[seq_id] = []
                pfam_dict[seq_id].append(pos)

    return pfam_dict


def long_couvert(pfam_dict):
    '''Description de la fonction: renvoit la longeur totale couvert par les domaines sur la proteine
    :permet de passer de la redondance des protein_id à un seul id par protein
    '''
    len_domains = {}
    for seq_id in pfam_dict:
        values = pfam_dict[seq_id]
        for start, stop in values:
            ld = stop - start
            len_domains[seq_id] = len_domains.get(seq_id, 0) + ld

    return len_domains

def taux_couvert(dict_proteome, len_domains):
    '''Description de la fonction: prend comme entrees deux dico avec les mm keys
    dans dict_proteome: values = len_fasta
    dans len_domains: values = len_domains
    renvoit:
    Taux de couvertures pour le proteomes
    nombre de protein ayant au moins un domaine sur le nombre total de protein
    '''
    nmb_domain1 = 0
    nmbr_prot = 0
    list_taux= []
    len_proteome_domain = 0
    len_proteome_fasta  = 0
    for seq_id in dict_proteome:
        nmbr_prot =nmbr_prot+ 1
        len_fasta = dict_proteome[seq_id]
        len_proteome_fasta = len_proteome_fasta + len_fasta
        if seq_id in len_domains:
            print(seq_id)
            nmb_domain1 = nmb_domain1+1
            list_taux.append(seq_id)
            ld = len_domains[seq_id]
            len_proteome_domain = len_proteome_domain + ld
    taux = (len_proteome_domain/len_proteome_fasta)*100
    print(taux)
    print(nmb_domain1,'sur',nmbr_prot)
    return taux, nmb_domain1, nmbr_prot
        # else:
        #     print( seq_id, 'No domain found in Pfam\n')

# def taux_couvert_proteome():

if __name__ == '__main__':
    path_proteome = 'example.fasta'
    fasta = read_fasta (path_proteome)
    pfam_output = 'Pfam_example.pfam'
    table_pfam_dict = read_pfam_outp(pfam_output)
    len_domainsTot = long_couvert(table_pfam_dict)
    taux = taux_couvert(fasta, len_domainsTot)
