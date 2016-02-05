#!/usr/bin/env python
# -*- coding: utf-8 -*-


''' Pour chaque proteome: calcul le taux de couverture pour chaque proteome
'''
import os, sys,re


def read_fasta(path_proteome):
    ''''Description de la fonction : lit un fichier .fasta et renvoit un dict avec key = protein_name et val = sequence
    '''
    proteome_dict = {}
    with open(path_proteome) as f:
        proteome_dict.setdefault(filename , {})
        protein_dict = {}
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                name_protein = line[1:].split()[0]

                if name_protein in protein_dict:
                    print("Error, protein {} already exists".format(name_protein))
                    sys.exit(1)

                protein_dict[name_protein] = 0
            else:
                len_seq = len(line.strip())
                protein_dict[name_protein] += len_seq
        proteome_dict[filename]=protein_dict
    return proteome_dict


def read_pfam_outp (path_pfam_ouf):
    '''Description de la fonction : renvoit un dict pour :key = protein_id, value = [(start, stop)] du domaine trouve sur pfam
    :un protein id est present en fonction du nombre de domaines qu'il porte '''
    all_pfam = {}
    with open(path_pfam_ouf) as pf:
        all_pfam.setdefault(filename, {})
        pfam_dict = {}
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
        all_pfam[filename] = pfam_dict
    return all_pfam

def long_couvert(all_pfam):
    '''Description de la fonction: renvoit la longeur totale couvert par les domaines sur la proteine
    :permet de passer de la redondance des protein_id à un seul id par protein
    '''
    all_domains = {}
    for filename_pfam in all_pfam:
        all_domains.setdefault(filename_pfam, {})
        len_domains = {}
        pfam_dict = all_pfam[filename_pfam]
        for seq_id in pfam_dict:
            values = pfam_dict[seq_id]
            for start, stop in values:
                ld = stop - start
                len_domains[seq_id] = len_domains.get(seq_id, 0) + ld
        all_domains[filename_pfam]= len_domains
    return all_domains

def taux_couvert(proteome_dict, all_domains):
    '''Description de la fonction: prend comme entrees deux dico avec les mm keys
    dans dict_proteome: values = len_fasta
    dans len_domains: values = len_domains
    renvoit:
    Taux de couvertures pour le proteomes
    nombre de protein ayant au moins un domaine sur le nombre total de protein
    '''
    taux_proteome_couvert = {}
    freq_1domain = {}
    for filename_pfam, filename_fasta in zip (all_domains, proteome_dict):
        protein_dict = proteome_dict[filename_fasta]
        len_domains = all_domains[filename_pfam]
        # print (filename_fasta, 'prot =',  'protein_dict')
        # print (filename_pfam, 'domain =' , 'len_domains', '\n=================================')
        #initialise les dicos pour le taux de couv de chaque proteome et la frequence du nbr de domaine = 1
        taux_proteome_couvert.setdefault(filename_fasta, 0)
        freq_1domain.setdefault(filename_pfam, 0)
        nmb_domain1 = 0
        nmbr_prot = 0
        len_proteome_domain = 0
        len_proteome_fasta  = 0

        for seq_id, dom_id in zip(protein_dict, len_domains):
            nmbr_prot =nmbr_prot+ 1
            len_fasta = protein_dict[seq_id]
            len_proteome_fasta = len_proteome_fasta + int(len_fasta)

            nmb_domain1 = nmb_domain1+1
            ld = len_domains[dom_id]
            len_proteome_domain = len_proteome_domain + int(ld)

        taux = (len_proteome_domain/len_proteome_fasta)*100
        freq_domain1 = nmb_domain1/nmbr_prot
        print (filename_fasta, 'prot =',  nmbr_prot)
        print (filename_pfam, 'domain =' , nmb_domain1, '\n=================================')

        taux_proteome_couvert[filename_fasta] = taux
        freq_1domain[filename_pfam] = freq_domain1
    return taux_proteome_couvert, freq_1domain


if __name__ == '__main__':

    directory = 'proteomes_test/'
    liste_fasta = os.listdir(directory)

    pfam_out = 'results_pfam/'
    list_pfam_out = os.listdir(pfam_out)

    for filename in liste_fasta:
        if filename.endswith('.fasta'):
            path_proteome = os.path.join(directory, filename)
            filename = filename.split('.')[0]
            fasta = read_fasta (path_proteome)

    for filename in list_pfam_out:
        if filename.endswith('.pfam'):
            path_pfam_ouf = os.path.join(pfam_out, filename)
            filename = filename.split('.')[0]

        table_pfam_dict = read_pfam_outp(path_pfam_ouf)
        len_domainsTot = long_couvert(table_pfam_dict)
        taux = taux_couvert(fasta, len_domainsTot)
