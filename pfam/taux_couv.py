#!/usr/bin/env python
# -*- coding: utf-8 -*-


''' Pour chaque proteome:
    calcul le taux de couverture des domaines
    et la frequence des proteines avec au moins 1 domaine present

Parametres:

-i dir (.fasta) dir (.pfam) -o []


usage
==============
python taux_couv.py
'''
import os, sys,re
import matplotlib
import matplotlib.pyplot as plt

def read_fasta(path_proteome):
    ''''Description de la fonction : lit un fichier .fasta et renvoit un dict avec key = proteome_name et val = dict('protein_name': len_sequence)
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
    '''Description de la fonction : renvoit un dict pour :key = proteome_name, value = dict('protein_name':[(start, stop)]) du domaine trouve sur pfam
    :NB: un protein_id est present en fonction du nombre de domaines qu'il porte '''
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
    :NB: permet de passer de la redondance des protein_id de la fonction precedente à un seul id par protein et comme value : la couverture total des domaine/protein
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
    '''Description de la fonction: prend comme entrees deux dico avec les mm keys = proteome name
    dans proteome_dict : values = dict {'protein_name': seq_len}
    dans all_domains: values = dict {'protein_name': len_domain}
    renvoit:
    Taux de couvertures pour les proteomes
    Et la frequence de protein(ayant au moins un domaine) sur le nombre total de protein
    '''
    taux_proteome_couvert = {}
    freq_1domain = {}
    for proteome_name in proteome_dict:
        #recupere les longeurs des sequences proteiques
        protein_dict = proteome_dict[proteome_name]
        taux_proteome_couvert.setdefault(proteome_name, 0)
        #recupere les longeurs des sequences des domaines
        len_domains = all_domains[proteome_name]
        freq_1domain.setdefault(proteome_name, 0)

        nmb_domain1 = 0
        nmbr_prot = 0
        len_proteome_domain = 0
        len_proteome_fasta  = 0
        #Calcul de la longueur totale pour chaque  proteome
        for seq_id in protein_dict:
            nmbr_prot =nmbr_prot+ 1
            len_fasta = protein_dict[seq_id]
            len_proteome_fasta = len_proteome_fasta + int(len_fasta)

            for seq_id in len_domains:
                nmb_domain1 = nmb_domain1+1
                ld = len_domains[seq_id]
                len_proteome_domain = len_proteome_domain + int(ld)
        print(len_proteome_fasta)

        taux = (len_proteome_domain/len_proteome_fasta)*100
        freq_dom1 = nmb_domain1/nmbr_prot
        taux_proteome_couvert[proteome_name] = taux
        freq_1domain[proteome_name] = freq_dom1
    return taux_proteome_couvert, freq_1domain


def plotting(taux_proteome_couvert, freq_1domain):
    '''renvoit un histogramme pour le taux de couverture de proteome
    Un histogramm pour la frequence de la presence d'un domaine
    '''
    list_de_taux = []
    list_freq_domain1 = []
    prot_name = []
    for proteome_name in taux_proteome_couvert:
        prot_name.append(prot_name)
        taux = taux_proteome_couvert[proteome_name]
        list_de_taux.append(taux)
    for proteome_name in freq_1domain:
        freq = freq_1domain[proteome_name]
        list_freq_domain1.append(freq)

if __name__ == '__main__':

    directory = 'proteomes_test/'
    liste_fasta = os.listdir(directory)

    pfam_out = 'results_pfam/'
    list_pfam_out = os.listdir(pfam_out)
    my_two_dict = [liste_fasta, list_pfam_out]
    for filename in liste_fasta:
        if filename.endswith('.fasta'):
            path_proteome = os.path.join(directory, filename)
            filename = filename.split('.')[0]

    for filename in list_pfam_out:
        if filename.endswith('.pfam'):
            path_pfam_ouf = os.path.join(pfam_out, filename)
            filename = filename.split('.')[0]

            fasta = read_fasta (path_proteome)
            table_pfam_dict = read_pfam_outp(path_pfam_ouf)
            len_domainsTot = long_couvert(table_pfam_dict)
            taux_prot, taux_dom = taux_couvert(fasta, len_domainsTot)
            plot_hist = plotting(taux_prot, taux_dom)
