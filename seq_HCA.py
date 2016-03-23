#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This programm is for matching the found hydrophobic cluster on the sequence
"""

from collections import defaultdict

def read_fasta (fasta_in):
    """return a dict which key = protein_id and value = list of amino acid in the sequence
    """
    with open(fasta_in, 'r') as fasta:
        dico_fasta = {}
        for line in fasta:
            line = line.strip()
            amino_acids = []
            if line.startswith('>'):
                header = line[1:].split()[0]
                # 2- essaye quelque chose comme ca:
                # dico_fasta[header] = []
            else:
                sequence = line.strip()
                seq_lenght = len(sequence)
                for aa in sequence:
                    amino_acids.append(aa)
                    # 3- dico_fasta[header].append(aa)
            # 1- si ta sequence est sur plusieurs lignes tu vas perdre la composition en acide amine des ligne precedentes
            dico_fasta[header] = amino_acids
    return dico_fasta


def read_pyHCA_outF (hca_in):
    """return domains match in sequences by pyHCA program
    defaultdict(list) : key = protein_id , value = tpl(domain_start, domain_end, domain_lenght)
    """
    #il faut faire tres attention avec le defaultdict
    #si la cle n'existe pas il la cree
    # donc quand tu parcours un dictionnaire si tu ne testes pas bien que la clef est presente ou absente il ne te le dira pas
    # ce qui peux cacher des erreurs / problemes / bugs
    # prefere un dictionnaire normal
    sequence_lenghts = defaultdict(list) # to delete if not use
    domain_hca = defaultdict(list)
    with open(hca_in, 'r') as hca:
        for line in hca:
            line = line.strip()
            domain = ()
            if line.startswith('>'):
                #header, seq_length = line[1:].split()
                #seq_length = int(seq_length)
                header = line[1:].split()[0]
                seq_lenght = line[1:].split()[1]
                sequence_lenghts[header].append(seq_lenght)
            if line.startswith('domain'):
                #name, domain_start, domain_end = line.strip().split()
                #ATTENTION n'oublie pas d'enlever 1 a la position start
                #dans segHCA comme dans pfam comme dans CDD on ne compte pas a partir de 0 mais a partir de 1
                # pas le stop qui est inclusif
                #il faut donc enlever 1 pour la que la premiere position d'un resultat de HCA/PFAM/CDD match l'index 0 d'une sequence en python
                #domain_start = int(domain_start) - 1
                domain_start  = line.strip().split()[1]
                domain_end    = line.strip().split()[2]
                domain_lenght = int(domain_end) - int(domain_start)
                domain   = (domain_start, domain_end, domain_lenght)
                domain_hca[header].append(domain)
    return domain_hca

def read_cdd_outF (cdd_in):
    """return dimain match in sequences by CDD scan
    defaultdict(list) : key = protein_id , value = tpl(domain_start, domain_end, domain_lenght)
    """
    domain_cdd = defaultdict(set)
    with open(cdd_in, 'r') as cdd:
        for line in cdd:
            line = line.strip()
            domain = ()
            if line.startswith('Q#'):
                element = line.split()
                header = element[2][1:]
                domain_start = element[6]
                domain_end   = element[7]
                # ATTENTION index start -1 (pas le stop)
                domain_lenght= int(domain_end) - int(domain_start)
                domain = (domain_start, domain_end, domain_lenght)
                domain_cdd[header].add(domain)

    return domain_cdd


def match_domain (dico_fasta, domain_hca, domain_cdd):
    # questions a te poster a toi meme:
    # - qu est ce que tu veux recuperer a la fin de la fonction?
    # - pour faire quoi ensuite?
    # - quel format sera le plus pratique pour reutiliser la sortie?
    """find matchs in HCA and CDD domain_list
    """
    # trop de defaultdict, c'est inutile un dictionnaire normal est mieux
    specific_protein_hca, specific_protein_cdd = defaultdict(list), defaultdict(list)
    hca, cdd = defaultdict(list), defaultdict(list)

    orphan_domains, same_domains = defaultdict(list), defaultdict(list) # len(hca_domain) >= 10
    overlap = 100
    print('CDD _ Data =', domain_cdd)
    print('HCA _ Data =', domain_hca)
    
    # tu n'as pas besoin de refaire des dictionnaires, il faut juste que tu recuperes les clefs 
    # qui sont a la fois dans domain_hca et dans domain_cdd
    # prot_has_both_dom = set(domain_hca.keys()).intersection(set(domain_cdd.keys()))
    # ca t'evite deux boucles "for"
    for protein in list(domain_hca):
        position_hca = domain_hca[protein]
        number_domain_hca = len(position_hca)

        if protein not in domain_cdd:
            specific_protein_hca[protein].append(domain_hca[protein])
            mes= ('{prot} not have domain in CDD').format(prot=protein)
            # print('NOTE :', mes)
            # print('======================>')
        else:
            hca[protein].append(domain_hca[protein])

    for proteine in domain_cdd:
        position_cdd = domain_cdd[proteine]
        number_domain_cdd = len(position_cdd)

        if proteine not in domain_hca:
            specific_protein_cdd[proteine].append(domain_cdd[proteine])
            mes = ('{prot} not have domain in pyHCA').format(prot=proteine)
            # print('NOTE :', mes)
        else:
            cdd[proteine].append(domain_cdd[proteine])

    print('the domain specific in CDD', specific_protein_cdd)
    print('the domain specific in pyHCA', specific_protein_hca)
    for protein in hca:
        hca_pos_domain = hca[protein]
        cdd_pos_domain = cdd[protein]
        # pourquoi faire un zip?
        for (h, c) in zip(hca_pos_domain, cdd_pos_domain):
            # ca c'est bizarre:
            for start, stop , lenght in (dom for dom in h):
                position_domaine = (start, stop , lenght)
                if lenght >= 10:

                    for start_cdd, stop_cdd , lenght_cdd in (dom_cdd for dom_cdd in c):
                        dom = (start_cdd, stop_cdd , lenght_cdd)
                        # a verifier
                        if -(overlap) <= (int(start) - int(start_cdd)) <= overlap or -(overlap) <= (int(stop) - int(stop_cdd)) <= overlap:
                            # ma formule est:
                            # new_start = max(start, start_cdd)
                            # new_stop = min(stop, stop_cdd)
                            # if new_start < new_stop and (new_stop - new_start) >= authorized_overlap:
                            # ____## it's an overlap
                            diff = (lenght, lenght_cdd)
                            same_domains[protein] = diff
                            mes =  ('I find the same domain between the CDD psition {dom1} and the HCA position {dom2} in {prot}').format(dom1 = dom, dom2 = position_domaine, prot=protein)
                            print(mes)
                        else:
                            orphan_domains[protein].append(position_domaine)
                            orphan_domains_cdd[proteine].append(dom)
                            print('No shared')
    print('orphan domains = ', orphan_domain)
            # for dom in h:
            #     for i, pos1 in enumerate(dom):
            #         print(pos1)
            #         if i == len(h) - 1 and int(pos1) >= 10:
            #             print(pos1)
                #     for j, pos2 in enumerate(domain_c for domain_c in cdd):
                #         print(pos2)
                    #     if -50 <= int(pos1) - int(pos2) <= 50:
                    #         mes =  ('I find one same in {prot} domain').format(prot=protein)
                    #         # print(mes)
                    # if len(position_hca) > len(position_cdd):
                    #     difference_hca = len(position_hca) - len(position_cdd)
                    #     mes =('I found {diff} more domain(s) in pyHCA').format(diff=difference_hca)
                    #     # print(mes)
                    # else:
                    #     difference_cdd = len(position_cdd) - len(position_hca)
                    #     mes = ('I found {diff} more domain(s) in CDD').format(diff=difference_cdd)
                    #     print(mes)

if __name__ == '__main__':

    fasta_in = 'fasta_test.txt'
    hca_in   = 'test_hca.txt'
    cdd_in = 'cdd_test.txt'
    fasta = read_fasta(fasta_in)
    hca = read_pyHCA_outF(hca_in)
    cdd = read_cdd_outF(cdd_in)
    match_domain (fasta, hca, cdd)
