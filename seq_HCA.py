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
            else:
                sequence = line.strip()
                seq_lenght = len(sequence)
                for aa in sequence:
                    amino_acids.append(aa)
            dico_fasta[header] = amino_acids
    return dico_fasta


def read_pyHCA_outF (hca_in):
    """return domains match in sequences by pyHCA program
    defaultdict(list) : key = protein_id , value = tpl(domain_start, domain_end, domain_lenght)
    """

    sequence_lenghts = defaultdict(list) # to delete if not use
    domain_hca = defaultdict(list)
    with open(hca_in, 'r') as hca:
        for line in hca:
            line = line.strip()
            domain = ()
            if line.startswith('>'):
                header = line[1:].split()[0]
                seq_lenght = line[1:].split()[1]
                sequence_lenghts[header].append(seq_lenght)
            if line.startswith('domain'):
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
                domain_lenght= int(domain_end) - int(domain_start)
                domain = (domain_start, domain_end, domain_lenght)
                domain_cdd[header].add(domain)
    for x in domain_cdd:
        print(x,len(domain_list[x]), domain_list[x])
    return domain_cdd


def match_domain (dico_fasta, domain_hca, domain_cdd):
    """find matchs in HCA and CDD domain_list
    
    """


if __name__ == '__main__':

    fasta_in = 'fasta_test.txt'
    hca_in   = 'test_hca.txt'
    cdd_in = 'cdd_test.txt'
    read_fasta(fasta_in)
    read_pyHCA_outF(hca_in)
    read_cdd_outF(cdd_in)
