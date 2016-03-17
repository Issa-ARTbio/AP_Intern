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
            amino_acids = []
            if line.startswith('>'):
                header = line[1:].strip().split()[0]
            else:
                sequence = line.strip()
                seq_lenght = len(sequence)
                for aa in sequence:
                    amino_acids.append(aa)
            dico_fasta[header] = amino_acids
        # for x in dico_fasta:
        #     print(x,len(dico_fasta[x]), dico_fasta[x])
    return dico_fasta


def read_pyHCA_outF (hca_in):
    """return domains match in sequences by pyHCA program
    """

    sequence_lenghts = defaultdict(list) # to delete if not use
    domain_list = defaultdict(list)
    with open(hca_in, 'r') as hca:
        for line in hca:
            domain_line = ()
            if line.startswith('>'):
                header = line[1:].strip().split()[0]
                seq_lenght = line[1:].strip().split()[1]
                sequence_lenghts[header].append(seq_lenght)
            if line.startswith('domain'):
                domain_start  = line.strip().split()[1]
                domain_end    = line.strip().split()[2]
                domain_lenght = int(domain_end) - int(domain_start)
                domain_line   = (domain_start, domain_end, domain_lenght)
                domain_list[header].append(domain_line)
    for x in domain_list:
        print(x,len(domain_list[x]), domain_list[x])
        # print(dico_domain)
        # print(sequence_lenghts)
    return domain_list


# def match_domain (dico_fasta, domain_list):
#     """
#     """


if __name__ == '__main__':

    fasta_in = 'fasta_test.txt'
    hca_in   = 'test_hca.txt'

    read_fasta(fasta_in)
    read_pyHCA_outF(hca_in)
