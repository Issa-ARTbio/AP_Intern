#!/usr/bin/env python
# -*- coding: utf-8 -*-

''' Annotation des proteines sur les sp biomines trouves a partir du hmmscan des profils des DOs
'''




import os
from collections import defaultdict


def read_hhmscan_output(file_in):
    ''' read the output of hmmscan and return dict with the found protein_id
    '''
    dict_count = defaultdict(list)
    with open(file_in, 'r') as f:
        for line in f:
            if line.startswith('#'):
                pass
            else:
                line = line.strip()
                element = line.split()
                cluster_id = element[0]
                cluster = cluster_id.split('.')[0]
                protein_id=element[3]
                hmm_len = element[2]
                protein_len = element[5]
                hmm_start_pos = element[15]
                alignment_start_pos = element[17]
                tpl = (protein_id, protein_len, )
                if protein_id not in dict_count[cluster]:
                    dict_count[cluster].append(protein_id)
    return dict_count


def read_cdd_outF (cdd_in):
    """summerize domain match in sequences by CDD scan in a output file
    """
    domain_cdd = defaultdict(set)
    with open(cdd_in, 'r') as g :
        for line in g:
            line = line.strip()
            domain = ()
            if line.startswith('Q#'):
                element = line.split('\t')
                header = element[0].split()[2][1:]

                if header.find('['):
                    protein_id= header.split('[')[0][:-2].split()[0]
                sp = element[0].split('[')[2][:-2]
                if sp.find(']'):
                    organism = sp.split(']')[0]
                print (organism)
                domain_start = element[3]
                domain_start = int(domain_start) - 1
                domain_end   = int(element[4])
                acc = str(element[7])
                short_name = str(element[8])
                definition = str(element[11])
                domain_lenght= domain_end - domain_start

if __name__ == '__main__':

    file_in = '/home/issa/Documents/stage/cd-hit/dom_architec/Synechococcus_sp_PCC_6312.fasta.out'
    cdd_in = '/home/issa/Documents/stage/cd-hit/dom_architec/Synechococcus_sp_PCC_6312.fasta.cdd'
    # directory = '/home/issa/Documents/stage/cd-hit/dom_architec'
    # liste_files = os.listdir(directory)
    hmm_matches = read_hhmscan_output(file_in)
    cdd = read_cdd_outF (cdd_in)
