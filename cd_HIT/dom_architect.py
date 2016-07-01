#!/usr/bin/env python
# -*- coding: utf-8 -*-

''' Annotation des proteines sur les sp biomines trouves a partir du hmmscan des profils des DOs
'''




import os
from collections import defaultdict


# def read_hhmscan_output(file_in):
#     ''' read the output of hmmscan and return dict with the found protein_id
#     '''
#     dict_count = defaultdict(list)
#     with open(file_in, 'r') as f:
#         for line in f:
#             if line.startswith('#'):
#                 pass
#             else:
#                 line = line.strip()
#                 element = line.split()
#                 cluster_id = element[0]
#                 cluster = cluster_id.split('.')[0]
#                 protein_id=element[3]
#                 hmm_len = element[2]
#                 protein_len = element[5]
#                 hmm_start_pos = element[15]
#                 alignment_start_pos = element[17]
#                 alignment_stop_pos = element[18]
#                 tpl = (protein_id, protein_len, alignment_start_pos, alignment_stop_pos)
#                 if tpl not in dict_count[cluster]:
#                     dict_count[cluster].append(tpl)
#     # print(dict_count)
#     return dict_count


def read_cdd_outF (cdd_in):
    """summerize domain sequences match in output file of CDD scan.
    """
    annot_cdd = defaultdict(set)
    with open(cdd_in, 'r') as g :
        for line in g:
            line = line.strip()
            domain = ()
            if line.startswith('Q#'):
                element = line.split('\t')
                header = element[0].split()[2][1:]
                sp = element[0].split()[4]
                # if sp.find('['):
                #     org = sp.split('[')[1][:-1]
                definition = str(element[11])
                annot_cdd[header] = definition
                print(header, '\t', sp, '\t', definition)
    return annot_cdd

# def annotation(directory, list_prot, output):
#     '''
#     '''
#     with open(output+str('table_dom_archi')+'.txt', 'w') as outf:
#         for filename in list_prot:
#             if filename.endswith('.fasta'):
#                 cdd_in = os.path.join(directory, filename+'.cdd')
#                 file_in = os.path.join(directory, filename+'.out')
#                 proteome_name = filename.split('.')[0]
#
#                 hmm_matches = read_hhmscan_output(file_in)
#                 cdd = read_cdd_outF(cdd_in)
#                 # print(proteome_name, hmm_matches)
#
#                 for cluster in hmm_matches:
#                     # outf.write(str(cluster)+'\n'+str(proteome_name)+'\n')
#
#                     nmb_prot = hmm_matches[cluster]
#                     for protein_id, protein_len, start_pos, stop_pos in nmb_prot:
#                         if protein_id in cdd:
#                             annot = cdd[protein_id]
#                             print(cluster, proteome_name, protein_id)
#                             print(protein_id, protein_len, start_pos, stop_pos, annot)
#                             outf.write(cluster+'\t'+proteome_name+'\t'+protein_id +'\t'+ protein_len +'\t'+ str(start_pos)+'\t'+ str(stop_pos)+'\t'+ annot+'\n')
#                             print('ANNOT: ', annot)
#                             print('***************************')





if __name__ == '__main__':
    cdd_in = '/home/issa/Documents/stage/cluster_biom6/cluster_33/HMA/3D/archtec.cdd'
    read_cdd_outF(cdd_in)
    # directory = '/home/issa/Documents/stage/cd-hit/dom_architec/'
    # liste_files = os.listdir(directory)
    # output = '/home/issa/Documents/stage/cd-hit/dom_architec_output/'
    # annotation(directory, liste_files, output)
