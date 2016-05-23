#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os, csv
from collections import defaultdict


def read_file(path_proteome):
    '''Read
    '''

    dict_count = defaultdict(list)
    with open(path_proteome) as f:
        for line in f:
            if line.startswith('#'):
                pass
            else:
                line = line.strip()
                element = line.split()
                cluster_id = element[0]
                cluster = cluster_id.split('.')[0]
                protein_id=element[3]
                if protein_id not in dict_count[cluster]:
                    dict_count[cluster].append(protein_id)


    return dict_count


def write_csv(directory, list_prot, output):

    with open (output+str('DO_tab.csv'), 'w') as outf:
        for filename in list_prot:
            if filename.endswith('.out'):
                path_proteome = os.path.join(directory, filename)
                proteome_name = filename.split('.')[0]
                dico = read_file(path_proteome)
                outf.write(str(proteome_name)+'\t')
                print('**********',proteome_name,'**********')
                for x in dico:
                    l = len(dico[x])
                    outf.write((str(x)+'\t'+str(l)+'\t'))
                outf.write('\n')
                print ('>', proteome_name , x, l)
if __name__ == '__main__':
    # main(output)
    directory = '/home/issa/Documents/stage/cd-hit/tmp/results_hmmscan/'
    list_prot = os.listdir(directory)

    output = '/home/issa/Documents/stage/cd-hit/tmp/'
    write_csv(directory, list_prot, output)
