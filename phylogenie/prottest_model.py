#!/usr/bin/env python
# -*- coding: utf-8 -*-


'''Extraire le model de sortie de prottest
'''

import sys, os, re

def read_file (prottest_outf, model_retrieve):


    dico = {}
    with open(prottest_outf, 'r') as f, open (model_retrieve, 'w') as out_model:
        lines = f.read().splitlines()
        for i in range(len(lines)):
            line = lines[i]
            if line.startswith(my_search):
                out_model.write (file_name)
                tmp = line.split()
                model = tmp[-1]

                confidence = lines[i+1]
                t = confidence.split()
                CI = t[-1]
                dico[file_name] = str(model)
            # out_model.write ('|'+str(model))
                # print (file_name, '|', model, '|', CI)
    cluster_mod = []
    for cluster in dico:

        if dico[cluster] == 'LG':
            cluster_mod.append(cluster)
            print(cluster)
            # DAYHOFF, DCMUT, JTT, MTREV, WAG, RTREV, CPREV, VT, BLOSUM62, MTMAM, LG, MTART, MTZOA, PMB, HIVB, HIVW, JTTDCMUT, FLU, GTR

    # print(len(cluster_mod))


if __name__ == '__main__':

    prottest_outf = '/home/issa/Documents/stage/muscle/prottest/results/'
    liste_inp = os.listdir(prottest_outf)
    outf = 'sortie_model_prottest.txt'
    my_search = 'Best model according to LnL:'
    for file_name in liste_inp:
        if file_name.endswith('.prt'):
            protest = os.path.join(prottest_outf, file_name)
            file_name = file_name.split('.')[0]

        read_file (protest, outf)
