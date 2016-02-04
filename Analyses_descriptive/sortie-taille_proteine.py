#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" renvoit un fichier de sortie'.dat' avec 3 colonnes( Nom du Proteome - Nom Protein - Len Protein)
"""

import os, sys, argparse
from Box_plot_P20 import read_fasta

def protein_name_lenght (proteome_dir, out_file):
    with open(out_file, 'w') as outf:
        # outf.write('P'
        for filename in list_proteome_name:
            if filename.endswith('.fasta'):
                path_proteome = os.path.join(proteome_dir, filename)
                filename = filename.split('.')[0]
                list_proteome_name.append(filename)
                dico_protein = read_fasta(path_proteome)

                for prot in dico_protein :
                    seq = dico_protein[prot]
                    lenght_seq = len(seq)
                    outf.write( filename + '\t' +  prot + '\t' + str(lenght_seq) + '\n')
                    # print (seq, lenght_seq)
proteome_dir = '/home/issa/Documents/STAGE/Data/Init_data/proteomes/'
list_proteome_name = os.listdir(proteome_dir)
out_file = 'proteome_protein_name_lenght.dat'

protein_name_lenght(proteome_dir,  out_file)
