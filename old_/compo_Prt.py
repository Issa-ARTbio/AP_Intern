# -*- coding: utf-8 -*-

import os
import os.path
import sys
from path import path

def comp_proteome (file_in):

    # proteome_file = os.listdir(proteome_dir)
    # proteome_path = [os.path.join(proteome_dir, filename) for filename in  os.listdir(proteome_dir)]
    list_aa = ['A', 'R', 'N', 'D', 'C', 'Q', 'E','G', 'H', 'I', 'L', 'K', 'M',
    'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    # for file_in in proteome_file:
    #     if file_in.endswith(".fasta"):
    #         print (file_in)
    # with open(file_in,'r') as f:
    f= open (file_in,'r')
    Prot = {}
    for line in f:
        line = line[:-1]
        # print (line)
        if line [0] != '>':
            for aa in line:
                Prot[aa] = Prot.get(aa, 0) +1
    f.close()
    line =  proteome_dir + 'proteome_dir = '
    for aa in list_aa:
        line += aa + ': ' + str(Prot.get(aa, 0)) + ' '
    # print (line)
    return (Prot)

file_in = '/home/issa/Documents/STAGE/Data/Données_Initiales/Proteomes_test/Anabaena_sp_PCC_7108.fasta'
proteome_dir = '/home/issa/Documents/STAGE/Data/Données_Initiales/Proteomes_test/'

result = []
for prt_file in proteome_dir:
    f_prt = comp_proteome (prt_file)
    result.append(f_prt)
print (result)
