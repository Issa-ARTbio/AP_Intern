# -*- coding: utf-8 -*-

import os
import os.path
import sys
from path import path
# aa = {'Ala' : 'A', 'aRg': 'R', 'asN': 'N', 'asp': 'D', 'Cys': 'C', 'Gln': 'Q', 'Gly': 'E', 'Gly' : 'G' , 'His': 'H', 'Ile' : 'I' , 'Leu' : 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro' : 'P', 'Pyl' : 'O', 'Sec':'U', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'};

def dir_path (l):

    fp = []
    for root, directories, files in os.walk(l):
        # print ('root---------=', root,'Directory _______= ', directories, 'files ##########=', files)
        for filename in files:
            filepath = os.path.join(root, filename)
            fp.append(filepath)
    for file_in in fp:
        if file_in.endswith(".fasta"):
            f = open (file_in , 'r')
            print(f)
            f.close()


l = '/home/issa/Documents/STAGE/Data/Données_Initiales/proteomes/'
dir_path (l)
