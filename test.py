#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''Renvoit la FREQUENCE MINIMALE, MEDIUM et MAXIMALE pour chaque aa et le proteome name

usage:
=====
$ python3 Freq_AA.py -i(.dat)
'''

import os, sys, argparse, re
import matplotlib
import matplotlib.pyplot as plt
import numpy as np


def Freq_AA (proteom_composition_aa):


    dico_freq = dict()
    with open(proteom_composition_aa) as f:
        line = f.readline()
        line = line.strip()
        liste_aa = line.split('\t')[1:]
        for line in f:
            line = line.strip()
            tmp = line.split('\t')
            proteome_name = tmp[0]
            proteome_name = proteome_name.split('.')[0]
            values = tmp[1:]
            for aa in liste_aa:
              my_list = [aa]
              for proteome in proteome_name:
                for i in values:
                  freq = i.rstrip().split(' ')[1]
                  freq_min = min(freq)
                  freq_max = max(freq)
                  print (proteome_name,  freq_min)






proteom_composition_aa = '/home/issa/Documents/STAGE/Firsts_Analyses/Resultats/proteome_composition_AA.csv'
frequence_aa = Freq_AA (proteom_composition_aa)
