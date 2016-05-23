#!/usr/bin/env python
# -*- coding: utf-8 -*-



import os
from collections import defaultdict

def read_fasta(path_proteome, proteome_name):

    dict_proteome = {}
    with open(path_proteome) as f:
        dict_proteome[proteome_name]=0
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                pass
            else:
                dict_proteome[proteome_name] += len(line.strip())
    return (dict_proteome)

def read_fasta_DO (file_in, proteome_name):
    dict_orphan = {}
    with open(file_in) as orp:
        dict_orphan[proteome_name] = 0
        for line in orp:
            line = line.strip()
            if line.startswith('>'):
                pass
            else:
                dict_orphan[proteome_name]+=len(line.strip())
    return dict_orphan


def main():
    '''pour chaque proteome, retourne la taille du proteome, la taille de tous les DOs et le taux de couverture en résidues
    '''

    for filename in list_proteome:
        if filename.endswith('.fasta'):
            path_proteome = os.path.join(directory, filename)
            file_in = os.path.join(directory, filename+'.pos')
            proteome_name = filename.split('.')[0]

            fasta = read_fasta(path_proteome, proteome_name)
            do = read_fasta_DO(file_in, proteome_name)
            for proteome in do:
                taux = (do[proteome]/fasta[proteome])
                print(proteome, fasta[proteome], do[proteome], taux)

if __name__ == '__main__':

    directory = '/home/issa/Documents/stage/CDD_pyHCA/references/sortie_Cdd/'
    list_proteome = os.listdir(directory)
    main()
