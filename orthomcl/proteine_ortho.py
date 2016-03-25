#!/usr/bin/env python
# -*- coding: utf-8 -*-


def read_proteinOrtho_output(file_in):

    with open(file_in, 'r') as PO:
        nb=0
        for line in PO:
            line = line.strip()
            # if line.startswith('>'):
            nb+=1
        print(nb)



if __name__ == '__main__':
    file_in = '/home/issa/Documents/stage/proteinOrtho/MyresultPO.proteinortho'
    # fasta = '/home/issa/Documents/stage/CDD_pyHCA/CDD/analyse/test_couverture/Cyanobacterium_YellowstoneA.fasta'
    read_proteinOrtho_output(fasta)
#Cyanobacterium_YellowstoneA 2760
