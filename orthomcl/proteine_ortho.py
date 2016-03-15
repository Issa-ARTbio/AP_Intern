#!/usr/bin/env python
# -*- coding: utf-8 -*-






def read_proteinOrtho_output(file_in):

    with open(file_in, 'r') as PO:
        nb=0
        for line in PO:
            line = line.strip()
            nb+=1
        print(nb)





if __name__ == '__main__':
    file_in = '/home/issa/Documents/stage/proteinOrtho/MyresultPO.proteinortho'
    read_proteinOrtho_output(file_in)
