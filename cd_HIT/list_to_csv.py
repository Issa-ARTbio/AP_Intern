#!/usr/bin/env python
# -*- coding: utf-8 -*-


import csv
import pyexcel as pe

def read_listDO(file_in, file_out, hmm):


    proteome=set()
    with open(file_in) as f, open(file_out, 'w') as outf:
        cwriter = csv.writer(outf, delimiter=' ')
        n = ''
        for profil in hmm:
            n+=profil+' '
        outf.write(n+'\n')
        for line in f:
            tmp_list = ''
            if line.startswith('**********'):
                proteome_name=line.split()[1]
                proteome.add(proteome_name)
                # print(proteome_name)
                outf.write(proteome_name+'\n')
                # if line.split()[0] == (x for x in hmm):
                #     # profil_hmm = line.split()[0]
                #     nb_hit = line.split()[1]
                #     tmp_list.append(nb_hit)
                #     print(profil_hmm, nb_hit)
                # for x in tmp_list:
                #     outf.write(x, end=" ")
                # else:
                #     outf.write('\n')







if __name__ == '__main__':

    hmm = ['Glycine_zipper', 'HMA_domain' , 'Biomin_2822_DO1' , 'Biomin_2822_DO2' , 'Biomin_12450', 'Biomin_16547' , 'Biomin_22533' , 'Biomin_22544', 'Biomin_32587']
    file_in = '/home/issa/Documents/stage/cd-hit/tmp/liste_DO.txt'
    file_out = '/home/issa/Documents/stage/cd-hit/tmp/tableau.csv'
    read_listDO(file_in, file_out, hmm)
