#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os

def concate_fasta(dir_in, fasta_files, outputF):

    """make ONE fasta file from the concatenation of all fasta in the directory

    """
    with open(outputF, 'w') as outF:
        for filename in fasta_files:
            if filename.endswith('.orp'):
                fasta_in = os.path.join(dir_in, filename)
                proteome_name = filename.split('.')[0]
                with open (fasta_in) as fasta:
                    for line in fasta:
                        line = line.strip()
                        if line.startswith('>'):
                            proteine_name = str(line)
                            outF.write(proteine_name+'\n')
                        else:
                            outF.write(line+'\n')

def format_header(file_out, fasta_in):
    """Change the header in the header in the output file of CDD scan
    """

    #  for filename in fasta_files:
    #      if filename.endswith('.fasta'):
    #          fasta_in = os.path.join(dir_in, filename)
    #          file_out = os.path.join(dir_out,filename)
    #          proteome_name = filename.split('.')[0]
    with open(file_out, 'w') as outF, open (fasta_in) as fasta:
         for line in fasta:
             line = line.strip()
             if line.startswith('>'):
                 tmp = line.split('|')
                 protein_id = tmp[3]
                 sp_id = tmp[4].split('[')[1][:-1]
                 outF.write('>'+protein_id+'|'+str(sp_id)+'\n')
                 print(protein_id, sp_id)
             else:
                 seq = line
                 outF.write(line+'\n')
if __name__ == '__main__':

    dir_in = '/home/issa/Documents/stage/cd-hit/fasta/'
    fasta_files = os.listdir(dir_in)
    outputF = '/home/issa/Documents/stage/cd-hit/data2.fasta'
    concate_fasta(dir_in, fasta_files, outputF)
    # file_out = '/home/issa/Documents/stage/cluster_biom6/cluster_33/Glycine_zipper/renamed_re_run_0_AS.fasta.aln'
    # fasta_in = '/home/issa/Documents/stage/cluster_biom6/cluster_33/Glycine_zipper/re_run_0_AS.fasta.aln'
    # format_header (file_out, fasta_in)
