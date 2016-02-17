#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" renvoit un fichier de sortie '.fasta' avec seulement les noms de sequences formattes pour etre lisible sur orthomcl(BlastParser)

-i dir/*.fasta -o dir/*.fasta

usage
===============
python adjust_seq.py
"""

import os, sys, argparse
import subprocess


def adjust_name (proteome_dir, out_dir):
    ''' Description de la fonction : lit et creer en mmeme temps un dossier avec des fichiers fasta et renvoit dans le dossier de sortie les fasta noms des fichiers formattes
    '''
  for filename in list_proteome_name:
      if filename.endswith('.fasta'):
        path_proteome = os.path.join(proteome_dir, filename)
        filename = filename.split('.')[0]
        sub_id = filename.split('_')
        identifiant = sub_id[0][0]+'_'+sub_id[-1]
        out_file = os.path.join(out_dir, identifiant+'.fasta')

        with open(path_proteome) as proteome, open(out_file, u'w') as outf:

          for line in proteome:
            line.rstrip()

            if line.startswith('>'):
              line = line[1:]
              protein_name = line.split()[0]
              print (protein_name)
              line = '>'+identifiant+'|'+str(protein_name)+'\n'
              outf.write(line)

            else:
              outf.write(line)

proteome_dir = '/home/issa/Documents/STAGE/Data/Init_data/proteomes/'
list_proteome_name = os.listdir(proteome_dir)
out_dir = '/home/issa/Documents/STAGE/Data/Init_data/proteomes_orthomcl/'
if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

adjust_name(proteome_dir,  out_dir)
# subprocess.call(["echo", i], shell=True)
