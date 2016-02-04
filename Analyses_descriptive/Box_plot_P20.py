#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Renvoit un plot de la taille des proteomes ( 20 proteomes au max)

-i dir/*.fasta -o plot.png

"""


import os, sys, argparse,re
import matplotlib
import matplotlib.pyplot as plt

def read_fasta(path_proteome):
    ''''Description de la fonction' : lit un fichier .fasta et renvoit un dict avec key = protein_name et val = sequence
    '''
    dict_proteome = {}
    list_protein_name = []
    with open(path_proteome) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                name_protein = line[1:].split()[0]
                # list_protein_name.append(name_protein)
                if name_protein in dict_proteome:
                    print("Error, protein {} already exists".format(name_protein))
                    sys.exit(1)
                dict_proteome[name_protein] = ""
            else:
                dict_proteome[name_protein] += line.strip()

    return (dict_proteome)

def compute_length (dict_proteome):
    ''' 'Description de la fonction':  a partir du dict, renvoit la liste de la longeur de chaque sequence
    '''
    list_protein_lenght  = []
    for protein_name in dict_proteome.keys():
        len_prot = len(dict_proteome [protein_name])
        list_protein_lenght.append (len_prot)
    return list_protein_lenght

if __name__ == '__main__': #Execute ce qui est en dessous que quand le programme est execute. Importe plus facilement les fonctions

    proteome_dir = '/home/issa/Documents/STAGE/Init_data/proteomes/'
    out_file = 'outf.dat'
    l = os.listdir(proteome_dir)

    #Lecture du dossier contenant les proteomes
    liste_Totale = []
    list_proteome_name = []
    for filename in l:
        if filename.endswith('.fasta'):
            path_proteome = os.path.join(proteome_dir, filename)
            filename = filename.split('.')[0]
            list_proteome_name.append(filename)
            dico_protein = read_fasta(path_proteome)
            long_protein = compute_length (dico_protein)
            liste_Totale.append(long_protein)

    # Construction d'un plot avec 20 listes de valeurs au max

    #taille max d'un plot = 20
    for i in range(0, len(liste_Totale), 20):
        data20 = liste_Totale[i : i +20]
        name_20 = list_proteome_name[i : i +20]
        # Construction du plot
        fig, ax = plt.subplots()
        plt.boxplot(data20, 0, 'gD')
        ax.set_xticklabels (name_20, rotation='vertical')
        for i in ax.set_xticklabels(name_20, rotation= 'vertical'):
            if re.search('Synechococcus_sp_PCC_6312', str(i)) :
                i.set_color('green')
            elif re.search('Synechococcus_calcipolaris' , str(i)):
                i.set_color('green')
            elif re.search('Thermosynechococcus_elongatus_BP1' , str(i)):
                i.set_color('green')


            elif re.search('Gloeomargarita_lithophora' , str(i)):
                i.set_color('red')
            elif  re.search ('Cyanothece_sp_PCC_7425', str(i)):
                i.set_color('red')
            elif  re.search ('Chroococcidiopsis_thermalis_PCC_7203', str(i)):
                i.set_color('red')
        plt.grid(True)
        plt.ylabel('Sequences', fontsize=15, color='blue')
        plt.xlabel('Proteomes', fontsize=15, color='red')
        plt.title(u"Distribution de la Taille des sequences dans les Proteomes", fontsize=17, fontdict={'family': 'monospace'})
        # plt.savefig(filename+u'boxplotProteome.png')
        plt.show()
