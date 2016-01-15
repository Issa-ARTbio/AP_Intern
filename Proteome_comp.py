# -*- coding: utf-8 -*-

import os
import os.path
import sys

''' Pour chaque protome dans le dossier, renvoie le nombre et et pourcentage des acides aminés

'''

# aa = {'Ala' : 'A', 'aRg': 'R', 'asN': 'N', 'asp': 'D', 'Cys': 'C', 'Gln': 'Q', 'Gly': 'E', 'Gly' : 'G' , 'His': 'H', 'Ile' : 'I' , 'Leu' : 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro' : 'P', 'Pyl' : 'O', 'Sec':'U', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'};

import os
import os.path
import sys
from path import path

def comp_proteome (proteome_dir, list_aa):
    Taille_proteome = 0 # Initialise la taille de chaque protéome
    # output = open(output_file, 'w')
    with open(proteome_dir) as f:
        Prot = {}
        for line in f:
            line = line[:-1]
            if line [0] != '>':
                for aa in line:
                    Taille_proteome= Taille_proteome+1
                    Prot[aa] = Prot.get(aa, 0) +1 #Ajoute 0 si la clé est absente et 1 si elle est présente dans le dico Prot
        line = os.path.basename(proteome_dir)+ '\t'
        for aa in list_aa:
            line += str(Prot.get(aa, 0)) # Ajoute le nombre total de l'aa s'il est présente 0
            percentage = Prot.get(aa, 0)*100 / Taille_proteome
            line = line + ' ' + str(percentage) + '\t'
        print (line)
        # output.write(line)
        # output.close()

# output_file = '/home/issa/Documents/STAGE/Data/Données_Initiales/Proteomes_test/'.
proteome_dir = '/home/issa/Documents/STAGE/Data/Données_Initiales/Proteomes_test/'
list_aa = ['A', 'R', 'N', 'D', 'C', 'Q', 'E','G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
print ( 'Proteome\t'+"\t".join(list_aa) )
for filename in os.listdir(proteome_dir):
    if filename.endswith(".fasta"):
        path_proteome = os.path.join(proteome_dir, filename)
        comp_proteome(path_proteome, list_aa)

# proteome_file = os.listdir(proteome_dir)
