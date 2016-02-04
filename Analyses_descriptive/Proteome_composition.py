#!/usr/bin/env python
# -*- coding: utf-8 -*-
''' Pour chaque proteome dans le dossier, renvoie le nombre et pourcentage des acides aminés

usage:
python AP_Intern/Proteome_comp.py  -i proteomes/*.fasta -o composition_proteomes.dat
'''

# on peut tout importer sur une ligne (de façon raisonnable)
import os, sys, argparse

list_aa = ['A', 'R', 'N', 'D', 'C', 'Q', 'E','G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

def get_cmd():
    """ lit les arguments en ligne de commande
    """
    # initialise un parser d'arguments:
    # https://docs.python.org/3/library/argparse.html
    parser = argparse.ArgumentParser()
    # ajoute un argument
    parser.add_argument("-i", action="store", dest="inputfiles", help="list of proteomes", nargs="+")
    parser.add_argument("-o", action="store", dest="output", help="path to the output file")
    params = parser.parse_args()
    # print help if no arguments
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    # on accede aux parametres en faisant params.<nom du parametre>,
    # example: params.inputfiles retourne la liste des proteomes
    return params

def comp_proteome(path_proteome, outf):
    """ 'Description de la fonction' Lit un proteome et extrait sa composition en acide amines

    Parameters:
    ===========
    path_proteome: string
        the system path to the proteome in fasta format
    outf: file
        opened file to write results in
    """
    # en python on prefere reserver les noms de variables avec une majuscule pour les objects
    # Initialise la taille de chaque proteome
    taille_proteome = 0
    new_line = [os.path.basename(path_proteome)]
    composition = {}
    with open(path_proteome) as f:
        for line in f:
            # enleve les characters de fin de sequences
            line = line.strip()
            if not line.startswith('>'):
                for aa in line:
                    taille_proteome += 1
                    #Ajoute 0 si la clef est absente et 1 si elle est présente dans le dico composition
                    composition[aa] = composition.get(aa, 0) +1
        for aa in list_aa:
            # Ajoute le nombre total de l'aa s'il est presente, sinon, renvoie 0
            nb_aa = composition.get(aa, 0)
            percentage = nb_aa*100 / taille_proteome
            new_line.append("{} {}".format(nb_aa, percentage))
        outf.write("{}\n".format("\t".join(new_line)))
def main():
    """ start programme
    """
    # read command line
    params = get_cmd()
    with open(params.output, "w") as outf:
        # 'format' permet de formater un string, remplace '%s' '%f' dans la nouvelle syntax python
        outf.write("Proteome\t{}\n".format("\t".join(list_aa)))
        for filename in params.inputfiles:
            if filename.endswith(".fasta"):
                comp_proteome(filename, outf)



if __name__ == "__main__":
    main()
