#!/usr/bin/env python
# -*- coding: utf-8 -*-




''' Extract fasta sequences of proteins found from the hmmscan of ODs
'''

import os
from collections import defaultdict

def read_hhmscan_output(file_in):
    ''' read the output of hmmscan and return dict with the found protein_id
    '''
    dict_count = defaultdict(list)
    with open(file_in, 'r') as f:
        for line in f:
            if line.startswith('#'):
                pass
            else:
                line = line.strip()
                element = line.split()
                cluster_id = element[0]
                cluster = cluster_id.split('.')[0]
                protein_id=element[3]
                if protein_id not in dict_count[cluster]:
                    dict_count[cluster].append(protein_id)
    return dict_count

def read_fasta(path_proteome):

    ''''Description de la fonction : lit un fichier .fasta et renvoit un dict avec key = proteome_name et val = dict('protein_name': sequence)
    '''
    with open(path_proteome) as f:
        protein_dict = {}
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                name_protein = line[1:].split()[0]

                if name_protein in protein_dict:
                    print("Error, protein {} already exists".format(name_protein))
                    sys.exit(1)

                protein_dict[name_protein] = ''
            else:
                protein_dict[name_protein] += line.strip()
    return protein_dict
def match_sequence(directory, list_prot, output):

    for filename in list_prot:
        if filename.endswith('.fasta'):
            path_proteome = os.path.join(directory, filename)
            file_in = os.path.join(directory, filename+".out")
            proteome_name = filename.split('.')[0]

            hmm = read_hhmscan_output(file_in)
            fasta = read_fasta(path_proteome)

            with open(output+str(proteome_name)+'.do.fasta', 'w') as outf:
                outf.write(str(proteome_name)+'\n')
                for cluster in hmm:
                    proteins_id = hmm[cluster]
                    outf.write('****** '+str(cluster)+'\t'+len(proteins_id)+'\n')
                    for protein in proteins_id:
                        if protein in fasta:
                            seq = fasta[protein]
                            outf.write('>'+str(protein)+'\n'+str(seq)+'\n')
                            print(protein)
                            print(seq)
                    outf.write('\n')

if __name__ == '__main__':
    # main(output)
    directory = '/home/issa/Documents/stage/cd-hit/biominerales_clusters/sequences_DO_Biomin/'
    list_prot = os.listdir(directory)
    output = '/home/issa/Documents/stage/cd-hit/biominerales_clusters/sequences_DO_Biomin/output/'
    match_sequence(directory, list_prot, output)
