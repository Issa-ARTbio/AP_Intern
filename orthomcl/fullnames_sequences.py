#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''remplacer les noms des sequences formattes par leurs noms en entiers
'''




import sys, os, re

def read_All_fasta(directory, liste_fasta):
    '''extract fasta sequence corresponding to the names in core file
    '''
    dict_proteome = {}
    list_protein_name = []
    for filename in liste_fasta:
        if filename.endswith('.fasta'):
            path_proteome = os.path.join(directory, filename)
            proteome_name = filename.split('.')[0]
            for line in open(path_proteome):
                line = line.strip()
                if line.startswith('>'):
                    name = line[1:]
                    tmp = name.split()
                    name_protein = tmp[0]
                    full_proteome_name = tmp[1]
                    seq_name = full_proteome_name+'|'+name_protein
                    if seq_name in dict_proteome:
                        print("Error, protein {} already exists".format(name_protein))
                        sys.exit(1)
                    dict_proteome[seq_name] = ""
                else:
                    dict_proteome[seq_name] += line.strip()

    return dict_proteome

def give_full_name(liste_fasta, proteomes_name):
    ''''Description de la fonction' : retransforme les noms des Proteomes (sous format lisible par orthomcl) aux noms complets initiaux
    '''

    # proteomes_name= [thermalis.replace('C_thermalis7203','C_7203') for thermalis in proteomes_name]
    match_name={}
    cnt = 0
    for filename in liste_fasta:
        if filename.endswith('.fasta'):
          filename = filename.split('.')[0]
          sub_id = filename.split('_')
          identifiant = sub_id[0][0]+'_'+sub_id[-1]
          match_name[identifiant]= filename
          cnt += 1
    match_name['C_thermalis7203']= match_name['C_7203']
    return match_name

def given_fullnames (cluster_out, cluster_in, all_fasta,match_name):

    with open(cluster_in, 'r') as cluster, open(cluster_out, 'w' ) as outfile:

        for line in cluster:
            line = line.strip()
            if line.startswith('>'):


                name_protein = line[1:].split('|')
                protein_name = name_protein[1]
                proteome_id  = name_protein[0]
            else:
                seq = line.strip()

                name = match_name[proteome_id]+"|"+protein_name
                tmpr = name.split('|')
                if name in all_fasta:
                    outfile.write('>'+name+'\n')
                    if all_fasta[name]== seq:
                        outfile.write(seq+'\n')
                    else:
                        print("Error : sequence is not corresponding for {} ".format(name))
                        print(name)
                        print(proteome_id, protein_name)
                        print(all_fasta[name])
                        print(seq)
                        sys.exit(1)
                else:
                    print("Error cannot find name", name)
                    sys.exit(1)
                    # Error cannot find name Microcoleus_vaginatus_FGP2|EGK89409.1
                    # Error cannot find name Fischerella_sp_JSC11|EHC11136.1
                    # Error cannot find name Cylindrospermopsis_raciborskii_CS505|EFA71456.1

        print("Fasta writing")

if __name__ == '__main__':
    directory = '/home/issa/Documents/stage/initial_data/proteomes/'
    liste_fasta = os.listdir(directory)
    all_fasta   = read_All_fasta(directory, liste_fasta)
    cluster_out = '/home/issa/Documents/stage/orthomcl/Intermediaires_clusters/wo_para/tmp/'
    cluster_dir = '/home/issa/Documents/stage/orthomcl/Intermediaires_clusters/wo_para/cluster/'
    clusters    =  os.listdir(cluster_dir)

    directory = '/home/issa/Documents/stage/orthomcl/proteomes_format_shortname/'
    liste_fa = os.listdir(directory)
    proteomes_name=[]
    for filename in liste_fa:
        if filename.endswith('.fasta'):
            proteome_name = filename.split('.')[0]
            proteomes_name.append(proteome_name)
    match = give_full_name(liste_fasta, proteomes_name)
    for filename in clusters:
            if filename.endswith('.fasta'):
                cluster_in = os.path.join(cluster_dir, filename)
                cluster_outf= os.path.join(cluster_out, filename)
                given_fullnames(cluster_outf, cluster_in, all_fasta, match)
    print('succefully done !!!')
