#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os, sys
from collections import defaultdict


def read_cdd_outF (cdd_in, read_out):
    """summerize domain match in sequences by CDD scan in a output file
    """
    domain_cdd = defaultdict(set)
    with open(cdd_in, 'r') as cdd, open (read_out, 'w') as f:
        for line in cdd:
            line = line.strip()
            domain = ()
            if line.startswith('Q#'):
                element = line.split('\t')
                header = element[0].split()[2][1:]
                if header.find('['):
                    protein_id= header.split('[')[0][:-2].split()[0]
                # sp = element[0].split('[')[2][:-2]
                # if sp.find(']'):
                #     organism = sp.split(']')[0]
                # print (organism)
                domain_start = element[3]
                domain_start = int(domain_start) - 1
                domain_end   = int(element[4])
                acc = str(element[7])
                short_name = str(element[8])
                definition = str(element[11])
                domain_lenght= domain_end - domain_start

                f.write('>'+'\t'+str(header)+'|'+str(protein_id)+'\t'+str(domain_start)+'\t'+str(domain_end)+'\t'+acc+'\t'+short_name+'\t'+definition+'\n')



if __name__ == '__main__':

    #
    # cdd_in = '/home/issa/Documents/stage/cluster_biom6/Glycine_zipper/run_0_CS.txt.cdd'
    # read_out = '/home/issa/Documents/stage/cluster_biom6/Glycine_zipper/parser_run_0_CS_cdd.txt'
    # read_cdd_outF(cdd_in, read_out)


    directory = '/home/issa/Documents/stage/cd-hit/biominerales_clusters/CDD_des_7CB'
    list_files_in = os.listdir(directory)

    for cluster in list_files_in:
        if cluster.endswith('.cdd'):
            cdd_in = os.path.join(directory, cluster)
            read_out = os.path.join(directory, cluster+'.re')
            read_cdd_outF(cdd_in, read_out)
