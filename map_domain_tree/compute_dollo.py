#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""

"""

import os
from collections import defaultdict
import ete3
from ete3 import PhyloTree, TreeStyle, Tree

def read_cdd_outF (cdd_in):
    """return domain match in sequences by CDD scan
    dict : key = protein_id , value(list of tuples) = tpl(domain_start, domain_end, domain_lenght)
    """
    domain_cdd = defaultdict(set)
    with open(cdd_in, 'r') as cdd:
        for line in cdd:
            line = line.strip()
            domain = ()
            if line.startswith('Q#'):
                element = line.split('\t')
                header = element[0].split()[2][1:]
                if header.find('['):
                    header= header.split('[')[0]
                domain_start = element[3]
                domain_start = int(domain_start) - 1
                domain_end   = int(element[4])
                domain_lenght= domain_end - domain_start
                domain = (domain_start, domain_end, domain_lenght)
                domain_cdd[header].add(domain)

    return domain_cdd


def assign_domain_tree (sp_tree, domain_cdd):

    tree = Tree(sp_tree)
    for node in tree.traverse():


if __name__ == '__main__':

    dir_dom_cdd = 'directory for cdd domains'
    sp_tree = '/home/issa/Documents/stage/raxml/specie_tree_Trimal/RAxML_bestTree.specie_TREE_trimal.tree'
