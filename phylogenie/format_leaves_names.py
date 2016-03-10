#!/usr/bin/env python
# -*- coding: utf-8 -*-

''' This program is for write in a file, the specie tree and one gene tree is the same leaves names formats.
Output files must be ready for HGT analyses by the tool HGT detection 3.2.
To see more : Boc, A., Philippe, H. and Makarenkov, V. (2010), Inferring and validating horizontal gene transfer events using bipartition dissimilarity, Systematic Biology.
'''

'''
usage
===========
python hgt_format.py

PARAMETERS:
Input:
directory_in contains newick format gene trees
specie_Tree: the reference specie tree to write first in each output
directory_out where output files will be store
'''
import re, os
from ete3 import Tree


def read_newick (file_in, directory):
    '''this function reads a -dir containing newick format and return a dict as key=name_file ans value= tree
    '''
    trees = []
    dico_tree = {}
    for filename in file_in:
        if filename.endswith('.tree'):
            path_tree = os.path.join(directory, filename)
            tree_cluster = filename.split('.')[1]

            tree = Tree(path_tree)
            for node in tree.traverse():
                if node.is_leaf():
                    name = node.name
                    if name.find("|") != -1:
                        split = name.split("|")
                        name  = split[0]
                        taxon = split[1]
                        node.add_feature("name", name)
            # print(tree_cluster, tree.write())
            dico_tree[tree_cluster] = tree.write()

    return dico_tree

def write_tree_file (specie_Tree, dico_tree, outf):
    '''For each entrie of the dict, writting in output format an file with the specie tree and the corresponding tree in the dict.
    the name of the output is the entrie of the dict
    '''
    with open(specie_Tree) as sp_tree:
        for line in sp_tree:
            line = line.strip()
        for cluster in dico_tree:
            tree = dico_tree[cluster]

            with open(outf+cluster+'.hgt', 'w') as outTree:
                outTree.write(line+'\n')
                outTree.write(tree+'\n')
        print (cluster, 'is written')




if __name__ == '__main__':

    directory_in = '/home/issa/Documents/stage/raxml/clusters_Gblock/bestTree/'
    file_in   =  os.listdir(directory_in)
    specie_Tree = '/home/issa/Documents/stage/raxml/specie_tree_Gblock/RAxML_bestTree.specie_TREE.tree'
    directory_out = '/home/issa/Documents/stage/hgt/Gblock_cluster_hgt/'
    my_tree = read_newick(file_in, directory_in)
    w = write_tree_file(specie_Tree, my_tree, directory_out)
