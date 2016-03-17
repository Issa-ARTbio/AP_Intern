#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''change protein name to proteome name in the cluster trees
'''
import re, os
from ete3 import Tree
import random
def read_newick (path_tree):

    trees = []
    # for filename in file_in:
    #     if filename.endswith('.aln'):
    #         path_tree = os.path.join(directory, filename)
    #         tree_cluster = filename.split('.')[1]

    tree = Tree(path_tree)
    for node in tree.traverse():
        if node.is_leaf():
            name = node.name
            if name.find("|") != -1:
                split = name.split("|")
                name  = split[0]
                taxon = split[1]
                node.add_feature("name", name)
    # print(tree)

    trees.append(tree.write())
    print(trees)
    return trees


# def write_tree_file(trees, concat_tree):
#     with open(concat_tree, 'w') as outtree:
#         for line in trees:
#             outtree.write(str(line)+'\n')


if __name__ == '__main__':
    # concat_tree = '/home/issa/Documents/stage/raxml/clusters_Trimal/bestTree/concatTree.tre'
    # directory = '/home/issa/Documents/stage/raxml/clusters_Trimal/bestTree'
    # file_in   =  os.listdir(directory)
    path_tree = '/home/issa/Documents/stage/raxml/clusters_Gblock/bestTree/RAxML_bestTree.cluster_133.fasta.aln-gb.phy.tree'
    my_tree = read_newick(path_tree)
    # my_tree = read_newick(file_in, directory)
    # w = write_tree_file(my_tree, concat_tree)
