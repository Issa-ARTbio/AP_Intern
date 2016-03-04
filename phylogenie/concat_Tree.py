#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''change protein name to proteome name in the cluster trees
'''
import re, os
from ete3 import Tree
import random
def read_newick (file_in, directory):

    trees = []
    for filename in file_in:
        if filename.endswith('.tree'):
            path_tree = os.path.join(directory, filename)
            tree_cluster = filename.split('.')[1]

            tree = Tree(path_tree)
            t = Tree("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;", format=1)
            for node in tree.traverse():
                if node.is_leaf():
                    name = node.name
                    if name.find("|") != -1:
                        split = name.split("|")
                        name  = split[0]
                        taxon = split[1]
                        node.add_feature("name", name)

        trees.append(tree.write())
    return trees


def write_tree_file(trees, concat_tree):
    with open(concat_tree, 'w') as outtree:
        for line in trees:
            # line = line.strip()
            print(line)
            outtree.write(str(line)+'\n')


if __name__ == '__main__':
    concat_tree = '/home/issa/Documents/stage/raxml/data/Gblock/bestTree/concatTree.tre'
    directory = '/home/issa/Documents/stage/raxml/data/Gblock/bestTree/'
    file_in   =  os.listdir(directory)

    tree = read_newick(file_in, directory)
    w = write_tree_file(tree, concat_tree)
