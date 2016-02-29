#!/usr/bin/env python
# -*- coding: utf-8 -*-



from ete3 import Tree


tree_Gblock = '/home/issa/Documents/stage/raxml/specie_tree_Gblock/RAxML_bestTree.specie_TREE.tree'
t_trimal = '/home/issa/Documents/stage/raxml/specie_tree_Trimal/RAxML_bestTree.specie_TREE_trimal.tree'


g = open(tree_Gblock)
t = open(t_trimal)
gblock, trimal = [], []
for line in g:
    gblock.append(line)
for line in t:
    trimal.append(line)
gb = Tree(tree_Gblock)
tr = Tree(t_trimal)
print (tr)
# rf, max_rf, common_leaves, parts_t1, parts_t2 = gb.robinson_foulds(tr)
# print ("RF distance is %s over a total of %s" %(rf, max_rf))
