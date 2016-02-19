#!/usr/bin/env python
# -*- coding: utf-8 -*-


'''Represente la distribution de différence de position dans la sortie Gblocks
'''

import os
import matplotlib
import matplotlib.pyplot as P
import numpy as np

def read_file(file_in):

    cluster_labels = []
    initial= []
    gblock = []
    with open(file_in, 'r') as f:
        lines = f.read().splitlines()
        for i in range(len(lines)):
            line = lines[i]
            if line.startswith(my_search):
                name = line.split()[5]
                gb = name.split('/')[-1]
                cluster_numb = gb.split('.')[0]
                cluster_labels.append(cluster_numb)
                number_pos_Gblock = line.split()[7]
                gblock.append(number_pos_Gblock)
                # num_initial = tmp[13]
                num_initial = line.split() [13]
                initial.append(num_initial)
                print(num_initial)
    # print(cluster_labels)
    return cluster_labels, initial, gblock

def plotting_gblock (cluster_labels, initial, gblock):


    n, bins, patches = P.hist(initial, normed=1, histtype='stepfilled')
    P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)

    P.show()


if __name__ == '__main__':

    directory = '/home/issa/Documents/stage/muscle/Gblock_results2/htm/'
    list_dir = os.listdir(directory)
    my_search = 'New number of positions in'
    for file_name in list_dir:
        if file_name.endswith('.htm'):
            pos = os.path.join(directory, file_name)
            file_name = file_name.split('.')[0]

        name, inition, gb = read_file (pos)
        plotting_gblock (name, inition, gb)
