#!/usr/bin/env python


#by T-B-F (https://github.com/T-B-F/myBio/blob/master/scripts/remove_paralogous.py)
"""
script to remove paralogous proteins from the list of orthologous.
For each potential paralogous proteins, compare the blast scores of each of them with the other proteins
Keep the most similar one
"""
import itertools
import argparse, os, sys
import numpy as np
import pickle
import random
from math import exp
from operator import mul


def read_blast(path, kept_proteins, mth="bitscore"):
    """ read the blast result file
    """
    result = {}
    species = set([])
    with open(path) as inf:
        for line in inf:
            tmp = line.split()
            sp0 = tmp[0].split("|")[0]
            sp1 = tmp[1].split("|")[0]
            species.add(sp0)
            species.add(sp1)
            name = (tmp[0],tmp[1])
            if tmp[0] in kept_proteins and tmp[1] in kept_proteins:
                if name not in result:
                    result[name] = -1e6
                if float(tmp[-1]) > result[name]:
                    if mth == "bitscore":
                        result[name] = float(tmp[-1])
                    else:
                        print("Error, no other methods implemented yet")
    return result, species


def read_orthomcl(path):
    """ read orthologous groups of proteins from orthomcl
    """
    ortholist = []
    species = set()
    with open(path) as inf:
        for line in inf:
            tmp = line.split()
            group_name = tmp[0]
            dgroup = {}
            for spprot in tmp[1:]:
                sp, prot = spprot.split("|")
                species.add(sp)
                dgroup.setdefault(sp, []).append(spprot)
            lgroup = []
            for sp in species:
                if sp in dgroup:
                    lgroup.append(dgroup[sp][:])
            ortholist.append(lgroup[:])
    return ortholist

def read_proteinortho(path):
    """ read proteinortho results
    """
    ortholist = []
    cnt = 0
    with open(path) as inf:
        # header
        line = inf.readline()
        species = line[1:].split()[3:]
        for line in inf:
            if line[0] != "#":
                tmp = line.split()
                group_name = "group_{}".format(cnt)
                cnt += 1
                dgroup = {}
                for i, proteins in enumerate(tmp[3:]):
                    if proteins != "*":
                        prots = proteins.split(",")
                        # store proteinortho group list
                        ortholist.append(prots[:])
    return ortholist


def get_score(proti, protj, blast_scores_ij):
    """ get the score between two proteins depending of a dictionary of blast
    scores
    """
    present = False
    score_ij = 0
    name = (proti, protj)
    if name in blast_scores_ij :
        score_ij = blast_scores_ij[name]
        present = True
    return score_ij, present

def score_path(protpath, blast_scores, indexes_combi, stop_if_missing=True):
    """ compute the score of the path and if all links are presents
    """
    all_present = True
    score = 0
    for i, j in indexes_combi:
        proti = protpath[i]
        protj = protpath[j]
        scoreij, all_presentij = get_score(proti, protj, blast_scores)
        scoreji, all_presentji = get_score(protj, proti, blast_scores)
        score += scoreij + scoreji
        all_present *= all_presentij * all_presentji
        if all_present == False and stop_if_missing:
            #print proti, protj
            return -score, all_present
    return -score, all_present



def init_temperature(orthogroup, blast_result, number_of_iteration = 500):
    """ initialise temperature by performing random walk at infinite temperature
    """
    size = len(orthogroup)
    indexes = np.arange(len(orthogroup))
    indexes_combi = list(itertools.combinations(indexes, 2))
    T = np.inf
    deltas = []
    # start random walk
    #cur_path = [random.choice(group) for group in orthogroup]
    #try:
    cur_path = [random.choice(group) for group in orthogroup]
    #except:
    #raise ValueError("Unable to find a clique with all species")
    cur_score, all_present = score_path(cur_path, blast_result, indexes_combi)
    cnt_try = 0
    while all_present == False and cnt_try < 200:
        #print "here"
        cur_path = [random.choice(group) for group in orthogroup]
        cur_score, all_present = score_path(cur_path, blast_result, indexes_combi)
        cnt_try += 1
    if cnt_try == 200 and all_present == False:
        raise RuntimeError("Unable to find a starting path")
    # start iterations
    for m in range(number_of_iteration):
        #print m
        idx = np.random.randint(size)
        while len(orthogroup[idx]) == 1 :
            # is the group only made of one index?
            idx = np.random.randint(size)
        new_node = random.choice(orthogroup[idx])
        while new_node == cur_path[idx]:
            # is the node the same than the previous choosen node?
            new_node = random.choice(orthogroup[idx])
        # create new path
        new_path = cur_path[:]
        new_path[idx] = new_node
        new_score, all_present = score_path(new_path, blast_result, indexes_combi)
        if all_present:
            # scores are negatif
            delta_score = new_score - cur_score
            val = min(1, exp( -(delta_score) / T ))
            uni = np.random.uniform()
            if uni <= val:
                cur_score = new_score
                cur_path = new_path[:]
                deltas.append(delta_score)
    if deltas == []:
        for m in range(number_of_iteration):
            cur_path = [random.choice(group) for group in orthogroup]
            cur_score, all_present = score_path(cur_path, blast_result, indexes_combi, stop_if_missing=False)
            deltas.append(cur_score)
    Tinit = 100000 * max(deltas)
    return Tinit


def sa_solution(orthogroup, blast_result, number_of_start = 10,
                number_of_iteration = 20000, metropolis_criterion = 1e-6,
                verbose=False):
    """ sometimes the number of possibilities is too big to try a brute
    force approach. In these cases, a simulated annealing approach is used.
    Start N time:
       randomly select a path, compute the score
       mark as best of all
       during M steps:
           select a node in the path and replace it
           compute the scores
           if the score is better or worse < threshold
               keep this path
           if the score is the best of all
               keep as best of all
    """
    size = len(orthogroup)
    indexes = np.arange(len(orthogroup))
    indexes_combi = list(itertools.combinations(indexes, 2))
    best_score = 1e9
    best_path = -1
    cnt = 0
    tot = 0
    #fixing T
    if verbose:
        print ("init T")
    try:
        Tinit = init_temperature(orthogroup, blast_result)
    except:
        raise RuntimeError("Unable to find a starting path")
    alpha = 0.99
    for n in range(number_of_start):
        cur_path = [random.choice(group) for group in orthogroup]
        cur_score, all_present = score_path(cur_path, blast_result, indexes_combi)
        cnt_try = 0
        while all_present == False and cnt_try < 200:
            cur_path = [random.choice(group) for group in orthogroup]
            cur_score, all_present = score_path(cur_path, blast_result, indexes_combi)
            cnt_try += 1
        if cnt_try == 200 and all_present == False:
            raise RuntimeError("Unable to find a starting path")
        if all_present and cur_score < best_score:
            best_score = cur_score
            best_path = cur_path[:]
        T = Tinit
        for m in range(number_of_iteration):
            if verbose:
                sys.stdout.write("{:2}/{:2}\t{:5}/{:5}\r".format(n, number_of_start, m, number_of_iteration))
                sys.stdout.flush()
            idx = np.random.randint(size)
            while len(orthogroup[idx]) == 1 :
                # is the group only made of one index?
                idx = np.random.randint(size)
            new_node = random.choice(orthogroup[idx])
            while new_node == cur_path[idx]:
                # is the node the same than the previous choosen node?
                new_node = random.choice(orthogroup[idx])
            # create new path
            new_path = cur_path[:]
            new_path[idx] = new_node
            new_score, all_present = score_path(new_path, blast_result, indexes_combi)
            if all_present:
                # scores are negatif
                delta_score = new_score - cur_score
                if new_score < cur_score:
                    val = 1
                else:
                    val = exp( -(delta_score) / T )
                #print val
                uni = np.random.uniform()
                if uni <= val:
                    cur_score = new_score
                    cur_path = new_path[:]
                    if delta_score > 0:
                        cnt += 1
                # only if the path is complete,check the scores
                if new_score < best_score:
                    best_score = new_score
                    best_path = new_path[:]
                tot += 1
            T = alpha * T
    if verbose:
        if tot > 0:
            print ("Metropolis rate : {}".format(cnt/float(tot)))
        else:
            print ("Metropolis rate : {}/{}".format(cnt,tot))
    return best_path

def bruteforce_solution(orthogroup, blast_result):
    """ run a brute force solution to find the best possible path
    """
    best_possibility = -1
    best_score = 1e9
    indexes = np.arange(len(orthogroup))
    indexes_combi = list(itertools.combinations(indexes, 2))
    for k, possibility in enumerate(itertools.product(*orthogroup)):
        score, all_present = score_path(possibility, blast_result, indexes_combi)
        # only keep best score
        if score < best_score and all_present:
            best_score = score
            best_possibility = possibility[:]
    return best_possibility

def find_best_possibility(cnt, orthogroup, blast_result, cutoff, verbose=False):
    """ remove paralogous sequences and ambigous groups by keeping only
    the best 1:1 set of proteins based on the sum of their reciprocal
    similarities
    Groups with too many possibilities (combinations) are not takend into
    account
    """
    best_possibility = -1
    sizes = [len(group) for group in orthogroup]
    try:
        tot = reduce(mul, sizes)
    except:
        print ("Error\t{}\t{}\n{}".format(sizes, orthogroup, best_possibility))
    # stop if too many combinations are possible
    if tot >= cutoff:
        try:
            best_possibility = sa_solution(orthogroup, blast_result, verbose=verbose)
        except:
            print ("Error Simulated Annealing, {}".format(orthogroup))
    else:
        # brute force
        best_possibility = bruteforce_solution(orthogroup, blast_result)
    return best_possibility


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", action="store", required=True, dest="orthores", help="orthomcl raw results")
    parser.add_argument("-b", action="store", required=True, dest="blastfile", help="blast file or pickle object of the blast dictionary")
    parser.add_argument("-c", action="store", dest="cutoff", default=1000000, type=int, help="above which a group is considered to have to many member and will be optimised by SA")
    parser.add_argument("-o", action="store", required=True, dest="outputfile", help="new output file")
    parser.add_argument("-m", action="store", dest="maxsp", help="number of species", type=int, default=0)
    parser.add_argument("-k", action="store_true", default=True, dest="keep_true_one2one", help="keep or discard true one 2 one without paralogous [default=True]")
    parser.add_argument("-v", action="store_true", default=False, dest="verbose", help="activate verbose mode")
    parser.add_argument("-f", action="store", required=True, choices=["orthomcl", "proteinortho"], dest="orthoformat", help="orthology format")
    params = parser.parse_args()

    max_nb_species = params.maxsp
    if params.verbose:
        print ("Maximal number of species : ", max_nb_species)

    if params.verbose:
        print ("Reading ortho")

    if params.orthoformat == "orthomcl":
        ortholist = read_orthomcl(params.orthores) #, species)
    else:
        ortholist = read_proteinortho(params.orthores) #, species)

    # only keep interesting groups
    kept_orthogroup = list()
    kept_proteins = set()
    for orthogroup in ortholist:
        nb_group = len(orthogroup)
        if nb_group == max_nb_species:
            kept_orthogroup.append(orthogroup)
            for group in orthogroup:
                for prot in group:
                    kept_proteins.add(prot)
    # cleanup
    del ortholist
    if params.verbose:
        print("Number of group to check {}".format(len(kept_orthogroup)))
        print("Number of proteins {}".format(len(kept_proteins)))

    if params.verbose:
        print ("Reading Blast")

    pathblast, blastname = os.path.split(params.blastfile)
    if os.path.basename(params.blastfile) == "blast_pickle.dat":
        if params.verbose:
            print ("from pickle")
        with open(params.blastfile, "rb") as inf:
            blast_scores, species = pickle.load(inf)
    else:
        pathout = os.path.join(pathblast, "blast_pickle.dat")
        blast_scores, species = read_blast(params.blastfile, kept_proteins)
        with open(pathout, "wb") as outf:
            pickle.dump((blast_scores, species), outf)

    # cleanup
    del kept_proteins
    if params.verbose:
        print ("Checking groups")

    with open(params.outputfile, "w") as outf:
        #outf.write("\t".join(species)+"\n")
        cnt = 0
        # for each orthologuos groups,
        for orthogroup in kept_orthogroup:
            cnt += 1
            nb_group = len(orthogroup)
            if cnt %10 == 0:
                if params.verbose:
                    print("orthogroup {}/{}\n".format(cnt, len(kept_orthogroup)))
                    sys.stdout.flush()
            if nb_group == max_nb_species:
                #print orthogroup
                #continue
                sizes = [len(group) for group in orthogroup]
                # if only one protein per species, the group is kept without modifications
                # else ( sum != max_nb_species ) the paralogs are filtered out of the group
                if sum(sizes) != max_nb_species:
                    best_possibility = find_best_possibility(cnt, orthogroup, blast_scores, params.cutoff, verbose=params.verbose)
                    # only if one path among orthologous proteins is found
                    if best_possibility != -1:
                        new_orthogroup = best_possibility[:]
                        line = "\t".join(new_orthogroup)
                        outf.write(line+"\n")
                    else :
                        # the group is discarded if no pass found
                        print ("Problem {} {}".format(str(cnt)," ".join([",".join(group) for group in orthogroup])), file=sys.stderr)
                elif params.keep_true_one2one:
                    # group with only 1 protein per species are kept without modification
                    line = "\t".join([",".join(group) for group in orthogroup])
                    outf.write(line+"\n")

    sys.exit(0)


if __name__ == "__main__":
    main()
