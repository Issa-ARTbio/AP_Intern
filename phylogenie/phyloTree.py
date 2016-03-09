import ete3
from ete3 import PhyloTree, TreeStyle, Tree

alg =  '/home/issa/Documents/stage/muscle/final_alignement/cluster_9.fasta.aln'

def get_example_tree():

    # Performs a tree reconciliation analysis
    gene_tree_nw = '/home/issa/Documents/stage/raxml/clusters_Trimal/bestTree/RAxML_bestTree.cluster_9.fasta.aln'
    species_tree_nw = '/home/issa/Documents/stage/raxml/specie_tree_Trimal/RAxML_bestTree.specie_TREE_trimal.tree'
    genetree = PhyloTree(gene_tree_nw)
    sptree = PhyloTree(species_tree_nw)
    recon_tree, events = genetree.reconcile(sptree)
    recon_tree.link_to_alignment(alg)
    return recon_tree, TreeStyle()

if __name__ == "__main__":
    # Visualize the reconciled tree
    t, ts = get_example_tree()
    t.show(tree_style=ts)
    #recon_tree.render("phylotree.png", w=750)
