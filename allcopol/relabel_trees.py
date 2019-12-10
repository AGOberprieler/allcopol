#!/usr/bin/env python3

import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("tree_file", help="Species trees with inferred pseudo-diploids, one newick string per line.")
    parser.add_argument("taxon", help="Name of the taxon which has been split into pseudo-diploids.")
    parser.add_argument("permutation_file", help="Alignment of clusters among multiple clusterings.")
    args = parser.parse_args()

    with open(args.tree_file, "r") as treefile:
        with open(args.permutation_file, "r") as permfile:
            for tree, perm in zip(treefile, permfile):
                if not (tree and perm):
                    break
                
                perm = perm.strip().split("\t")
                tree_relabeled = tree
                
                for i, cluster in enumerate(perm):
                    old = args.taxon  + "___" + "{:02}".format(int(cluster))
                    new = args.taxon  + "_P" + str(i+1)
                    tree_relabeled = tree_relabeled.replace(old, new)
                
                print(tree_relabeled, end="")

                
if __name__ == "__main__":
    main()
