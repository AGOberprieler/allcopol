#!/usr/bin/env python3
"""
Convert inferred allele mappings into a format that
can be used by CLUMPP or align_clusters.py.
"""

import argparse, operator, re, sys

def main():
    parser = argparse.ArgumentParser(description="Convert inferred allele mappings into a format that can be used by CLUMPP or align_clusters.py.")
    parser.add_argument("allele_mapping_file", help="Inferred allele mappings (PhyloNet format), one line per run.")
    parser.add_argument("species", help="Species which has been splitted into pseudodiploids.")
    args = parser.parse_args()

    species = args.species

    with open(args.allele_mapping_file, "r") as infile:
        txt = infile.read()

    mapping_strings = [m.strip() for m in txt.split("\n") if m.strip()]

    p = re.compile(species + "___[^;]*")

    # get max. number of clusters
    n_clusters = max([len(p.findall(mapping_strings[1])) for m in mapping_strings])


    for i, mstring in enumerate(mapping_strings):
        mappings_pseudodipl = sorted(p.findall(mstring)) # constant suffix length required
        membership = dict()

        for j, mp in enumerate(mappings_pseudodipl):
            alleles = mp.split(":")[1]
            for allele in alleles.split(","):
                membership[allele] = j

    #    for allele in sorted(membership.keys()):
    #        membership_coefficents = ["0"] * n_clusters
    #        membership_coefficents[membership[allele]] = "1"
    #        print(" ".join([str(i), str(i), "(x) 1 :"] + membership_coefficents))

        for k, allele in enumerate(sorted(membership.keys())):
            membership_coefficents = ["0"] * n_clusters
            membership_coefficents[membership[allele]] = "1"
            print(" ".join([str(k+1), str(k+1), "(x) 1 :"] + membership_coefficents))

        print("")

if __name__ == "__main__":
    main()
