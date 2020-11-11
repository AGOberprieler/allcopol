#!/usr/bin/env python3

import configargparse, numpy, random, time
import collections, copy, itertools, os, re, sys
from Bio import Phylo
from allcopol.utils import *


def main():

    cparser = configargparse.ArgParser()

    cparser.add_argument(
        "-C", "--config-file", 
        required = False, 
        is_config_file = True
    )
    
    # obligate parameters:
    
    cparser.add_argument(
        "-A", "--allele-mapping", 
        required = True
    )
    cparser.add_argument(
        "-G", "--gene-trees", 
        required = True, 
        help = "Gene tree file (1 newick line per tree)."
    )
    cparser.add_argument(
        "-S", "--tree-sample-size", 
        required = True, 
        type = int, 
        help = "Number of gene trees per locus."
    )
    cparser.add_argument(
        "-P", "--phylonet-path", 
        required = True, 
        help = "Path to PhyloNet jar file."
    )

    # optional parameters:
    
    cparser.add_argument(
        "-o", "--outgroup", 
        required = False, 
        nargs = '+', 
        default = [],
        help = "Root species trees at outgroup. Warning: If multiple taxa are specified, \
        it is not checked, whether non-outgroup taxa are monophyletic."
    )
    cparser.add_argument(
        "-p", "--polyploids", 
        required = False, 
        nargs = '+', 
        default = [],
        help = "If specified, only explicitly mentioned polyploid taxa will be analyzed \
        (if present in the mapping file). The option can be used to entirely skip the \
        optimization procedure (eg. \"-p none\")"
    )
    cparser.add_argument(
        "-m", "--max-procs", 
        required = False, 
        type = int,
        default = 1,
        help = "Maximum number of PhyloNet instances to be run in parallel. \
        Warning: (Over)parallelization causes significant computational overhead."
    )
    cparser.add_argument(
        "-t", "--tabu-tenure", 
        required = False, 
        type = int, 
        default = 2,
        help = "Length of the tabu list (iterations)."
    )
    cparser.add_argument(
        "-s", "--nb-sample-size", 
        required = False, 
        type = float, 
        default = "inf",
        help = "(Max.) number of evaluated candidate solutions per iteration. Default: \"inf\" \
        (complete neighborhood)"
    )
    cparser.add_argument(
        "-i", "--max-iter", 
        required = False, 
        type = int, 
        default = 500,
        help = "Maximum number of tabu search iterations."
    )
    cparser.add_argument(
        "-u", "--max-unimproved", 
        required = False, 
        type = float, 
        default = "inf",
        help = "Reinitialize search if no improvement has been found for MAX_UNIMPROVED iterations."
    )
    cparser.add_argument(
        "-w", "--working-directory", 
        required = False, 
        default = ""
    )
    cparser.add_argument(
        "-c", "--phylonet-command", 
        required = False, 
        default = "Infer_ST_MDC", 
        help = "Only tested with Infer_ST_MDC (default)!"
    )
    cparser.add_argument(
        "-f", "--remove-files", 
        required = False, 
        action = "store_true", 
        default = False, 
        help = "When finished, remove every file in WORKING_DIRECTORY except results.txt"
    )
    cparser.add_argument(
        "-j", "--java-options", 
        required = False, 
        default = "-Xmx1200m"
    )
    cparser.add_argument(
        "-b", "--display-branch-length", 
        required = False, 
        action = "store_true", 
        default = False
    )
    cparser.add_argument(
        "-a", "--print-ascii-tree", 
        required = False, 
        action = "store_true", 
        default = False
    )
    cparser.add_argument(
        "-r", "--report-interval", 
        required = False, 
        type = float, 
        default = "inf", 
        help = "Print intermediate results every report_interval iterations."
    )
    cparser.add_argument(
        "-e", "--evaluate-initial", 
        required = False, 
        action = "store_true", 
        default = False, 
        help = "Evaluate initial solution (iteration 0) and report results."
    )
    cparser.add_argument(
        "-d", "--seed", 
        required = False, 
        type = int, 
        default = -1, 
        help = "Initialize pseudorandom number generator, default: -1 (disabled)"
    )
    cparser.add_argument(
        "--diploids-only", 
        action = "store_true", 
        required = False, 
        default = False, 
        help = "Calculate species tree without reconstructed pseudo-diploids."
    )
    cparser.add_argument(
        "--combine-taxa", 
        action = "store_true", 
        required = False, 
        default = False, 
        help = "Combine pseudo-diploids inferred for different polyploid taxa into one species tree"
    )

    conf = cparser.parse_args()

    # correct list arguments
    conf.outgroup = [y.strip() for x in conf.outgroup for y in x.split(",") if y.strip()]
    conf.polyploids = [y.strip() for x in conf.polyploids for y in x.split(",") if y.strip()]

    # initialize pseudorandom number generators
    if conf.seed != -1:
        random.seed(conf.seed)
        numpy.random.seed(conf.seed)
    
    # named tuple for evaluated solutions
    evaluated = collections.namedtuple("evaluated", ["solution", "tree", "cost", "iteration"])

    # create working directory
    if not conf.working_directory:
        conf.working_directory = time.strftime("%Y_%h_%d")
    i = 1
    dirname = conf.working_directory
    while os.path.exists(dirname):
        dirname = conf.working_directory + "_" + str(i)
        i += 1
    os.mkdir(dirname)
    conf.working_directory = dirname

    # print/write parsed options
    d = vars(conf)
    
    with open(os.path.join(conf.working_directory, "results.txt"), "w") as f:
        for k in sorted(d.keys()):
            f.write(k + ": " + str(d[k]) + "\n")
            print(k + ": " + str(d[k]))
        f.write("----------------------\n")
        print("----------------------")

    time_start = time.time()

    IO_files = set()

    # read input:
    
    # get name, species membership, allele IDs and ploidy level of each accession
    species_list = read_allele_mapping(conf.allele_mapping)

    species_diploid = list(filter(lambda x: x.ploidy == 2, species_list))
    species_polyploid = list(filter(lambda x: x.ploidy > 2, species_list))

    # create allele mapping string for true diploids (accession-wise)
    mapping_orig = ";".join(acc.name + ":" + ",".join(acc.alleles) for sp in species_diploid for acc in sp.accessions)
    # species-wise
    mapping_orig_spec = ";".join(sp.name + ":" + ",".join(sp.get_alleles()) for sp in species_diploid)

    # read and preprocess gene trees
    input_data_trees = read_trees(conf.gene_trees, species_polyploid, conf.tree_sample_size)
    tree_IDs = ",".join("g" + "{0:07d}".format(i+1) for i in range(len(input_data_trees)))


    if conf.diploids_only:
        print ("\nDiploids only:")
        print("#" * 14, "\n")
        
        # prune gene trees (remove alleles of all polyploid species)
        ptrees = copy.deepcopy(input_data_trees)
        accessions_2b_pruned = [acc for sp in species_polyploid 
                                    for acc in sp.accessions ]
        
        for tree in ptrees:
            for acc in accessions_2b_pruned:
                prune_tree(tree, acc.alleles, "_m[0-9]+")

        # write pruned trees into file (unused)
        fname_ptrees = os.path.join(conf.working_directory,  "diploids_pruned.txt")
        Phylo.write(ptrees, fname_ptrees, "newick", plain=True)
        IO_files.add(fname_ptrees)
        
        PhyloNet_instruction = build_PhyloNet_instruction(tree_IDs, mapping_orig_spec, "", conf)
        # call PhyloNet
        PhyloNet_single_job(ptrees, PhyloNet_instruction, "diploids", conf, IO_files)
        # evaluate output
        stree, dist, _ = get_results_single("diploids", conf)
        print_results("Species tree diploids", stree, mapping_orig_spec, dist, conf)


    for species in species_polyploid:
        if len(conf.polyploids) > 0 and species.name not in conf.polyploids:
            continue
        
        print("\nSpecies:", species.name)
        print("#" * (9+len(species.name)), "\n")
        
        # prune gene trees (remove alleles of all but one polyploid species)
        ptrees = copy.deepcopy(input_data_trees)
        accessions_2b_pruned = [acc for sp in species_polyploid
                                    for acc in sp.accessions
                                    if sp != species]
        
        for tree in ptrees:
            for acc in accessions_2b_pruned:
                prune_tree(tree, acc.alleles, "_m[0-9]+")

        # write pruned trees into file (unused)
        fname_ptrees = os.path.join(conf.working_directory, species.name + "_pruned.txt")
        Phylo.write(ptrees, fname_ptrees, "newick", plain=True)
        IO_files.add(fname_ptrees)
        
        # allele names incl. marker suffix
        alleles_blockwise = []
        
        for accession in species.accessions:
            # loop across loci/markers
            for j in range(0, len(ptrees), conf.tree_sample_size):
                tree = ptrees[j]
                # search patterns matching preprocessed alleles of current accession and marker
                pattern_list = ["^" + re.escape(allele) + "_m[0-9]+$" for allele in accession.alleles]
                alleles_blockwise.append([leaf.name for p in pattern_list 
                                                for leaf in tree.get_terminals()
                                                if re.match(p, leaf.name)])
        
        # check number of alleles
        valid = True
        for i_marker, n_alleles in enumerate(alleles_blockwise):
            if len(n_alleles) > species.ploidy:
                print("Error: too many alleles at marker " + str(i_marker + 1) + ":", file=sys.stderr)
                print(", ".join([re.sub("_m[^_]+$", "", al) for al in alleles_blockwise[i_marker]]), file=sys.stderr)
                valid = False
        if not valid:
            sys.exit()
        
        n_alleles = [len(block) for block in alleles_blockwise]
        
        # initialization
        current_solution = evaluated(initial_solution(species.ploidy, n_alleles), "", float("inf"), 0)
        best_solution = current_solution
        best_solution_run = current_solution
        n_unimproved = 0
        tabu_list = []
        
        if conf.evaluate_initial:
            print("iteration: 0")
            partition = get_allele_partition(current_solution.solution, alleles_blockwise, species.ploidy)
            mapping_hyp = build_mapping(species.name, partition)
            PhyloNet_instruction = build_PhyloNet_instruction(tree_IDs, mapping_orig_spec, mapping_hyp, conf)
            # call PhyloNet
            PhyloNet_single_job(ptrees, PhyloNet_instruction, "initial", conf, IO_files)
            # evaluate output
            stree, dist, _ = get_results_single("initial", conf)
            current_solution = evaluated(current_solution.solution, stree, dist, 0)
            best_solution = current_solution
            
            print("current solution:", str(current_solution.cost), "/ best so far:", str(best_solution.cost))
            print_results(
                "Species tree " + species.name, 
                best_solution.tree, 
                mapping_orig_spec + ";" + mapping_hyp, 
                best_solution.cost, 
                conf
            )
        
        # optimize allele partitions
        for iteration in range(1, conf.max_iter+1):
            
            print("iteration:", str(iteration))
            
            neighborhood = list(neighborhood_feasible(current_solution.solution, species.ploidy))
            ssize = int(min(conf.nb_sample_size, len(neighborhood)))
            candidates, moves = zip(*random.sample(neighborhood, ssize))
            
            partitions = [get_allele_partition(parent_indices, alleles_blockwise, species.ploidy) 
                          for parent_indices in candidates]
            
            PhyloNet_instructions = []
            for p in partitions:
                # allele mapping string for hypothetical parents
                mapping_hyp = build_mapping(species.name, p)
                instr = build_PhyloNet_instruction(tree_IDs, mapping_orig_spec, mapping_hyp, conf)
                PhyloNet_instructions.append(instr)
            
            # call PhyloNet
            input_files = []
            PhyloNet_batch(ptrees, PhyloNet_instructions, species.name, "TS", input_files, conf, IO_files)
            
            # evaluate output
            strees, dist, _ = get_results_multi(input_files, len(candidates), conf)

            updated = False
            n_rejected = 0

            for j in numpy.argsort(dist):
                if not is_tabu(moves[j], tabu_list) or dist[j] < best_solution_run.cost:  # best_solution.cost:
                    current_solution = evaluated(candidates[j], strees[j], dist[j], iteration)
                    tabu_list = ([[reverse_move(m) for m in moves[j]]] + tabu_list) [:conf.tabu_tenure]
                    updated = True
                    
                    if current_solution.cost < best_solution_run.cost:
                        best_solution_run = current_solution
                        n_unimproved = 0
                    else:
                        n_unimproved += 1
                        
                    if current_solution.cost < best_solution.cost:
                        best_solution = current_solution
                        
                    break
                else:
                    n_rejected += 1
            
            if not updated:
                tabu_list = tabu_list[:-1]
            
            species.partition = get_allele_partition(best_solution.solution, alleles_blockwise, species.ploidy)
            if n_rejected > 0:
                print("    " + str(n_rejected) + " move(s) rejected")
            print("    current solution:", str(current_solution.cost), "/ best so far:", str(best_solution.cost))
            
            # reinitialization
            if n_unimproved >= conf.max_unimproved:
                print("\n---\nreinitialization\n---\n")
                current_solution = evaluated(initial_solution(species.ploidy, n_alleles), "", float("inf"), 0)
                best_solution_run = current_solution
                n_unimproved = 0
                tabu_list = []
            
            # report intermediate results
            if iteration % conf.report_interval == 0 and iteration < conf.max_iter:
                mapping_hyp = build_mapping(species.name, species.partition)
                print_results(
                    "Species tree " + species.name, 
                    best_solution.tree, 
                    mapping_orig_spec + ";" + mapping_hyp, 
                    best_solution.cost, 
                    conf
                )
            
        # report results for current polyploid
        print("\n---------------\n\nBest solution found in iteration", str(best_solution.iteration))
        
        mapping_hyp = build_mapping(species.name, species.partition)
        print_results(
            "Species tree " + species.name, 
            best_solution.tree, 
            mapping_orig_spec + ";" + mapping_hyp, 
            best_solution.cost, 
            conf
        )

    # species tree for multiple polyploids
    if conf.combine_taxa:
        print ("\nCombined polyploids:")
        print("#" * 20, "\n")
        
        # prune gene trees (remove alleles of all but one polyploid species)
        ptrees = copy.deepcopy(input_data_trees)
        accessions_2b_pruned = [acc for sp in species_polyploid 
                                    for acc in sp.accessions 
                                    if len(conf.polyploids) > 0 and sp.name not in conf.polyploids]
        
        for tree in ptrees:
            for acc in accessions_2b_pruned:
                prune_tree(tree, acc.alleles, "_m[0-9]+")

        # allele mapping string for hypothetical ancestors
        mapping_hyp = []
        for species in species_polyploid:
            if species.name in conf.polyploids or len(conf.polyploids) == 0:
                mapping_hyp.append(build_mapping(species.name, species.partition))
        mapping_hyp = ";".join(mapping_hyp)

        PhyloNet_instruction = build_PhyloNet_instruction(tree_IDs, mapping_orig_spec, mapping_hyp, conf)

        # call PhyloNet
        PhyloNet_single_job(ptrees, PhyloNet_instruction, "combined", conf, IO_files)

        # evaluate output
        stree, dist, _ = get_results_single("combined", conf)
        print_results(
            "Species tree combined", 
            stree,
            mapping_orig_spec + ";" + mapping_hyp, 
            dist, 
            conf
        )

    # final cleanup
    if conf.remove_files:
        for fname in IO_files:
            os.remove(fname)
                    
    time_stop = time.time()
    
    print("Total time elapsed:", "%.2f" % (time_stop-time_start), "s")

if __name__ == "__main__":
    main()
