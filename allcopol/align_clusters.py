#!/usr/bin/env python3
# version: 05.12.19
# Author: Ulrich Lautenschlager
# changes by Tankred Ott:
#    put functionality into run function to be able
#    to import the script into other applications

import argparse, collections, itertools, numpy, os, random, re, scipy.stats, sys, time

solution = collections.namedtuple("solution", ["comb", "cost"])


def mean_ent(cluster_assignments):
    """
    Evaluate cost function for one solution.
    """
    nrow, ncol, nrun = cluster_assignments.shape
    csum = cluster_assignments.sum(axis=2)
    m = numpy.mean([scipy.stats.entropy(row) for row in csum])
    return m


def neighborhood_sample(solution, ssize):
    """
    Returns list of moves + list of (subsampled) neighborhood solutions
    """
    n_runs = len(solution)
    n_clusters = len(solution[0])
    
    swaps_per_run = n_clusters * (n_clusters-1) // 2
    max_size = n_runs * swaps_per_run
    
    ssize = min(ssize, max_size)
    
    combinations = list(itertools.combinations(range(n_clusters), 2))
    selection = random.sample(range(max_size), ssize)
    
    nb = [0] * ssize
    
    for i_nb, i_sel in enumerate(selection):
        i_run = i_sel // swaps_per_run
        i_comb = i_sel % swaps_per_run
        
        i1, i2 = combinations[i_comb]
        swapped = solution[i_run][:]
        swapped[i1], swapped[i2] = swapped[i2], swapped[i1]
        nb[i_nb] = solution[:i_run] + [swapped] + solution[i_run+1:]
    
    return selection, nb


def run(indfile, n_iter=400, tabu_tenure=10, sample_size=None, quiet=True):
    time_start = time.time()

    clustering = []
    runs = []  # list of clusterings

    with open(indfile, "r") as fin:
        for line in fin:
            # ignore additional info
            line = re.sub(r".*:", "", line)
            # get cluster assignments
            all = re.findall(r"[0-9]*[.]?[0-9]+", line)
            if all:
                clustering.append(all)
            else:
                if clustering:
                    runs.append(clustering)
                    clustering = []
    if clustering:
        runs.append(clustering)
        
    cluster_assignments = numpy.dstack([numpy.array(r, dtype="f8") for r in runs])
    n_rows, n_clusters, n_runs = cluster_assignments.shape
    if not quiet:
        print(cluster_assignments.shape)

    max_sample_size =  n_runs * n_clusters * (n_clusters-1) // 2
    if not sample_size:
        sample_size = max_sample_size
    if not quiet:
        print("Evaluate", sample_size, "solutions per iteration.")

    tabu_list = []
    
    # initial solution
    current_solution = [random.sample(range(n_clusters), n_clusters) for i in range(n_runs)]
    best_solution = current_solution
    cost = [0] * sample_size
    cost_best = float("inf")
    
    for iteration in range(n_iter):
        
        moves, neighborhood = neighborhood_sample(current_solution, sample_size)
        
        for j, candidate in enumerate(neighborhood):
            cluster_permuted = numpy.dstack([cluster_assignments[:,candidate[k],k] for k in range(n_runs)])
            cost[j] = mean_ent(cluster_permuted)
        
        updated = False
        
        for j in numpy.argsort(cost):
            if cost[j] < cost_best:
                i_best = iteration
                current_solution = neighborhood[j]
                best_solution = current_solution
                cost_best = cost[j]
                tabu_list = [moves[j]] + tabu_list
                tabu_list = tabu_list[:tabu_tenure]
                updated = True
                break
            else:
                if moves[j] not in tabu_list:
                    current_solution = neighborhood[j]
                    tabu_list = [moves[j]] + tabu_list
                    tabu_list = tabu_list[:tabu_tenure]
                    updated = True
                    break
        
        if not updated:
            tabu_list = tabu_list[:-1]
        
        if not quiet:
            print("best solution of iteration ", str(iteration), ": ", str(cost[j]))
            print("best solution so far: ", str(cost_best))
            #print("variance: ", str(numpy.var(cost)))
            print("------------------------")
            sys.stdout.flush()
        

    time_stop = time.time()
    if not quiet:
        print("Total time elapsed:", "%.2f" % (time_stop-time_start), "s")
        print("Best solution found in iteration", str(i_best))
    
    final_perm = numpy.array(best_solution, dtype="uint8") + 1

    cluster_permuted = numpy.dstack([cluster_assignments[:, best_solution[k], k] for k in range(n_runs)])
    q_matrix_avg = numpy.divide(cluster_permuted.sum(axis=2), n_runs)

    return final_perm, q_matrix_avg


def save(indfile, final_perm, q_matrix_avg):
    spl = os.path.splitext(indfile)

    fout = spl[0] + ".permutations"
    numpy.savetxt(fout, final_perm, fmt="%u", delimiter="\t")
    
    fout = spl[0] + ".clustering"
    #numpy.savetxt(fout, cluster_permuted.sum(axis=2), fmt="%u", delimiter="\t")
    numpy.savetxt(fout, q_matrix_avg, fmt="%1.6f", delimiter="\t")


def main():
    parser = argparse.ArgumentParser(description="Align clusters among multiple clusterings.")
    parser.add_argument("indfile", help="Membership coefficents (CLUMPP-compatible format).")
    parser.add_argument("-n", "--n-iter", default=400, type=int, help="Number of iterations.")
    parser.add_argument("-s", "--sample-size", type=int, help="Number of evaluated neighborhood solutions per iteration (default: all).")
    parser.add_argument("-t", "--tabu-tenure", default=10, type=int, help="Number of forbidden reverse moves.")
    args = parser.parse_args()

    final_perm, q_matrix_avg = run(
        args.indfile,
        n_iter=args.n_iter,
        tabu_tenure=args.tabu_tenure,
        sample_size=args.sample_size,
        quiet=False
    )
    
    print("List of permutations:")
    sys.stdout.flush()
    numpy.savetxt(sys.stdout.buffer, final_perm, fmt="%u", delimiter="\t")
    
    print("------------------------")
    
    print("Q-matrix:")
    sys.stdout.flush()
    #numpy.savetxt(sys.stdout.buffer, cluster_permuted.sum(axis=2), fmt="%u", delimiter="\t")
    numpy.savetxt(sys.stdout.buffer, q_matrix_avg, fmt="%1.6f", delimiter="\t")

    # write output files

    spl = os.path.splitext(args.indfile)

    fout = spl[0] + ".permutations"
    numpy.savetxt(fout, final_perm, fmt="%u", delimiter="\t")
    
    fout = spl[0] + ".clustering"
    #numpy.savetxt(fout, cluster_permuted.sum(axis=2), fmt="%u", delimiter="\t")
    numpy.savetxt(fout, q_matrix_avg, fmt="%1.6f", delimiter="\t")
    

if __name__ == "__main__":
    main()
