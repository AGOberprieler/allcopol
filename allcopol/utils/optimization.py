import collections, itertools, random


def initial_solution(ploidy, n_alleles):
    """
    Creates random feasible solution.
    
    Args:
        ploidy: determines the number of pseudo-diploids (n/2)
        n_alleles: list, specifies how many alleles have to be 
                   partitioned (one entry per accession and marker)
    Returns:
        pseudo-diploid indices for each allele as list of lists
        
    Example:
        >>> initial_solution(4, [4,4])
        [[0,1,1,0], [1,1,0,0]]
    """
    start = []
    
    for n in n_alleles:
        parents = [i for i in range(ploidy//2) for j in range(2)]
        parents = random.sample(parents, n)
        start.append(parents)
    
    return start


def neighborhood_feasible(solution, ploidy):
    """
    Generates all feasible solutions in the neighborhood of an existing solution.
    
    Args:
        solution: current feasible solution
        ploidy
        
    Yields:
        new candidate solution along with 
        corresponding moves (locus, allele, FROM, TO)
    
    Example:
        >>> nb = neighborhood_feasible([[0,1,1,0], [1,1,0,0]], 4)   
        >>> nb.__next__()  # (new_solution, (move1, move2))
        ([[1, 0, 1, 0], [1, 1, 0, 0]], ((0, 0, 0, 1), (0, 1, 1, 0)))
    """
    # shift parent (if possible)
    for i, locus in enumerate(solution):
        c = collections.Counter(locus)
    
        # find accessible parents (i.e. less than 2 alleles)
        for j in range(ploidy//2):
            if c[j] < 2:
                # shift alleles from other parents
                for k, parent in enumerate(locus):
                    if parent != j:
                        locus_new = solution[i][:k] + [j] + solution[i][k+1:]
                        solution_new = solution[:i] + [locus_new] + solution[i+1:]
                        moves = ((i, k, parent, j),)  # locus, allele, FROM, TO
                        yield(solution_new, moves)
                
    # swap parents
    for i, locus in enumerate(solution):
        # group alleles by parent
        allele_groups = [[j for j in range(len(locus)) if locus[j] == parent] 
                            for parent in set(locus)]
        for comb in itertools.combinations(allele_groups, 2):
            # build pairs of alleles from different parents
            for pair in itertools.product(*comb):
                i0 = pair[0]
                i1 = pair[1]
                locus_new = locus[:]
                locus_new[i0], locus_new[i1] = locus_new[i1], locus_new[i0]
                solution_new = solution[:i] + [locus_new] + solution[i+1:]
                # (locus, allele, FROM, TO)
                moves = ((i, i0, solution[i][i0], solution[i][i1]), 
                         (i, i1, solution[i][i1], solution[i][i0]))
                yield(solution_new, moves)


def reverse_move(move):
    """
    Swaps FROM and TO from move.
    
    Args:
        move: tuple of indices (locus, allele, FROM, TO)
        
    Returns:
        reverse move
    """
    return move[0], move[1], move[3], move[2]


def is_tabu(moves, tabu_list):
    """
    Checks if moves contain any tabu-active allele-parent assignment.
    
    Args:
        moves: list [(locus, allele, FROM, TO), ...]
        
    Returns:
        tabu status
    """
    forbidden_assignments = [(locus, allele, TO) for (locus, allele, FROM, TO) in itertools.chain(*tabu_list)]
    tabus = [(m[0], m[1], m[3]) in forbidden_assignments for m in moves]
    
    if any(tabus):
        return True
    else:
        return False


def get_allele_partition(parent_indices, alleles_blockwise, ploidy):
    """
    Concatenates alleles names for each pseudo-diploid into strings.
    
    Args:
        parent_indices: pseudo-diploid indices for each allele 
                        as list of lists
        alleles_blockwise: suffixed allele names as list of lists,
                           one entry per accession and marker
        ploidy
    
    Returns:
        list of allele mapping strings for each pseudo-diploid
    
    Example:
        >>> get_allele_partition(
        ...     [[0,1,1], [0,0]],
        ...     [["a1_m0", "a2_m0", "a3_m0"], ["a1_m1", "a2_m1"]],
        ...     4
        ... )
        ['a1_m0,a1_m1,a2_m1', 'a2_m0,a3_m0']
    """
    partition = [[] for i in range(ploidy//2)]
    
    for i, block in enumerate(alleles_blockwise):
        for j, allele in enumerate(block):
            partition[parent_indices[i][j]].append(allele)
    
    # one string per parent
    partition = [",".join(parent) for parent in partition]
    return partition

