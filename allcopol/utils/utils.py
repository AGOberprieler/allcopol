import io, os, re, sys
from Bio import Phylo

# classes for storage of accession/taxon-specific information


class Accession:
    def __init__(self, name, alleles):
        self.name = name
        self.alleles = alleles
        self.partition = []


class Species:
    def __init__(self, species, ploidy, accession, alleles):
        self.name = species
        self.ploidy = ploidy
        self.accessions = [Accession(accession, alleles)]
        self.partition = []

    def get_alleles(self):
        """
        Summarizes allele IDs from all accessions of the species.
        """
        return sorted(set(sum((acc.alleles for acc in self.accessions), [])))


def clearStringIO(s):
    s.seek(0)
    s.truncate(0)


def rename_polyploids(input_data, species_polyploid, tree_repeats):
    """
    Appends locus/marker index (e.g. "_m2") to allele IDs from polyploids.
    
    Args:
        input_data: list of newick strings
        species_polyploid: list of Species objects
        tree_repeats: number of trees per locus
    
    Yields:
        newick strings with modified leaf labels
    """
    
    # forbidden characters within a node label: 
    # whitespace, parentheses, square brackets, comma, colon, semicolon
    non_label_chars = "([\\(\\)\\[\\],;:\\s])"  # enclosed within substitution group
    for i, line in enumerate(input_data):
        j = i//tree_repeats  # locus/marker index
        for sp in species_polyploid:
            for acc in sp.accessions:
                for allele in acc.alleles:
                    search_pattern = non_label_chars + re.escape(allele) + non_label_chars
                    replacement_string = r"\g<1>" + allele + "_m" + str(j) + r"\g<2>"
                    line = re.sub(search_pattern, replacement_string, line)
        yield line


def read_trees(fname, species_polyploid, tree_repeats):
    """
    Reads input gene trees.
    
    Args:
        fname: tree file (one newick string per line)
        species_polyploid: list of Species objects
        tree_repeats: number of trees per locus
                
    Returns:
        list of parsed trees
    """
    
    # remove trailing whitespace from each line
    with open(fname, "r") as f:
        input_data = [line.rstrip() for line in f.readlines()]

    # count and remove trailing newlines
    nt = 0
    for line in reversed(input_data):
        if line:
            break
        else:
            nt += 1
    if nt > 0:
        del input_data[-nt:]

    input_data = rename_polyploids(input_data, species_polyploid, tree_repeats) # append marker IDs
    input_data = io.StringIO("\n".join(input_data))  # create file handle
    input_data_trees = list(Phylo.parse(input_data, "newick"))

    return input_data_trees


def prune_tree(tree, alleles, suffix=""):
    """
    Prunes all leaves from a given tree which completely match one of a 
    given list of allele names. The matching pattern can be extended by 
    an optional suffix (regular expressions allowed).
    
    Args:
        tree: tree object to be pruned in place
        alleles: list of IDs
        suffix
    """
    for a in alleles:
        # search pattern matching preprocessed allele ID
        pattern = "^" + re.escape(a) + suffix + "$"
        # prune all matching leaves
        leaves_2b_pruned = [leaf.name for leaf in tree.get_terminals() if re.match(pattern, leaf.name)]
        for leaf in leaves_2b_pruned:
            tree.prune(name=leaf)


def build_mapping(prefix, allele_combinations):
    """
    Builds an allele mapping string for pseudo-diploids. For each 
    non-empty subset of alleles one hypothetical species is created.
    
    Args:
        prefix: controls name of created pseudo-diploids
        allele_combinations: list of strings, IDs already concatenated
                             for each pseudo-diploid
    
    Returns:
        allele mapping string
    """
    if len(allele_combinations) > 99:
        sys.exit("Error: maximal number of ancestors exceeded")
    
    allele_combinations = filter(None, allele_combinations)  # remove empty entries
    
    # create new species IDs incl. allele mapping
    species_list = [ "".join([ prefix, "___", "{0:02d}".format(j+1), ":", a]) 
                     for j, a in enumerate(allele_combinations) ]
    allele_map = ";".join(species_list)
    
    return allele_map


def build_PhyloNet_instruction(tree_IDs, mapping_orig, mapping_hyp, conf):
    """
    Creates an instruction for PhyloNet. Allele mappings for either true
    and hypothetical species are combined.
    
    Args:
        tree_IDs: specifies gene trees to be used
        mapping_orig: allele mapping string (diploids)
        mapping_hyp: allele mapping string (pseudo-diploids)
        conf: namespace of program options
    
    Returns:
        Instruction for species tree inference.
    """
    mapping_new = ";".join(filter(None, [mapping_orig, mapping_hyp]))
    if not mapping_new:
        sys.exit("Error: cannot build PhyloNet instruction, empty allele mapping")
        
    instruction = "".join([conf.phylonet_command, "(", tree_IDs, ") -a\n<", mapping_new, ">", ";\n"])
    return instruction


def print_results(title, tree, mapping, dist, conf):
    """
    Prints, writes, and optionally plots results.
    
    Args:
        title
        tree: inferred species tree (string)
        mapping: optimized allele mapping (string)
        dist: number of extra lineages
        conf: namespace of program options
    
    Warning:
        If the species tree is rerooted using an outgroup, it is not 
        checked wether the non-outgroup taxa are monophyletic!
    """
    
    # store newick string in place
    newick_line = io.StringIO()
    newick_line.write(tree)
    newick_line.seek(0)

    # modify tree
    final_tree = Phylo.read(newick_line, "newick")
    
    if conf.outgroup:
        final_tree.root_with_outgroup({'name':conf.outgroup[0]})
    if not conf.display_branch_length:
        clearStringIO(newick_line)
        # plain=True: branch lengths removed
        Phylo.write(final_tree, newick_line, "newick", plain=True)
        newick_line.seek(0)
        final_tree = Phylo.read(newick_line, "newick")

    outfile = open(os.path.join(conf.working_directory, "results.txt"), "a")
    result = "".join(["\n", title, ":\n", newick_line.getvalue(), "\nAllele mapping:\n", 
                      mapping, "\nTotal number of extra lineages: ", str(dist), "\n"])
    outfile.write(result)
    print(result)

    if conf.print_ascii_tree == True:
        Phylo.draw_ascii(final_tree)
        Phylo.draw_ascii(final_tree, outfile)
    
    outfile.close()


def read_allele_mapping(fname):
    """
    Reads name, species membership, allele names and ploidy level of 
    each accession.
    
    Args:
        fname: name of allele mapping file
    
    Returns:
        list of Species objects
    """
    species_list = []
    
    with open(fname, "r") as f:
        for i, line in enumerate(f, start=1):
            cols = line.split("\t")
            if len(cols) == 4:
                name = cols[0].strip()
                species = cols[1].strip()
                alleles = [ acc.strip() for acc in cols[2].split(",") if acc.strip() ]
                ploidy = int(cols[3].strip())
                
                if ploidy % 2:
                    sys.exit("Error: cannot handle uneven ploidy level")
                
                if name and species and alleles and ploidy:
                    try:
                        sp_ind = [sp.name for sp in species_list].index(species)  # species already known
                    except ValueError:
                        sp_ind = -1  # new species
                    
                    if sp_ind >= 0:
                        if ploidy != species_list[sp_ind].ploidy:
                            sys.exit("Error: conflicting ploidy levels within " + species)
                        species_list[sp_ind].accessions.append(Accession(name, alleles))  # update existing species
                    else:
                        species_list.append(Species(species, ploidy, name, alleles))  # create new species object
                    
                else:
                    sys.exit("Error: empty field in line " + str(i) + " of " + fname)
            
            elif len(cols) == 1:
                if cols[0].strip():
                    sys.exit("Error: single field in line " + str(i) + " of " + fname + " (4 expected)")
                else:
                    continue  # ignore empty lines
            else:
                sys.exit("Error: " + str(len(cols)) + " fields in line " + str(i) + " of " + fname + " (4 expected)")
    
    if species_list:
        return species_list
    else:
        sys.exit("Error: no valid entries in " + fname)

