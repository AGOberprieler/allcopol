# tested for PhyloNet 3.6.8
from Bio import Phylo
import concurrent, io, itertools, math, os, random, re, shutil, subprocess, sys, time

# call PhyloNet:


def PhyloNet_single_job(trees, PhyloNet_line, step, conf, IO_files):
    """
    Calls PhyloNet for single species tree reconstruction.
    
    Args:
        trees: gene trees in newick format (list of strings)
        PhyloNet_line: instruction for tree inference
        step: name of PhyloNet call
        conf: namespace of program options
        IO_files: set, stores names of all existing PhyloNet input/output
                  files for later cleanup
    
    Side effects:
        PhyloNet output is written, IO_files are updated.
    """
    fname_input = os.path.join(conf.working_directory, step + ".nex")
    fname_output = os.path.join(conf.working_directory, step + ".txt")
    
    with open(fname_input, "w") as f:
        f.write("#NEXUS\n\nBEGIN TREES;\n\n")
        for i, tree in enumerate(trees, 1):
            newick_line = io.StringIO()
            Phylo.write(tree, newick_line, "newick", plain=True)
            f.write("".join(["Tree g", "{0:07d}".format(i), " =\n", newick_line.getvalue()]))
        
        f.write("\nEND;\n\n\nBEGIN PhyloNet;\n\n")
        f.write(PhyloNet_line)
        f.write("\nEND;\n")

    # call PhyloNet
    os.system("".join(["java ", conf.java_options, " -jar ", conf.phylonet_path, 
                       " ", fname_input, " > ", fname_output]))
    
    # update file list
    IO_files.update([fname_input, fname_output])


def PhyloNet_batch(trees, PhyloNet_lines, name, step, input_files, conf, IO_files):
    """
    Calls PhyloNet for multiple species tree reconstructions.
    
    Args:
        trees: gene trees in newick format (list of strings)
        PhyloNet_lines: tree inference instructions
        name: name polyploid taxon
        step: name of PhyloNet calls
        input_files: list, stores names of currently processed 
                     PhyloNet input files
        conf: namespace of program options
        IO_files: set, stores names of all existing PhyloNet input/output
                  files for later cleanup
    
    Side effects:
        PhyloNet output is written, input_files and IO_files are updated.
    """
    input_files.clear()

    fname_head = os.path.join(conf.working_directory, name + "_" + step + "_head.nex")
    
    # write new head if necessary
    if fname_head not in IO_files:
        with open(fname_head, "w") as f:
            f.write("#NEXUS\n\nBEGIN TREES;\n\n")
            # write trees
            for i, tree in enumerate(trees, 1):
                newick_line = io.StringIO()
                Phylo.write(tree, newick_line, "newick", plain=True)
                f.writelines(["Tree g", "{0:07d}".format(i), " =\n", newick_line.getvalue()])
            f.write("\nEND;\n\n\n\nBEGIN PhyloNet;\n\n")
    
    chunk_size = max(len(PhyloNet_lines) // conf.max_procs, 1)
    n_chunks = math.ceil(len(PhyloNet_lines) / chunk_size)

    # split PhyloNet input into chunks
    for i in range(n_chunks):
        fname_chunk = os.path.join(conf.working_directory, name + "_" + step + "_stack{0:04d}.nex".format(i))
        shutil.copyfile(fname_head, fname_chunk)  # copy trees

        with open(fname_chunk, "a") as f:
            f.writelines(PhyloNet_lines[i*chunk_size : (i+1)*chunk_size])  # (out of range indices will be omitted)
            f.write("\n\nEND;\n")
        input_files.append(fname_chunk)

    processes = set()

    for infile in input_files:
        outfile = os.path.splitext(infile)[0] + ".txt"

        with open(outfile, "w") as fout:
            processes.add(subprocess.Popen(
                ["java", conf.java_options, "-jar", conf.phylonet_path, infile], stdout=fout
            ))

        if len(processes) >= conf.max_procs:
            os.wait()
            processes.difference_update(
                [p for p in processes if p.poll() is not None])

    # Check if all the child processes were closed
    for p in processes:
        if p.poll() is None:
            p.wait()
    
    # update file list
    IO_files.add(fname_head)
    IO_files.update(input_files)
    IO_files.update([os.path.splitext(fname)[0] + ".txt" for fname in input_files])


# get results:

def parse_output_chunk(data_chunk):
    """
    Parses PhyloNet output for one species tree reconstruction.
    
    Args:
        data chunk: list of four consecutive lines 
                    from PhyloNet output file
    
    Returns:
        species tree, input tree IDs, allele mapping, 
        number of extra lineages
    """
    match_tree = re.search("^(\\(.+;)\n?", data_chunk[2])
    match_ids = re.search("\\(([ ,g0-9]*)\\)", data_chunk[1])
    match_map = re.search("<(.+)>", data_chunk[1])
    match_dist = re.search("Total number of extra lineages:(.+)\n?", data_chunk[3])

    # check for valid matches and additional empty line
    if match_tree and match_ids and match_map and match_dist and data_chunk[0]=="\n":
        chunk_tree = match_tree.group(1)
        chunk_ids = match_ids.group(1)
        chunk_map = match_map.group(1)
        chunk_dist = int(float(match_dist.group(1)))
    else:
        sys.exit("error: cannot parse PhyloNet output!")
    
    return chunk_tree, chunk_ids, chunk_map, chunk_dist


def get_results_single(step, conf):
    """
    Evaluates PhyloNet output (single species tree reconstruction).
    
    Args:
        step: name of PhyloNet call
        conf: namespace of program options
    
    Returns:
        species tree, number of extra lineages, allele mapping
    """
    fname = os.path.join(conf.working_directory, step + ".txt")

    with open(fname, "r") as f:
        data_chunk = f.readlines()

    if not data_chunk:
        sys.exit("error: cannot parse PhyloNet output!")
    
    tree, ids, mapping, dist = parse_output_chunk(data_chunk)
    
    return tree, dist, mapping


def get_results_multi(input_files, n_PhyloNet_instructions, conf):
    """
    Evaluates PhyloNet output (multiple species tree reconstructions).
    
    Args:
        input_files: names of input files, which were passed to PhyloNet
        n_PhyloNet_instructions: number of species tree reconstructions
        conf: namespace of program options
    
    Returns:
        tuple of three equal-length lists:
        species trees, numbers of extra lineages, allele mappings 
    """
    chunk_size = max(1, n_PhyloNet_instructions // conf.max_procs)
    
    distances = [0] * n_PhyloNet_instructions
    trees = [0] * n_PhyloNet_instructions
    mappings = [0] * n_PhyloNet_instructions
    
    i = 0
    for fname in input_files:
        with open(os.path.splitext(fname)[0] + ".txt", "r") as f:
            for j in range(chunk_size):
                # read 4 lines from file
                data_chunk = [line for line in itertools.islice(f, 4)]
                
                if not data_chunk:
                    if (fname != input_files[-1]) or (j != n_PhyloNet_instructions % chunk_size):
                        sys.exit("error: cannot parse PhyloNet output!")
                    break
                
                trees[i], _, mappings[i], distances[i] = parse_output_chunk(data_chunk)
                i += 1
    return trees, distances, mappings

