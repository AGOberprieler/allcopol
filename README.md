# AllCoPol

AllCoPol is a collection of tools for the analysis of polyploids, which allow to 
infer ancestral allele combinations as well as corresponding subgenome 
phylogenies.

## Installation

AllCoPol is hosted at the Python Package Index (PyPI), so it can be easily 
installed via pip:

```bash
python3 -m pip install allcopol
```

To run PhyloNet, Java 1.7 or later has to be installed.

## Contained tools

### allcopol

This is the main tool of the package implementing heuristic optimization of
ancestral allele combinations. To get a complete list of program options, type

```bash
allcopol --help
```


The program requires at least four arguments, 
specifying the input files (`-A`, `-G`), the number of supplied gene trees per marker 
(`-S`), and the path to a PhyloNet jar file (`-P`), 
which can be obtained from https://bioinfocs.rice.edu/phylonet (newest tested version: 3.8.0).
Besides, the tabu tenure (`-t`) and the number of iterations (`-i`) are crucial 
 for proper optimization. 

The allele mapping file (`-A`) is a tab-delimited text file with one line 
per accession and four columns: accession, taxon, allele IDs (comma separated), 
and ploidy level.
The gene tree file (`-G`) consists of newick strings supplied as one tree per line. 
For the trees, which are assumed to be rooted, topologies are sufficient while 
edge lengths, support values, etc. will be ignored. Multiple gene trees per 
marker can be supplied as consecutive lines in the tree file, e.g.

```text
<tree1 for marker1>
<tree2 for marker1>
<tree3 for marker1>
<tree1 for marker2>
<tree2 for marker2>
<tree3 for marker2>
...
```

Because the program cannot know from the input file which trees belong to the 
same marker, the `-S` option has to be set correctly (for the example above: `-S 3`).


_Minimal example:_

mapping.nw (input):

```text
acc1	sp1	A_1,A_2	2
acc2	sp1	B_1,B_2	2
acc3	sp2	C_1,C_2	2
acc4	sp2	D_1,D_2	2
acc5	sp3	E_1,E_2	2
acc6	sp4	F_1,F_2	2
acc7	sp5	G_1,G_2,H_1,H_2	4
```

trees.nw (input):

```text
(((E_2,(B_1,C_2)),(G_1,F_2)),((H_1,H_2),((A_2,A_1),((G_2,D_2),(F_1,((C_1,B_2),(D_1,E_1)))))));
((((F_2,(G_2,(G_1,F_1))),(A_1,A_2)),((((C_1,C_2),B_1),B_2),(((D_2,D_1),E_1),E_2))),(H_1,H_2));
((((G_2,((F_1,F_2),G_1)),A_1),((H_2,H_1),(D_1,E_1))),(((D_2,E_2),(C_2,(B_2,(C_1,B_1)))),A_2));
(((B_2,(C_2,B_1)),(H_2,H_1)),((A_2,A_1),(((G_2,G_1),(F_2,F_1)),((E_1,D_1),((D_2,E_2),C_1)))));
(((A_1,A_2),(H_2,H_1)),(((E_1,D_1),((F_1,F_2),(G_2,G_1))),(((E_2,D_2),(B_1,(C_2,B_2))),C_1)));
```

command:

```bash
allcopol -A mapping.txt -G trees.nw -S 1 -P PhyloNet_3.8.0.jar -t 5 -i 20
```

Setting the tabu tenure to zero and using reinitialization (-u), it is also 
possible to perform random restart hillclimbing instead of tabu search.
While this avoids extensive parameter tuning, it usually requires a higher 
number of iterations to obtain satisfactory solutions:

```bash
allcopol -A mapping.txt -G trees.nw -S 1 -P PhyloNet_3.8.0.jar -t 0 -u 1 -i 100
```

If runtime is limiting, the number of evaluated solutions per iteration can be 
limited via the -s option. Note that this may be at the expense of a lower final
solution quality.

_Parameter tuning:_

For single analyses, as far as runtime is not limiting, the simplest optimization strategy consists in repeated hillclimbing runs, allowing a high total number of iterations as described above.

In contrast, using a tabu list of adequate length allows to escape from local optima and therefore does not require to restart from random solutions.
This is usually more efficient because it allows to spend more time in promising regions of the search space.
Unfortunately, there is no simple rule how to specify the tabu tenure which is critical for proper function.
If it is set too small, the search can get stuck in local optima, if too high, many good solutions may be missed, which slows down the optimization progress.

For single reconstructions, the user could simply try out some different values and keep the most parsimonious result.

If many time-consuming (replicate) analyses have to be performed, it may be advisable to systematically measure the performance of different parameter settings before.
In this case, the obtained numbers of extra lineages for each setting should be averaged over multiple runs to account for the involved stochasticity.
Although suitable parameter settings are not easily transferable from one problem to another, it is often possible to find a good compromise which satisfies multiple similar analyses.
For this purpose, the number of extra lineages averaged over some exemplary inputs may be taken into account.

Whatever strategy is used, the applied heuristics cannot guarantee to find a global optimum.
In practise, it may be necessary to weigh up between solution quality and computational burden (mainly dependent on -i and, optionally, -s).
The preferred tabu tenure not only depends on the input, but also on the number of iterations, which should be multiple times higher, and the sample size parameter.

From our experience, the final topology is often found long before the lowest number of extra lineages, so it may not not necessary to put too much effort into finetuning.
If, on the other hand, repeated analyses end up with highly variable numbers of extra lineages, the used settings should questioned.



### create_indfile

This script takes a number of allele mapping strings (one per line, obtained 
by multiple runs of `allcopol` based on the same\* input files) as input and 
prints a matrix representation of the inferred allele partitions. The latter can
be used as input for CLUMPP or `align_clusters`.

\* The used gene trees may vary, but the underlying markers and their order in 
the tree files must be identical.

_Example:_

mappings.txt (input):

```text
sp2:C_1,C_2;sp3:D_1,D_2;sp1___01:B_1_m0,B_2_m0,B_1_m1,B_2_m1;sp1___02:A_1_m0,A_2_m0,A_1_m1,A_2_m1
sp2:C_1,C_2;sp3:D_1,D_2;sp1___01:A_1_m0,B_2_m0,B_1_m1,B_2_m1;sp1___02:A_2_m0,A_1_m1,A_2_m1,B_1_m0
sp2:C_1,C_2;sp3:D_1,D_2;sp1___01:A_1_m0,A_2_m0,A_1_m1,A_2_m1;sp1___02:B_1_m0,B_2_m0,B_1_m1,B_2_m1
sp2:C_1,C_2;sp3:D_1,D_2;sp1___01:A_1_m0,A_2_m0,A_1_m1,B_2_m1;sp1___02:B_1_m0,B_2_m0,B_1_m1,A_2_m1
```

command:

```bash
create_indfile mappings.txt sp1 > example.indfile
```

The second argument `sp1` specifies the name of the polyploid taxon, whose 
subgenomes have been inferred.

example.indfile (output):

```text
1 1 (x) 1 : 0 1
2 2 (x) 1 : 0 1
3 3 (x) 1 : 0 1
4 4 (x) 1 : 0 1
5 5 (x) 1 : 1 0
6 6 (x) 1 : 1 0
7 7 (x) 1 : 1 0
8 8 (x) 1 : 1 0

1 1 (x) 1 : 1 0
2 2 (x) 1 : 0 1
...
```

### align_clusters

This tool can be used to match clusters (subgenomes) among multiple
reconstructions. To avoid getting stuck in a local optimum, a tabu list is 
applied, whose size can be specified via the `-t` option. Compared to `allcopol`, 
the heuristic used for this step seems to be relatively robust with respect to 
its parameters.
`-n` sets the total number of optimization iterations.

Applied to the matrix representation written above, the command

```bash
align_clusters -n 50 -t 2 example.indfile
```

creates the output file example.permutations:

```text
2	1
2	1
1	2
1	2
```

and a second file containing the averaged cluster coefficients 
(example.clustering).

UPDATE: `align_clusters` can be replaced by the newer and much faster [Crimp](https://github.com/ulilautenschlager/crimp/) tool, e.g. by running `./crimp -n 20 example.indfile`.
In this case, the relevant output file will be named "example.indfile.permutations".
For small inputs as in this example, it is also possible to perform exact rather than heuristic optimization using the `-x` option, e.g. with `./crimp -x example.indfile`.
By default, Crimp optimizes a slightly different objective function than the original `align_clusters`. To stick with the entropy-based one, you can use the `-e` option, however, scores will differ by a factor of log(2) because `align_clusters` uses natural logarithms whereas Crimp uses binary ones.

### relabel_trees

Using the output of `align_clusters`, the species trees obtained by multiple
runs of `allcopol` can be relabeled to mitigate label switching.
`relabel_trees` expects three arguments, a file containing the inferred species
trees, the name of the analyzed polyploid taxon and a permutation file as 
written by `align_clusters`.

_Example:_

sp_trees.nw (input):

```text
(sp3,(sp1___02,(sp1___01,sp2)));
(sp3,(sp1___02,(sp1___01,sp2)));
(sp3,(sp1___01,(sp1___02,sp2)));
(sp3,(sp1___01,(sp1___02,sp2)));
```

The command

```bash
relabel_trees sp_trees.nw sp1 example.permutations
```

prints the relabeled trees:

```text
(sp3,(sp1_P1,(sp1_P2,sp2)));
(sp3,(sp1_P1,(sp1_P2,sp2)));
(sp3,(sp1_P1,(sp1_P2,sp2)));
(sp3,(sp1_P1,(sp1_P2,sp2)));
```

Now that the subgenomes are labeled according to their homology, 
conventional consensus methods can be applied to the trees.

