# FAMSA

[![GitHub downloads](https://img.shields.io/github/downloads/refresh-bio/famsa/total.svg?style=flag&label=GitHub%20downloads)](https://github.com/refresh-bio/FAMSA/releases)
[![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/famsa.svg?style=flag&label=Bioconda%20downloads)](https://anaconda.org/bioconda/famsa)
[![GitHub Actions CI](../../actions/workflows/main.yml/badge.svg)](../../actions/workflows/main.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Algorithm for large-scale multiple sequence alignments (400k proteins in 2 hours and 8GB of RAM)

## Quick start

```bash
git clone https://github.com/refresh-bio/FAMSA
cd FAMSA && make

# align sequences with default parameters (single linkage tree)
./famsa ./test/adeno_fiber/adeno_fiber sl.aln

# align sequences using UPGMA tree with 8 computing threads, store the result in the GZ archive
./famsa -gt upgma -t 8 -gz ./test/adeno_fiber/adeno_fiber upgma.aln.gz

# export a neighbour joining guide tree to the Newick format
./famsa -gt nj -gt_export ./test/adeno_fiber/adeno_fiber nj.dnd

# align sequences with the previously generated guide tree
./famsa -gt import nj.dnd ./test/adeno_fiber/adeno_fiber nj.aln

# export distance matrix to CSV format (lower triangular) 
./famsa -dist_export ./test/adeno_fiber/adeno_fiber dist.csv

# export pairwise identity (PID) matrix to CSV format (square) 
./famsa -dist_export -pid -square_matrix ./test/adeno_fiber/adeno_fiber pid.csv
```


## Installation and configuration

FAMSA comes with a set of [precompiled binaries](https://github.com/refresh-bio/FAMSA/releases) for Windows, Linux, and OS X. They can be found under Releases tab. 
The software is also available on [Bioconda](https://anaconda.org/bioconda/famsa):
```
conda install -c bioconda famsa
```
For detailed instructions how to set up Bioconda, please refer to the [Bioconda manual](https://bioconda.github.io/user/install.html#install-conda).
FAMSA can be also built from the sources distributed as:

* Visual Studio 2015 solution for Windows,
* MAKE project for Linux and OS X (g++-5 required, g++-8 recommended).

At the top of the makefile there are several switches controlling building process. These are:
* STATIC_LINK - enable static linking (default: false); may be helpful when binary portability is desired,
* NO_AVX - prevent from using AVX and AVX2 extensions (default: false),
* NO_AVX2 - prevent from using AVX2 extensions (default: false),

Note, that FAMSA by default takes advantage of AVX and AVX2 CPU extensions. Pre-built binary detetermines supported instructions at runtime, thus it is multiplatform. However, one may encounter a problem when building FAMSA version on a CPU without AVX and/or AVX2. For this purpose NO_AVX and NO_AVX2 switches are provided. 

The latest speed improvements in FAMSA limited the usefullness of the GPU mode. Thus, starting from the 1.5.0 version, there is no support of GPU in FAMSA. If maximum throughput is required, we encourage using new medoid trees feature (`-medoidtree` switch) which allows processing gigantic data sets in a reasonable time (e.g., the familiy of 3 million ABC transporters was analyzed in five minutes) . 


## Usage

`famsa [options] <input_file> <output_file>`

Positional parameters:
* `input_file` - input file in FASTA format (pass STDIN when reading from standard input)
* `output_file` - output file (pass STDOUT when writing to standard output); available outputs:
    * alignment in FASTA format,
    * guide tree in Newick format (`-gt_export` option specified),
	* distance matrix in CSV format (`-dist_export` option specified).

Options:
*  -help - show advanced options
* `-t <value>` - no. of threads, 0 means all available (default: 0)
* `-v` - verbose mode, show timing information (default: disabled)

* `-gt <sl | upgma | import <file>>` - the guide tree method (default: sl):
    * `sl` - single linkage,
    * `upgma` - UPGMA,
    * `import <file>` - import from a Newick file.
* `-medoidtree` - use MedoidTree heuristic for speeding up tree construction (default: disabled)
* `-medoid_threshold <n_seqs>` - if specified, medoid trees are used only for sets with `n_seqs` or more
* `-gt_export` - export a guide tree to output file in the Newick format
* `-dist_export` - export a distance matrix to output file in CSV format
* `-square_matrix` - generate a square distance matrix instead of a default triangle
* `-pid` - generate percent identity instead of distance
* `-gz` - enable gzipped output (default: disabled)
* `-gz-lev <value>` - gzip compression level [0-9] (default: 7)


### Guide tree import and export

FAMSA has the ability to import/export alignment guide trees in Newick format. From version 1.5.0, the interface for tree management has changed:
* `-gt_export` option does not require an additional parameter with file name. Instead, the name of the output file (`<output_file_name>`) is used. Note, that the output alignment is not produced.
* `-gt_import` option is deprecated. Instead, use a syntax for specifying a tree type `-gt import <file>` where `file` is a Newick file name.

E.g., in order to generate a UPGMA tree from the *input.fasta* file and store it in the *tree.dnd* file, run:
```
famsa -gt upgma -gt_export input.fasta tree.dnd
``` 
To align the sequences from *input.fasta* using the tree from *tree.dnd* and store the result in *out.fasta*, run:
```
famsa -gt import tree.dnd input.fasta out.fasta
```  

Below one can find example guide tree file for sequences A, B, and C:
```
(A:0.1,(B:0.2,C:0.3):0.4);
```
Note, that when importing the tree, the branch lengths are not taken into account, though they have to be specified in a file for successful parsing. When exporting the tree, all the branches are assigned with length 1, thus only the structure of the tree can be restored (we plan to output real lengths in the future release):
```
(A:1.0,(B:1.0,C:1.0):1.0);
```
## Algorithms
The major algorithmic features in FAMSA are:
* Pairwise distances based on the longest common subsequence (LCS). Thanks to the bit-level parallelism and utilization of SIMD extensions, LCS can be computed very fast. 
* Single-linkage guide trees. While being very accurate, single-linkage trees can be established without storing entire distance matrix, which makes them suitable for large alignments. Although, alternative guide tree algorithms like UPGMA and neigbour joining ale also provided.
* The new heuristic based on K-Medoid clustering for generating fast guide trees. Medoid trees can be calculated in O(NlogN) time and work with all types of subtrees (single linkage, UPGMA, NJ). The heuristic can be enabled with `-medoidtree` switch. Medoid trees allows ultra-scale alignments (e.g., the family of 3 million ABC transporters was processed in five minutes).


## Experimental results



## Citing
[Deorowicz, S., Debudaj-Grabysz, A., Gudyś, A. (2016) FAMSA: Fast and accurate multiple sequence alignment of huge protein families. 
Scientific Reports, 6, 33964](https://www.nature.com/articles/srep33964)
