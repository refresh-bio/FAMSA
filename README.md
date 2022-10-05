# FAMSA

[![GitHub downloads](https://img.shields.io/github/downloads/refresh-bio/famsa/total.svg?style=flag&label=GitHub%20downloads)](https://github.com/refresh-bio/FAMSA/releases)
[![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/famsa.svg?style=flag&label=Bioconda%20downloads)](https://anaconda.org/bioconda/famsa)
[![Biocontainer downloads](https://img.shields.io/endpoint?url=https%3A%2F%2Fmmseqs.com%2Fbiocontainer.php%3Fcontainer%3Dfamsa)](https://biocontainers.pro/tools/famsa)
[![GitHub Actions CI](../../workflows/GitHub%20Actions%20CI/badge.svg)](../../actions/workflows/main.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![Join the chat at https://gitter.im/refresh-bio/FAMSA](https://badges.gitter.im/refresh-bio/FAMSA.svg)](https://gitter.im/refresh-bio/FAMSA?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

![x86-64](https://img.shields.io/static/v1?label=%E2%80%8B&message=x86-64&color=yellow&logo=PCGamingWiki&logoColor=white)
![ARM](https://img.shields.io/static/v1?label=%E2%80%8B&message=ARM&color=yellow&logo=Raspberry%20Pi&logoColor=white)
![Apple M1](https://img.shields.io/static/v1?label=%E2%80%8B&message=Apple%20M1&color=yellow&logo=Apple&logoColor=white)
![Windows](https://img.shields.io/badge/%E2%80%8B-Windows-00A98F?logo=windows)
![Linux](https://img.shields.io/static/v1?label=%E2%80%8B&message=Linux&color=00A98F&logo=linux&logoColor=white)
![macOS](https://img.shields.io/badge/%E2%80%8B-macOS-00A98F?logo=apple)
[![PyPI](https://img.shields.io/pypi/v/pyfamsa?label=PyFAMSA)](https://pypi.org/project/pyfamsa)

Progressive algorithm for large-scale multiple sequence alignments.

## New features in FAMSA 2
* Fast guide tree heuristic called Medoid Tree (`-medoidtree` switch) for ultra-scale alignments:
  * the entire Pfam-A v33.1 in its largest NCBI variant (over 18 thousand families, 60 GB of raw FASTA files) was analyzed in 8 hours,	
  * the family PF00005 of 3 million ABC transporters was aligned in 5 minutes and 24 GB of RAM.
* Remarkable time and memory optimizations - SLINK has been replaced with Prim’s minimum spanning tree algorithm when constructing default (single linkage) guide trees. NOTE: This may change quality results slightly compared to FAMSA 1 due to different ties resolution.
* Neighbour joining guide trees (`-gt nj` option). NOTE: Neighbour joining trees are calculated with a use of original *O*(*N*<sup>3</sup>) algorithm, thus their applicability on large sets is limited (unless they are used as subtrees with Medoid Tree heuristic).
* Option for compressing output aligment to gzip (`-gz` switch).
* Compatibility with ARM64 8 architecture (including Apple M1).
* Duplicate removal - redundant sequences are by default removed prior the alignment and restored afterwards (feature introduced in revision 2.1.0). This can change output alignments when a family contains duplicates. The old behaviour can be obtained by using `-keep-duplicates` switch.
* Profile-profile alignments (available by specifying two input FASTA files; introduced in revision 2.2.0).


## Quick start

```bash
git clone https://github.com/refresh-bio/FAMSA
cd FAMSA && make

# align sequences with default parameters (single linkage tree)
./famsa ./test/adeno_fiber/adeno_fiber sl.aln

# align sequences using UPGMA tree with 8 computing threads, store the result in a gzip archive
./famsa -gt upgma -t 8 -gz ./test/adeno_fiber/adeno_fiber upgma.aln.gz

# export a neighbour joining guide tree to the Newick format
./famsa -gt nj -gt_export ./test/adeno_fiber/adeno_fiber nj.dnd

# align sequences with the previously generated guide tree
./famsa -gt import nj.dnd ./test/adeno_fiber/adeno_fiber nj.aln

# align sequences with an approximated medoid guide tree and UPGMA subtrees
./famsa -medoidtree -gt upgma ./test/hemopexin/hemopexin upgma.medoid.aln

# export a distance matrix to the CSV format (lower triangular) 
./famsa -dist_export ./test/adeno_fiber/adeno_fiber dist.csv

# export a pairwise identity (PID) matrix to the CSV format (square) 
./famsa -dist_export -pid -square_matrix ./test/adeno_fiber/adeno_fiber pid.csv

# profile-profile alignment without refining output 
./famsa -refine_mode off ./test/adeno_fiber/upgma.no_refine.part1.fasta ./test/adeno_fiber/upgma.no_refine.part2.fasta pp.fasta
```


## Installation and configuration

FAMSA comes with a set of [precompiled binaries](https://github.com/refresh-bio/FAMSA/releases) for Windows, Linux, and macOS. They can be found under Releases tab. 
The software is also available on [Bioconda](https://anaconda.org/bioconda/famsa):
```
conda install -c bioconda famsa
```
For detailed instructions how to set up Bioconda, please refer to the [Bioconda manual](https://bioconda.github.io/user/install.html#install-conda). 
A user-friendly [PyFAMSA](https://github.com/althonos/pyfamsa) module authored by [Martin Larralde](https://github.com/althonos/) allows running analyzes directly from Python.
Finally, FAMSA can be built from the sources distributed as:

* Visual Studio 2019 solution for Windows,
* MAKE project for Linux and OS X (g++-5 required, g++-8 recommended).

FAMSA can be built for x86-64 and ARM64 8 architectures (including Apple M1 based on ARM64 8.4 core) and takes advantage of AVX2 (x86-64) and NEON (ARM) CPU extensions. The default target platform is x86-64 with AVX2 extensions. This, however, can be changed by setting `PLATFORM` variable for `make`:

```bash
make PLATFORM=none    # unspecified platform, no extensions
make PLATFORM=sse4    # x86-64 with SSE4.1 
make PLATFORM=avx     # x86-64 with AVX 
make PLATFORM=avx2    # x86-64 with AVX2 (default)
make PLATFORM=native  # x86-64 with AVX2 and native architecture
make PLATFORM=arm8    # ARM64 8 with NEON  
make PLATFORM=m1      # ARM64 8.4 (especially Apple M1) with NEON 
```   

Note, that x86-64 binaries determine the supported extensions at runtime, which makes them backwards-compatible. For instance, the AVX executable will also work on SSE-only platform, but with limited performance. An additional `make` option can be used to force static linking (may be helpful when binary portability is desired): `make STATIC_LINK=true`

The latest speed improvements in FAMSA limited the usefullness of the GPU mode. Thus, starting from the 1.5.0 version, there is no support of GPU in FAMSA. If maximum throughput is required, we encourage using new medoid trees feature (`-medoidtree` switch) which allows processing gigantic data sets in short time (e.g., the familiy of 3 million ABC transporters was analyzed in five minutes). 


## Usage

`famsa [options] <input_file> [<input_file_2>] <output_file>`

Positional parameters:
* `input_file`, `input_file_2` - input files in FASTA format (first input can be replaced with STDIN string to read from standard input); action depends on the number of input files:
    * one input - multiple sequence alignment (input gaps, if present, are removed prior the alignment),
	* two inputs - profile-profile aligment (gaps are preserved).
* `output_file` - output file (pass STDOUT when writing to standard output); available outputs:
    * alignment in FASTA format,
    * guide tree in Newick format (`-gt_export` option specified),
	* distance matrix in CSV format (`-dist_export` option specified).

Options:
* `-help` - show advanced options
* `-t <value>` - no. of threads, 0 means all available (default: 0)
* `-v` - verbose mode, show timing information (default: disabled)

* `-gt <sl | upgma | nj | import <file>>` - the guide tree method (default: sl):
    * `sl` - single linkage,
    * `upgma` - UPGMA,
    * `nj` - neighbour joining,
    * `import <file>` - import from a Newick file.
* `-medoidtree` - use MedoidTree heuristic for speeding up tree construction (default: disabled)
* `-medoid_threshold <n_seqs>` - if specified, medoid trees are used only for sets with `n_seqs` or more
* `-gt_export` - export a guide tree to output file in the Newick format
* `-dist_export` - export a distance matrix to output file in CSV format
* `-square_matrix` - generate a square distance matrix instead of a default triangle
* `-pid` - calculate percent identity (the number of matching residues divided by the shorter sequence length) instead of distance
* `-keep-duplicates` - keep duplicated sequences during alignment (default: disabled - duplicates are removed prior and restored after the alignment)
* `-gz` - enable gzipped output (default: disabled)
* `-gz-lev <value>` - gzip compression level [0-9] (default: 7)
* `-refine_mode <on | off | auto>` - refinement mode (default: `auto` - the refinement is enabled for sets <= 1000 seq.)


### Guide tree import and export

FAMSA has the ability to import/export alignment guide trees in Newick format. E.g., in order to generate a UPGMA tree from the *input.fasta* file and store it in the *tree.dnd* file, run:
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
* Single-linkage guide trees. While being very accurate, single-linkage trees can be established without storing entire distance matrix, which makes them suitable for large alignments. Although, the alternative guide tree algorithms like UPGMA and neighbour joining are also provided.
* The new heuristic based on K-Medoid clustering for generating fast guide trees. Medoid trees can be calculated in *O*(*N* log*N*) time and work with all types of subtrees (single linkage, UPGMA, NJ). The heuristic can be enabled with `-medoidtree` switch and allow aligning millions of sequences in minutes.

## Experimental results
The analysis was performed on our extHomFam 2 benchmark produced by combining Homstrad (March 2020) references with Pfam 33.1 families (NCBI variant). The data set was deposited at Zenodo: [https://zenodo.org/record/6524237](https://zenodo.org/record/6524237). The following algorithms were investigated:

| Name  | Version  | Command line  |
|---|---|---|
| Clustal&Omega;  | 1.2.4 |  `clustalo --threads=32 -i <input> -o <output>` |
| Clustal&Omega; iter2  | 1.2.4   | `clustalo --threads=32 --iter 2 -i <input> -o <output>` |
| MAFFT PartTree  |  7.453 | `mafft --thread 32 --anysymbol --quiet --parttree <input> -o <output>` |
| MAFFT DPPartTree  |  7.453 |  `mafft --thread 32 --anysymbol --quiet --dpparttree <input> -o <output>` |
| Kalign3 | 3.3.2 | `kalign -i <input> -o <output>` | 
| FAMSA  | 1.6.2  | `famsa -t 32 <input> <output>`  |
| FAMSA 2 | 2.0.1  | `famsa -t 32 -gz <input> <output>`  |
| FAMSA 2 Medoid | 2.0.1  | `famsa -t 32 -medoidtree -gt upgma -gz <input> <output>`  |


The tests were performed with 32 computing threads on a machine with AMD Ryzen Threadripper 3990X CPU and 256 GB of RAM. For each extHomFam 2 subset we measured a fraction of properly aligned columns (TC score) as well as a total running time and a maximum memory requirements. The results are presented in the figure below. Notches at boxplots indicate 95% confidence interval for median, triangle represent means. The missing series for some algorithm-set pairs indicate that the running times exceeded a week. Kalign3 failed to process 10 families (5 in second, 3 in fourth, and 2 in the largest subset). FAMSA 2 alignments were stored in gzip format (`-gz` switch). 

![extHomFam-v2-TC-comparison](https://user-images.githubusercontent.com/14868954/171652224-af88d980-5b49-4dcc-95e7-4de5dc152fb3.png)


The most important observations are as follows: 
* FAMSA 2 was superior in terms of accuracy to all the competitors. Only on the smallest families (*N* < 10k) Clustal&Omega; kept up with our algorithm.
* The advantage of FAMSA 2 increased with the number of sequences and reached 20-30 percent points for (100k, 250k] subset. 
* FAMSA 2 with medoid trees offered astonishing throughput (a familiy PF00005 of 3 million ABC transporters was aligned in 5 minutes) with accuracy only slightly inferior to that of the default single linkage trees.
* None of the competing algorithms was able to complete all the families in the largest [250k, 3M) subset.
* The memory requirements of FAMSA 2 allow ultra-scale analyzes at a desktop computer (24 GB for 3M sequences).

## Datasets

Benchmark data sets developed and used in the FAMSA study:
* extHomFam: [https://doi.org/10.7910/DVN/BO2SVW](https://doi.org/10.7910/DVN/BO2SVW)
* extHomFam 2: [https://zenodo.org/record/6524237](https://zenodo.org/record/6524237)

## Citing
[Deorowicz, S., Debudaj-Grabysz, A., Gudyś, A. (2016) FAMSA: Fast and accurate multiple sequence alignment of huge protein families. 
Scientific Reports, 6, 33964](https://www.nature.com/articles/srep33964)
