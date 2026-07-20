# FAMSA2

[![GitHub downloads](https://img.shields.io/github/downloads/refresh-bio/famsa/total.svg?style=flag&label=GitHub%20downloads)](https://github.com/refresh-bio/FAMSA/releases)
[![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/famsa.svg?style=flag&label=Bioconda%20downloads)](https://anaconda.org/bioconda/famsa)
[![Build and tests](../../workflows/Build%20and%20tests/badge.svg)](../../actions/workflows/main.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

![x86-64](https://img.shields.io/static/v1?label=%E2%80%8B&message=x86-64&color=yellow&logo=PCGamingWiki&logoColor=white)
![ARM](https://img.shields.io/static/v1?label=%E2%80%8B&message=ARM&color=yellow&logo=Raspberry%20Pi&logoColor=white)
![Apple M1](https://img.shields.io/static/v1?label=%E2%80%8B&message=Apple%20M1&color=yellow&logo=Apple&logoColor=white)
![Windows](https://img.shields.io/badge/%E2%80%8B-Windows-00A98F?logo=windows)
![Linux](https://img.shields.io/static/v1?label=%E2%80%8B&message=Linux&color=00A98F&logo=linux&logoColor=white)
![macOS](https://img.shields.io/badge/%E2%80%8B-macOS-00A98F?logo=apple)
[![PyPI](https://img.shields.io/pypi/v/pyfamsa?label=PyFAMSA)](https://pypi.org/project/pyfamsa)
[![WebInterface](https://img.shields.io/badge/%20-Web%20interface-blue.svg?logo=data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCA2MDAgNjAwIj48cGF0aCBmaWxsPSIjZmZmIiBkPSJNNDUwIDQwMEgxODBjLTYwIDAtMTEwLTQ1LTExMC0xMDAgMC01MiA0My05NSA5Ny0xMDAgMjEtNzQgODgtMTMwIDE2OC0xMzAgNzkgMCAxNDYgNTQgMTY4IDEyNSA2NiAzIDEyMiA1MCAxMjIgMTA1IDAgNTgtNTEgMTAwLTE3NSAxMDB6Ii8+PC9zdmc+)](https://refresh-bio.github.io/tools/famsa/)

FAMSA2 is a progressive algorithm for large-scale multiple sequence alignments:
* the entire Pfam-A v37.0 (~22 thousand families, ~62 million sequences) was analyzed in 8 hours,	
* the family PF00005 of 3 million ABC transporters was aligned in 5 minutes and 18 GB of RAM.

## Overview and features

![extHomFam-SP-comparison](./img/overview.png)

## Online version
The simplest way to align proteins is to run FAMSA directly in the browser without installation. 
Because it is compiled to WebAssembly, all computations run locally on the device &ndash; your data never leaves your machine.
The <a href="https://refresh-bio.github.io/tools/famsa/">online version</a> supports input files up to 100,000 sequences.

## Quick start

```bash
git clone https://github.com/refresh-bio/FAMSA --recursive
cd FAMSA && gmake

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

* Visual Studio 2022 solution for Windows,
* GNU Make project for Linux and macOS (gmake 4.3 and gcc/g++ 11 or newer required).

FAMSA can be built for x86-64 and ARM64 8 architectures (including Apple M1 based on ARM64 8.4 core) and takes advantage of AVX2 (x86-64) and NEON (ARM) CPU extensions. The default target platform is x86-64 with AVX2 extensions. This, however, can be changed by setting `PLATFORM` variable for `gmake`:

```bash
gmake PLATFORM=none    # unspecified platform, no extensions
gmake PLATFORM=sse4    # x86-64 with SSE4.1 
gmake PLATFORM=avx     # x86-64 with AVX 
gmake PLATFORM=avx2    # x86-64 with AVX2 (default)
gmake PLATFORM=avx512  # x86-64 with AVX512 
gmake PLATFORM=native  # x86-64 with native architecture
gmake PLATFORM=arm8    # ARM64 8 with NEON  
gmake PLATFORM=m1      # ARM64 8.4 (especially Apple M1) with NEON 
```   

Note, that x86-64 binaries determine the supported extensions at runtime, which makes them backwards-compatible. For instance, the AVX executable will also work on SSE-only platform, but with limited performance. An additional `gmake` option can be used to force static linking (may be helpful when binary portability is desired): `gmake STATIC_LINK=true`

### macOS support 

Due to usage of C++ features not supported by clang compiler, FAMSA requires g++ for compilation. As `g++` command on macOS systems maps to `clang`, one needs to install gcc/g++ and specify its binary using gmake arguments. For instance, if gcc/g++ 11 is available (the lowest supported version), please use the following command line:
```bash
gmake -j CXX=g++-11 CC=gcc-11
```


## Usage

`famsa [options] <input_file> [<input_file_2>] <output_file>`

Positional parameters:
* `input_file`, `input_file_2` - input files in FASTA format (optionally gzipped); first input can be replaced with STDIN string to read from standard input; action depends on the number of input files:
    * one input - multiple sequence alignment (input gaps, if present, are removed prior the alignment),
    * two inputs - profile-profile aligment (gaps are preserved).
* `output_file` - output file (pass STDOUT when writing to standard output); available outputs:
    * alignment in FASTA format,
    * guide tree in Newick format (`-gt_export` option specified),
	* distance matrix in CSV format (`-dist_export` option specified).

Options:
* `-help` - show advanced options
* `-t <value>` - no. of threads, NOTE: exceeding number of physical (not logical) cores decreases performance, 0 indicates half of all the logical cores (default: 0)
* `-v` - verbose mode, show timing information (default: disabled)

* `-gt <sl | upgma | nj | import <file>>` - the guide tree method (default: `sl`):
    * `sl` - single linkage,
    * `upgma` - UPGMA,
    * `nj` - neighbour joining,
    * `import <file>` - import from a Newick file.
* `-medoidtree` - use medoid tree for fast approximated guide trees (default: disabled)
* `-medoid_seeds <value>` - number of seeds *k* in medoid trees (default: 100)
* `-medoid_threshold <value>` - minimum subset size *M* to apply medoid tree heuristic (default: 2000)

* `-gt_export` - export a guide tree to output file in the Newick format
* `-dist_export` - export a distance matrix to output file in CSV format
* `-square_matrix` - generate a square distance matrix instead of a default triangle
* `-pid` - calculate percent identity (the number of matching residues divided by the shorter sequence length) instead of distance
* `-keep_duplicates` - keep duplicated sequences during alignment (default: disabled - duplicates are removed prior and restored after the alignment)
* `-gz` - enable gzipped output (default: disabled)
* `-gz_lev <value>` - gzip compression level [0-9] (default: 7)
* `-trim_columns <fraction>` - remove columns with less than `fraction` of non-gap characters
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
* The new heuristic based on K-Medoid clustering for generating fast guide trees to align millions of sequences in minutes. Medoid trees can be calculated in *O*(*MN* log*N*) time with *M* being the maximum subtree size. Medoid trees work with all types of subtrees (single linkage, UPGMA, NJ). The heuristic can be enabled with `-medoidtree` and is controlled by two parameters: the number of seed sequences *k* (`-medoid_seeds`) and the maximum subtree size *M* (`-medoid_threshold`).


## Experimental results
The analysis was performed on our extHomFam v37.0 benchmark produced by combining Homstrad references with Pfam v37.0 families (see Data sets section). The following algorithms were investigated:

| Name  | Version  | Command line  |
|---|---|---|
| Clustal&Omega;  | 1.2.4 |  `./clustalo --threads=32 -i <input> -o <output>` |
| MAFFT DPPartTree  |  7.526 |  `./mafft --thread 32 --anysymbol --quiet --dpparttree <input> -o <output>` |
| Kalign3 | 3.4.1 | `./kalign --type protein -n 32 -f fasta -i <input> -o <output>` | 
| Muscle5 | 5.3 | `./muscle -super5 <input> -output <output> --threads 32` |
| T-Coffee regressive | 13.46.0.919e8c6b |  `./clustalo --threads=32 -i <input> --guidetree-out <guide_tree> --force -o /dev/null`<br>`./t_coffee -reg -reg_method mafftsparsecore_msa -reg_tree <guide_tree> -seq <input> -reg_nseq 100 -reg_thread 32 -outfile <output>` |
| FAMSA  | 1.1  | `./famsa -t 32 <input> <output>`  |
| FAMSA2 | 2.4.1  | `./famsa -t 32 -gz <input> <output>`  |
| FAMSA2 Medoid | 2.4.1  | `./famsa -t 32 -medoidtree -gz <input> <output>`  |


The tests were performed with 32 computing threads on a machine with AMD Epyc 9554 CPU and 1152 GiB (approx. 1237 GB) of RAM. We measured a fraction of properly aligned residue pairs and columns (SP and TC scores, respectively) as well as a total running time and a peak memory usage. The results are presented in the figure below. Notches at boxplots indicate 95% confidence interval for median, triangle represent means. FAMSA2 alignments were stored in gzip format (`-gz` switch). 

![extHomFam-SP-comparison](./img/extHomFam.png)


## Datasets

Datasets developed and used in the FAMSA2 study:
* extHomFam v37.0: structure-based benchmark combining Pfam v37.0 families with Homstrad reference alignments
  * full dataset: [https://doi.org/10.5281/zenodo.6524236](https://doi.org/10.5281/zenodo.6524236)
  * the largest 121 families downsampled: [https://doi.org/10.5281/zenodo.18418209](https://doi.org/10.5281/zenodo.18418209)
* afdb_clusters v1.0: AlphaFold-derived structure-based benchmark ([https://zenodo.org/records/16082639](https://zenodo.org/records/16082639))
* simulated_msa v1.0: simulated multiple sequence alignments with known phylogenies ([https://zenodo.org/records/15971353](https://zenodo.org/records/15971353))
* active_sites v1.0: enzyme domain sequences with annotated active sites from Pfam v37.1 ([https://zenodo.org/records/16023627](https://zenodo.org/records/16023627))

Older datasets:
* extHomFam 2: [https://zenodo.org/record/6524237](https://zenodo.org/record/6524237)
* extHomFam: [https://doi.org/10.7910/DVN/BO2SVW](https://doi.org/10.7910/DVN/BO2SVW)

## Citing

[Gudyś, A., Zieleziński, A., Notredame, C., Deorowicz, S. (2026) Fast and accurate multiple-protein-sequence alignment at scale with FAMSA2. Nature Biotechnology, https://doi.org/10.1038/s41587-026-03095-3](https://doi.org/10.1038/s41587-026-03095-3)

[Deorowicz, S., Debudaj-Grabysz, A., Gudyś, A. (2016) FAMSA: Fast and accurate multiple sequence alignment of huge protein families. 
Scientific Reports, 6, 33964](https://www.nature.com/articles/srep33964)

## Funding
This work was supported by the National Science Centre, Poland, project [DEC-2022/45/B/ST6/03032].

