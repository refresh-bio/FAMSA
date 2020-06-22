# FAMSA

[![GitHub downloads](https://img.shields.io/github/downloads/refresh-bio/famsa/total.svg?style=flag&label=GitHub%20downloads)](https://github.com/refresh-bio/FAMSA/releases)
[![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/famsa.svg?style=flag&label=Bioconda%20downloads)](https://anaconda.org/bioconda/famsa)
[![C/C++ CI](https://github.com/refresh-bio/FAMSA-dev/workflows/C/C++%20CI/badge.svg)](https://github.com/refresh-bio/FAMSA-dev/actions)

Algorithm for large-scale multiple sequence alignments (400k proteins in 2 hours and 8BG of RAM)


## Installation and configuration

FAMSA comes with a set of precompiled binaries for Windows, Linux, and OS X. They can be found under Releases tab. 
Starting from 1.5.0 version there is no support of GPU in FAMSA. Use one of the previous releases if you need that feature.

The software can be also built from the sources distributed as:

* Visual Studio 2015 solution for Windows,
* MAKE project (G++ 4.8 required) for Linux and OS X.

At the top of the makefile there are several switches controlling building process. These are:
* STATIC_LINK - enable static linking (default: false); may be helpful when binary portability is desired,
* NO_AVX - prevent from using AVX and AVX2 extensions (default: false),
* NO_AVX2 - prevent from using AVX2 extensions (default: false),

Note, that FAMSA by default takes advantage of AVX and AVX2 CPU extensions. Pre-built binary detetermines supported instructions at runtime, thus it is multiplatform. However, one may encounter a problem when building FAMSA version on a CPU without AVX and/or AVX2. For this purpose NO_AVX and NO_AVX2 switches are provided.

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

Advanced options:
* `-r <value>` - no. of refinement iterations (default: 100)
* `-fr` - force refinement (by default the refinement is disabled for sets larger than 1000 seq.)
* `-go <value>` - gap open cost (default: 14.85)
* `-ge <value>` - gap extension cost (default: 1.25)
* `-tgo <value>` - terminal gap open cost (default: 0.66)
* `-tge <value>` - terminal gap extenstion cost (default: 0.66)
* `-gsd <value>` - gap cost scaller div-term (default: 7)
* `-gsl <value>` - gap cost scaller log-term (default: 45)
* `-dgr` - disable gap cost rescaling (default: enabled)
* `-dgo` - disable gap optimization (default: enabled)
* `-dsp` - disable sum of pairs optimization during refinement (default: enabled)	

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


## Citing
[Deorowicz, S., Debudaj-Grabysz, A., Gudyś, A. (2016) FAMSA: Fast and accurate multiple sequence alignment of huge protein families. 
Scientific Reports, 6, 33964](https://www.nature.com/articles/srep33964)




