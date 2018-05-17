# FAMSA
Algorithm for large-scale multiple sequence alignments (400k proteins in 2 hours and 8BG of RAM)

## Installation and configuration

FAMSA comes with a set of precompiled binaries for Windows, Linux, and OS X. They can be found under Releases tab. 
The basic variant of the algorithm which executes entirely on CPU is contained in a file named *famsa*. 
For Windows and Linux systems there is an additional executable *famsa-gpu* which employs massively parallel processing (GPGPU)
with a use of OpenCL.

The software can be also built from the sources distributed as:

* Visual Studio 2015 solution for Windows,
* MAKE project (G++ 4.8 required) for Linux and OS X.

At the top of the makefile there are several switches controlling building process. These are:
* STATIC_LINK - enable static linking (default: false); may be helpful when binary portability is desired,
* COMPILE_GPU - build *famsa-gpu* binary (default: true); disable it when encounter OpenCL linking problems.

## Usage

`famsa [options] <input_file_name> <output_file_name>`

Positional parameters:
* `input_file_name` - input file in FASTA format or STDIN when reading from standard input
* `output_file_name` - output file in FASTA format or STDOUT when writing to standard output

Options:
* `-go <value>` - gap open cost (default: 14.85)
* `-ge <value>` - gap extension cost (default: 1.25)
* `-tgo <value>` - terminal gap open cost (default: 0.66)
* `-tge <value>` - terminal gap extenstion cost (default: 0.66)
* `-gsd <value>` - gap cost scaller div-term (default: 7)
* `-gsl <value>` - gap cost scaller log-term (default: 45)
* `-dgr` - disable gap cost rescaling (default: enabled)
* `-dgo` - disable gap optimization (default: enabled)
* `-dsp` - disable sum of pairs optimization during refinement (default: enabled)
* `-r <value>` - no. of refinement iterations (default: 100)
* `-fr` - disable auto refinement turning off (for sets larger than 1000 seq.)
* `-t <value>` - no. of threads, 0 means all available (default: 0)
* `-v` - verbose mode, show timing information (default: disabled)
* `-gt <sl, upgma, chained>` - choice of guide tree method: single linkage, UPGMA, chained (default: sl)
* `-gt_import <file_name>` - import guide tree in Newick format

When running *famsa-gpu* executable, two additional parameters must be specified:
* `-gpu_p <value>` - gpu platform id
* `-gpu_d <value>` - gpu device id

## Citing
[Deorowicz, S., Debudaj-Grabysz, A., Gudyś, A. (2016) FAMSA: Fast and accurate multiple sequence alignment of huge protein families. 
Scientific reports, 6, 33964](https://www.nature.com/articles/srep33964)




