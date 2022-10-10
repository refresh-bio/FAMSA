/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#ifndef _VERSION_H
#define _VERSION_H

#define FAMSA_VER		"2.2.2"
#define FAMSA_DATE		"2022-10-09"
#define FAMSA_AUTHORS	"S. Deorowicz, A. Debudaj-Grabysz, A. Gudys"

#endif

/*
Version history:
2.2.2 (2022-10-09):
- Fixed slowdown caused by the duplicate removal (feature added in 2.1.0).

2.2.1 (2022-10-05):
- Pairwise identity (-pid switch) properly calculated as the number of matching residues divided by the shorter sequence length.

2.2.0 (2022-10-05):
- Added possibility to align two pre-aligned profiles (two input files specified).

2.1.3 (2022-09-30):
- Fixed incorrect handling of single sequence sets (first residue cut) or sets containing only duplicates (hang).

2.1.2 (2022-08-04):
- Makefile updates improving cross-platform compilation.

2.1.1 (2022-08-01):
- Preserving non-standard amino acid symbols in the output alignment (instead of replacing with X).

2.1.0 (2022-07-26):
- Duplicated sequences are removed prior to the alignment and restored afterwards.

2.0.4 (2022-07-25)
- Fixed stack overflow exception when saving/loading large Newick trees 
	(boost-based algorithm replaced with own non-recursive approach).
	
2.0.3 (2022-05-27):
- The ordering of the input sequences preserved in the final alignment.

2.0.2 (2022-05-24): 
- Alignment allowed as an input (gaps are removed).

2.0.1 (2022-05-18):
- Several fixes for bioconda.

2.0.0-rc (2022-05-18)::
- Default algorithm for single linkage trees changed from SLINK to MST Prim,
- Small fixes.

1.16.0 (2022-04-29):
 - mimalloc added as a single source library (no external linking).

1.15.0 (2022-04-28):
1.14.0 (2022-04-22):
 - Highly optimized MST Prim implementation. 

1.13.0 (2022-04-20):
 - Added support of NEON SIMD extensions for ARM architectures. 

1.12.5 (2022-04-13)
 - Further memory optimizations.
 - Some refactoring.
 - Warnings removal.

1.12.4 (2022-04-11)
1.12.3 (2022-04-11)
- Futher memory optimizations.

1.12.2 (2022-04-09)
- Memory optimizations in gapped sequence representation.

1.12.1 (2022-04-08)
- Non-AVX and 32-bit compilation fixed.
- Proper handling of . and | symbols in sequence identifiers when importing Newick tree.

1.12.0 (2022-04-06)
- Parallel medoid trees.

1.11.0 (2022-04-01)
- Gzipped output.

1.10.0 (2022-03-24)
- Uniform distance measures.

1.9.0 (2022-03-11)
- Export of distance matrix significantly improved.

1.8.0 (2021-06-09)
- Added MST Prim algorithm for single linkage trees.

1.7.0 (2021-06-08)
- Parallel profile construction.
- Multiple optimizations.
- Single linkage draws resolution.

1.6.2 (2020-06-19)
- Clang compilation fixed.

1.6.1 (2020-06-18)
- Added parameter for automatic medoid tree usage. Some refactoring. Added license file.

1.6.0 (2020-06-18)
- Removed VCL and ASMLIB dependencies.
- Some low-level optimizations in LCS calculations.

1.5.20 (2020-05-26)
- PartTree always select assumed number of seeds (1.5.16 patch applied only to MedoidTree). 
- Small fix in Sackin index calculation.

1.5.19 (2020-04-18)
- Uniform distance computation in PartTree and MedoidTree.

1.5.18 (2020-04-17)
- Added modified UPGMA algorithm and distance correction (MAFFT-inspired).

1.5.17 (2020-04-16)
- Stats dumped to file in the verbose mode.
- Sackin index calculation right after tree construction.

1.5.16 (2020-04-10)
- PartTree and MedoidTree always select assumed number of seeds.

1.5.15 (2020-04-08)
- Fixed bug in calculating clustering cost and assignment update.

1.5.14 (2020-04-07)
- Bug in Neighbor Joining fixed.

1.5.13 (2020-04-06)
- Alternative method of combining children and parental trees in PartTree.
- Neighbor Joining algorithm added.

1.5.12 (2020-03-24)
- Fixed bug with -dist_export mode (sequences not ordered as in input FASTA file)

1.5.11 (2020-03-23)
- Possiblity to choose between PartTree and MedoidTree heuristic.
- K-medoid clustering with CLARANS heuristic.

1.5.8 (2020-03-18)
- Deterministic random generator added (to make Windows and Unix results the same).

1.5.7 (2020-03-17)
- Further memory improvements.

1.5.6 (2020-03-15)
- Serious refactorization.
- Bit vectors computed when needed and released afterwards.

1.5.5 (2020-03-13)
- Segmentation fault fix.

1.5.4 (2020-03-10)
- Interface for tree import/export changed a bit.

1.5.3 (2020-03-03)
- PartTree + UPGMA support.
- Size of a cluster in PartTree added as a command line parameter.
- Optimizations in the sequence representation.
- Added option for sequence shuffling.

1.5.0 (2020-03-02)
- PartTree mode added (only single linkage trees supported).
- Since sequences are sorted, the first one can be taken as the longest one.
- Serious refactoring.

1.4.0 (2020-02-27)
- GPU mode no longer supported.
- Possibility to calculate distance matrix or guide tree without doing alignment.
- Lots of refactoring.

1.3.2 (2020-02-21)
- Approved pull request "Fix * char emitted for unknown residues, emits X instead"
- Approved pull request "Support single input fasta files"
- Version.h file added

*/
