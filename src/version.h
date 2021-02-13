/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#ifndef _VERSION_H
#define _VERSION_H

#define FAMSA_VER		"1.6.2"
#define FAMSA_DATE		"2020-06-19"
#define FAMSA_AUTHORS	"S. Deorowicz, A. Debudaj-Grabysz, A. Gudys"

#endif

/*
Version history:

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
