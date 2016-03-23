/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

Version: 1.0
Date   : 2016-03-11
*/

#ifndef _DEFS_H
#define _DEFS_H

#include <cstdlib>
#include <cstdint>

#define FAMSA_VER		"1.0"
#define FAMSA_DATE		"2016-03-11"
#define FAMSA_AUTHORS	"S. Deorowicz, A. Debudaj-Grabysz, A. Gudys"

// Uncomment for huge alignments (e.g., no of sequences > 10^6 and final alignemnt length > 10^5), when
// the default int64_t type for storing alignment scores could be too small
//#define HUGE_ALIGNMENTS

//#define DEBUG_MODE

#define ALWAYS_3_DIRS

//#define NO_GAP_CORRECTION

#ifdef HUGE_ALIGNMENTS
typedef double score_t;
#else
typedef int64_t score_t;
#endif

typedef char symbol_t;
typedef int counter_t;

typedef unsigned long long bit_vec_t;

const int bv_size = sizeof(bit_vec_t) * 8;
const int bv_size128 = 64;					// length of a single word in AVX type used for bit-par LCS computation
const int bv_size256 = 64;					// length of a single word in AVX2 type used for bit-par LCS computation

#ifdef HUGE_ALIGNMENTS
const score_t infty = 1e30;
#else
const score_t infty = (1ll << 62);
const double cost_cast_factor = 1000.0;
#endif

const double infty_double = 1e30;

const symbol_t GAP			  = 30;			// value representing gap
const symbol_t GAP_OPEN       = 25;
const symbol_t GAP_EXT        = 26;
const symbol_t GAP_TERM_EXT   = 27;
const symbol_t GAP_TERM_OPEN  = 28;
const symbol_t UNKNOWN_SYMBOL = 23;

const size_t   NO_SYMBOLS				= 32;			// alphabet of preotein sequences (including gaps and special symbols)
const symbol_t GUARD					= (symbol_t) (NO_SYMBOLS) - 1;
const symbol_t NO_AMINOACIDS			= 24;
const symbol_t NO_VALID_AMINOACIDS		= 20;
const symbol_t NO_AMINOACIDS_AND_GAPS	= 30;

#endif