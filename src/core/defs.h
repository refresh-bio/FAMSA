/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#ifndef _DEFS_H
#define _DEFS_H

#include <cstdlib>
#include <cstdint>

// Uncomment for huge alignments (e.g., no of sequences > 10^6 and final alignemnt length > 10^5), when
// the default int64_t type for storing alignment scores could be too small
//#define HUGE_ALIGNMENTS

// Uncommenting enables additional parameters that were important for a development of FAMSA.
// They were also used to obtain some results presented in the FAMSA paper.
// Nevertheless, they are unimportant for end users, so avoid turning this mode on unless you are 
// really sure what you want to do. 
// Warning: FAMSA was not designed to be user-friendly in the developer mode.
//#define DEVELOPER_MODE


// ***** Internal defines
//#define DEBUG_MODE

#define ALWAYS_3_DIRS

//#define NO_GAP_CORRECTION

#define LOG_STATS

#ifdef HUGE_ALIGNMENTS
typedef double score_t;
#else
typedef int64_t score_t;
#endif

#define SIMD_NONE 0
#define SIMD_AVX1 1
#define SIMD_AVX2 2
#define SIMD_AVX512 3
#define SIMD_NEON 4

enum class instruction_set_t { none, sse, sse2, sse3, sse3s, sse41, sse42, avx, avx2 };

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


const symbol_t GAP			  = 30;			// value representing gap
const symbol_t GAP_OPEN       = 25;
const symbol_t GAP_EXT        = 26;
const symbol_t GAP_TERM_EXT   = 27;
const symbol_t GAP_TERM_OPEN  = 28;
const symbol_t UNKNOWN_SYMBOL = 22;

const size_t   NO_SYMBOLS				= 32;			// alphabet of protein sequences (including gaps and special symbols)
const symbol_t GUARD					= (symbol_t) (NO_SYMBOLS) - 1;
const symbol_t NO_AMINOACIDS			= 24;
const symbol_t NO_VALID_AMINOACIDS		= 20;
const symbol_t NO_AMINOACIDS_AND_GAPS	= 30;
const symbol_t NO_AA_SYMBOLS			= UNKNOWN_SYMBOL;		// no. of symbols that can be compared

#define MAX3(x, y, z)		(max(x, max(y, z)))
#define ABS(x)				((x) >= 0 ? (x) : -(x))


inline void *my_align(std::size_t alignment, std::size_t size,
	void *&ptr, std::size_t &space) {
	std::uintptr_t pn = reinterpret_cast< std::uintptr_t >(ptr);
	std::uintptr_t aligned = (pn + alignment - 1) & -alignment;
	std::size_t padding = aligned - pn;
	if (space < size + padding) return nullptr;
	space -= padding;
	return ptr = reinterpret_cast< void * >(aligned);
}

#if ((defined(_MSVC_LANG) && _MSVC_LANG >= 201703L) || __cplusplus >= 201703L)
#define FALL_THROUGH [[fallthrough]];
#else
#define FALL_THROUGH
#endif

#endif