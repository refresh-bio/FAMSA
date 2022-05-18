/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#ifndef _META_OPER_H
#define _META_OPER_H

//#include <functional>

//#define UNROLL(N, x)	((N > 0) ? ((x), UNROLL((N)-1, (x))) : (x))

template <size_t N> struct uint_{ };

// For loop (forward)
template <size_t N, typename Lambda>
static inline void IterFwd(const Lambda &oper, uint_<N>) {
	IterFwd(oper, uint_<N-1>());
	oper(N);
}

template <typename Lambda>
static inline void IterFwd(const Lambda &oper, uint_<0>) {
	oper(0);
}

// For loop (backward)
template <size_t N, typename Lambda>
inline void IterRev(const Lambda &oper, uint_<N>) {
	oper(N);
	IterRev(oper, uint_<N-1>());
}

template <typename Lambda>
inline void IterRev(const Lambda &oper, uint_<0>) {
	oper(0);
}

#ifdef _MSC_VER					// Visual C++
#include <intrin.h>
#define POPCNT(x)	(uint32_t) __popcnt64(x)
#endif

#ifdef __GNUC__
#define POPCNT(x)	(uint32_t) __builtin_popcountll(x)
#endif


#endif

// ***** EOF
