/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#ifndef _UTILS_H
#define _UTILS_H

#include "../core/defs.h"
#include <vector>

#ifdef _MSC_VER	
#include <immintrin.h>
#endif

void mem_clear(void* ptr, size_t size);

#if defined(SIMD_AVX) || defined(SIMD_AVX2) || defined(SIMD_AVX512)
void mem_clear_avx(void* ptr, size_t size);
#endif

#if defined(SIMD_AVX2) || defined(SIMD_AVX512)
void mem_clear_avx2(void* ptr, size_t size);
#endif

#if defined(SIMD_NEON)
void mem_clear_neon(void* ptr, size_t size);
#endif

// *******************************************************************
template<typename T>
T max4(T x1, T x2, T x3, T x4) 
{
	T p1 = (x1 > x2) ? x1 : x2;
	T p2 = (x3 > x4) ? x3 : x4;
	
	return (p1 > p2) ? p1 : p2;
}

// *******************************************************************
template<typename T>
T max8(T x1, T x2, T x3, T x4, T x5, T x6, T x7, T x8) 
{
	T m1 = max4(x1, x2, x3, x4);
	T m2 = max4(x5, x6, x7, x8);

	return (m1 > m2) ? m1 : m2;
}

// *******************************************************************
template<typename T>
void clear_vector(std::vector<T>& vec)
{
	std::vector<T>().swap(vec);
}

// *******************************************************************
template<typename T>
void delete_arr_ptr(T* &ptr)
{
	if (!ptr)
		return;

	delete[] ptr;
	ptr = nullptr;
}

// *******************************************************************
template<typename T>
void delete_ptr(T* &ptr)
{
	if (!ptr)
		return;

	delete ptr;
	ptr = nullptr;
}

// *******************************************************************
template<typename T>
void tpl_prefetch(T* ptr)
{
#ifdef _MSC_VER					// Visual C++
	_mm_prefetch((const char*) ptr, 2);
#endif
#ifdef __GNUC__
	//			__builtin_prefetch((&(*dist_vector)[pi[j + prefetch_offset]]), 1, 2);
	__builtin_prefetch(ptr, 1, 2);
#endif
}

#endif