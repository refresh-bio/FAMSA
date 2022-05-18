/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#include <cstring>
#include "../core/defs.h"

#if SIMD==SIMD_AVX2 || SIMD==SIMD_AVX512
void mem_clear_avx2(void* ptr, size_t size)
{
	memset(ptr, 0, size);
}
#endif


