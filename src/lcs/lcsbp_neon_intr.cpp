/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#ifdef __ARM_NEON

/* Neon instruction mappings:
https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html

https://developer.arm.com/architectures/instruction-sets/intrinsics/
https://developer.arm.com/documentation/den0018/a/NEON-Intrinsics/Constructing-a-vector-from-a-literal-bit-pattern
https://arm-software.github.io/acle/neon_intrinsics/advsimd.html


	__m128i						int64x2_t			uint64x2_t
	
	_mm_add_epi64				vaddq_s64			vaddq_u64
	
	_mm_and_si128				vandq_s64			vandq_u64
	
	_mm_cmpgt_epi64				vcgeq_s64			vcgeq_u64
	
X	_mm_cvtsi128_si32												// Copy the lower 32-bit integer in a to dst
			// Tylko do liczenia histogramu

	_mm_hadd_epi16				
	
X	_mm_loadu_si128													// Load 128-bits of integer data from memory into dst. mem_addr does not need to be aligned on any particular boundary.
		vld1q_dup_u16 - dla uint16x8_t

	_mm_min_epi16				vminq_s16
	
	_mm_or_si128				vorrq_s64			vorrq_u64
	
*	_mm_set_epi64x													// Set packed 64-bit integers in dst with the supplied values.
		vcombine_u64(vcreate_u64, vcreate_u64)

*	_mm_set1_epi64x				vdupq_n_s64			vdupq_n_u64		// Broadcast 64-bit integer a to all elements of dst. This intrinsic may generate the vpbroadcastq.
		
*	_mm_setzero_si128			vdupq_n_s64			vdupq_n_u64		// Return vector of type __m128i with all elements set to zero.
	
X	_mm_storeu_si128												// Store 128-bits of integer data from a into memory. mem_addr does not need to be aligned on any particular boundary.
			// Tylko do POPCNT - chyba zbêdne
	
	_mm_sub_epi64				vsubq_s64			vsubq_u64
	
	_mm_xor_si128				veorq_s64			veorq_u64

	__builtin_popcount			


*/

#include "lcsbp_neon_intr.h"
#include "../core/defs.h"

#include "algorithm"
#include <memory>

using namespace std;

// *******************************************************************
// Prepares (if necessary sufficient amount of memory for LCS calculation
void CLCSBP_NEON_INTR::prepare_X(uint32_t bv_len)
{
	uint32_t new_X_size = bv_len * sizeof(int64x2_t);

	if (new_X_size <= X_size)
		return;

	if (orig_X)
		free(orig_X);

	X_size = new_X_size;
	raw_X_size = X_size + 64;
	raw_X = malloc(raw_X_size);
	orig_X = raw_X;
	
	X = (int64x2_t*)my_align(64, X_size, raw_X, raw_X_size);
}

// *******************************************************************
void CLCSBP_NEON_INTR::calculate(CSequence* seq0, CSequence* seq1, CSequence* seq2, uint32_t* res, uint32_t bv_len, uint32_t max_len)
{
	int64x2_t V, tB, V2, sB;
	int64x2_t sign64_bit = vdupq_n_s64(1ull << 63);
	int64x2_t ones = vdupq_n_s64(~0ull);

	const Array<bit_vec_t>& bit_masks = seq0->bit_masks;

	for (size_t i = 0; i < bv_len; ++i)
		X[i] = ones;

	for (size_t i = 0; i < max_len; ++i)
	{
		sB = vdupq_n_s64(0);
		symbol_t c1 = seq1->data[i];
		symbol_t c2 = seq2->data[i];

		for (size_t j = 0; j < bv_len; ++j)
		{
			V = X[j];
			tB = vandq_s64(V, vcombine_s64(vcreate_s64(bit_masks[c1][j]), vcreate_s64(bit_masks[c2][j])));
			V2 = vsubq_s64(vaddq_s64(V, tB), sB);
			sB = (int64x2_t) vcgtq_s64(veorq_s64(V, sign64_bit), veorq_s64(V2, sign64_bit));
			X[j] = vorrq_s64(V2, vsubq_s64(V, tB));
		}
	}

	alignas(32) uint64_t p[2];
//#ifdef _MSC_VER					// Visual C++
	for (size_t i = 0; i < bv_len; ++i)
	{
		vst1q_s64((int64_t*) p, X[i]);

		res[0] += __builtin_popcountll(~p[0]);
		res[1] += __builtin_popcountll(~p[1]);
	}
/*#else
#ifdef __GNUC__
	for (size_t i = 0; i < bv_len; ++i)
	{
		_mm_storeu_si128((__m128i*) p, X[i]);

		res[0] += __builtin_popcount(~p[0]);
		res[1] += __builtin_popcount(~p[1]);
	}
#else
	for (size_t i = 0; i < bv_len; ++i)
	{
		_mm_storeu_si128((__m128i*) p, X[i]);

		for (int v = 0; v < 2; ++v)
			for (uint64_t T = ~p[v]; T; T &= T - 1)
				++res[v];
#endif
#endif*/
}

// *******************************************************************
// AVX variant of the bit-parallel LCS len calculation (processes 2 pairs of sequences in parallel)
void CLCSBP_NEON_INTR::Calculate(CSequence *seq0, CSequence *seq1, CSequence *seq2,
	uint32_t *dist)
{
	uint32_t max_len;
	max_len = max(seq1->length, seq2->length);

	uint32_t bv_len = (seq0->length + bv_size128 - 1) / bv_size128;

	prepare_X(bv_len);

	dist[0] = dist[1] = 0;

	switch (bv_len)
	{
	case 1:	CLCSBP_AVX_INTR_Impl<1>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 2: CLCSBP_AVX_INTR_Impl<2>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 3: CLCSBP_AVX_INTR_Impl<3>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 4: CLCSBP_AVX_INTR_Impl<4>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 5: CLCSBP_AVX_INTR_Impl<5>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 6: CLCSBP_AVX_INTR_Impl<6>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 7: CLCSBP_AVX_INTR_Impl<7>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 8: CLCSBP_AVX_INTR_Impl<8>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 9: CLCSBP_AVX_INTR_Impl<9>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 10: CLCSBP_AVX_INTR_Impl<10>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 11: CLCSBP_AVX_INTR_Impl<11>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 12: CLCSBP_AVX_INTR_Impl<12>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 13: CLCSBP_AVX_INTR_Impl<13>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 14: CLCSBP_AVX_INTR_Impl<14>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 15: CLCSBP_AVX_INTR_Impl<15>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 16: CLCSBP_AVX_INTR_Impl<16>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 17: CLCSBP_AVX_INTR_Impl<17>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 18: CLCSBP_AVX_INTR_Impl<18>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 19: CLCSBP_AVX_INTR_Impl<19>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 20: CLCSBP_AVX_INTR_Impl<20>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 21: CLCSBP_AVX_INTR_Impl<21>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 22: CLCSBP_AVX_INTR_Impl<22>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 23: CLCSBP_AVX_INTR_Impl<23>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 24: CLCSBP_AVX_INTR_Impl<24>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 25: CLCSBP_AVX_INTR_Impl<25>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 26: CLCSBP_AVX_INTR_Impl<26>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 27: CLCSBP_AVX_INTR_Impl<27>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 28: CLCSBP_AVX_INTR_Impl<28>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 29: CLCSBP_AVX_INTR_Impl<29>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 30: CLCSBP_AVX_INTR_Impl<30>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 31: CLCSBP_AVX_INTR_Impl<31>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 32: CLCSBP_AVX_INTR_Impl<32>::Calculate(seq0, seq1, seq2, dist, max_len, X);					break;
	default: calculate(seq0, seq1, seq2, dist, bv_len, max_len);		break;
	}
}

// *******************************************************************
uint32_t CLCSBP_NEON_INTR::HistogramLCS(const uint16_t* h0, const uint16_t* h1)
{
	// Currently no SIMD here

	uint32_t est_lcs = 0;

	for (uint32_t i = 0; i < NO_SYMBOLS; ++i)
		est_lcs += min(h0[i], h1[i]);

	return est_lcs;
}
#endif