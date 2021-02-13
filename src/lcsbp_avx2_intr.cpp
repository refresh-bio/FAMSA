/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#include "lcsbp_avx2_intr.h"
#include "defs.h"

#include <algorithm>
#include <memory>
#include <immintrin.h>

using namespace std;

// *******************************************************************
// Prepares (if necessary sufficient amount of memory for LCS calculation
void CLCSBP_AVX2_INTR::prepare_X(size_t bv_len)
{
	size_t new_X_size = bv_len * sizeof(__m256i);

	if (new_X_size <= X_size)
		return;

	if (orig_X)
		free(orig_X);

	X_size = new_X_size;
	raw_X_size = X_size + 64;
	raw_X = malloc(raw_X_size);
	orig_X = raw_X;

	X = (__m256i*)my_align(64, X_size, raw_X, raw_X_size);
}

// *******************************************************************
void CLCSBP_AVX2_INTR::prepare_mask_pairs(size_t bv_len, CSequence* seq0)
{
	if (seq0_prev == seq0)
		return;

	if (raw_precomp_masks)
		free(raw_precomp_masks);

	seq0_prev = seq0;

	size_precomp_masks = bv_len * sizeof(__m128i) * NO_SYMBOLS * NO_SYMBOLS;
	size_t raw_size_precomp_masks = size_precomp_masks + 64;

	auto p = raw_precomp_masks = malloc(raw_size_precomp_masks);

	precomp_masks = (__m128i*) my_align(64, size_precomp_masks, p, raw_size_precomp_masks);

	const Array<bit_vec_t>& bit_masks = *(seq0->bit_masks);
	uint64_t* bm = (uint64_t*)bit_masks[0];
	uint64_t bm_len = bit_masks.get_width();

	for (uint32_t i = 0; i < NO_SYMBOLS; ++i)
		for (uint32_t j = 0; j < NO_SYMBOLS; ++j)
		{
			__m128i* ptr = precomp_masks + i * NO_SYMBOLS * bv_len + j * bv_len;
			uint64_t* pbm1 = bm + i * bm_len;
			uint64_t* pbm2 = bm + j * bm_len;

			for (uint32_t k = 0; k < bv_len; ++k)
				*ptr++ = _mm_set_epi64x(*pbm2++, *pbm1++);
		}
}

// *******************************************************************
void CLCSBP_AVX2_INTR::calculate(CSequence* seq0, CSequence* seq1, CSequence* seq2, CSequence* seq3, CSequence* seq4, uint32_t* res, uint32_t bv_len, uint32_t max_len)
{
	__m256i V, tB, V2, sB, U1, U2;
	__m256i sign64_bit = _mm256_set1_epi64x(1ull << 63);
	__m256i ones = _mm256_set1_epi64x(~0ull);

	const Array<bit_vec_t>& bit_masks = *(seq0->bit_masks);

	for (size_t i = 0; i < bv_len; ++i)
		X[i] = ones;

	//	__int64* base = (__int64*) bit_masks[0];

	//	__m128i plus_one = _mm_set1_epi32(1);
	int bm_width = bit_masks.get_width();

	for (size_t i = 0; i < max_len; ++i)
	{
		sB = _mm256_setzero_si256();
		symbol_t c1 = seq1->data[i];
		symbol_t c2 = seq2->data[i];
		symbol_t c3 = seq3->data[i];
		symbol_t c4 = seq4->data[i];

		//		__m128i indices = _mm_set_epi32(c4 * bm_width, c3 * bm_width, c2 * bm_width, c1 * bm_width);

		for (size_t j = 0; j < bv_len; ++j)
		{
			V = X[j];
			//			U1 = _mm256_set_epi64x(bit_masks[c4][j], bit_masks[c3][j], bit_masks[c2][j], bit_masks[c1][j]);
			U2 = _mm256_set_m128i(precomp_masks[c3 * NO_SYMBOLS * bv_len + c4 * bv_len + j], precomp_masks[c1 * NO_SYMBOLS * bv_len + c2 * bv_len + j]);
			//			U2 = _mm256_i32gather_epi64(base, indices, 8);
//			indices = _mm_add_epi32(indices, plus_one);
//			tB = _mm256_and_si256(V, U1);
			tB = _mm256_and_si256(V, U2);
			//			tB = _mm256_and_si256(V, _mm256_set_epi64x(bit_masks[c4][j], bit_masks[c3][j], bit_masks[c2][j], bit_masks[c1][j]));
			V2 = _mm256_sub_epi64(_mm256_add_epi64(V, tB), sB);
			sB = _mm256_cmpgt_epi64(_mm256_xor_si256(V, sign64_bit), _mm256_xor_si256(V2, sign64_bit));
			X[j] = _mm256_or_si256(V2, _mm256_sub_epi64(V, tB));
		}
	}

	alignas(32) uint64_t p[4];
#ifdef _MSC_VER					// Visual C++

	for (size_t i = 0; i < bv_len; ++i)
	{
		_mm256_storeu_si256((__m256i*) p, X[i]);

		res[0] += POPCNT(~p[0]);
		res[1] += POPCNT(~p[1]);
		res[2] += POPCNT(~p[2]);
		res[3] += POPCNT(~p[3]);
	}
#else
#ifdef __GNUC__
	for (size_t i = 0; i < bv_len; ++i)
	{
		_mm256_storeu_si256((__m256i*) p, X[i]);

		res[0] += POPCNT(~p[0]);
		res[1] += POPCNT(~p[1]);
		res[2] += POPCNT(~p[2]);
		res[3] += POPCNT(~p[3]);
		/*		res[0] += __builtin_popcountll(~p[0]);
				res[1] += __builtin_popcountll(~p[1]);
				res[2] += __builtin_popcountll(~p[2]);
				res[3] += __builtin_popcountll(~p[3]);*/
	}
#else
	for (size_t i = 0; i < bv_len; ++i)
	{
		_mm256_storeu_si256((__m256i*) p, X[i]);
		for (int v = 0; v < 4; ++v)
			for (uint64_t T = ~p[v]; T; T &= T - 1)
				++res[v];
	}
#endif
#endif
}

// *******************************************************************
// AVX2 variant of the bit-parallel LCS len calculation (processes 4 pairs of sequences in parallel)
void CLCSBP_AVX2_INTR::Calculate(CSequence* seq0, CSequence* seq1, CSequence* seq2, CSequence* seq3, CSequence* seq4,
	uint32_t& dist1, uint32_t& dist2, uint32_t& dist3, uint32_t& dist4)
{
	size_t max_len;
	max_len = max(seq1->length, max(seq2->length, max(seq3->length, seq4->length)));

	size_t bv_len = (seq0->length + bv_size256 - 1) / bv_size256;

	prepare_X(bv_len);
	prepare_mask_pairs(bv_len, seq0);

	uint32_t res[4] = { 0 };

	switch (bv_len)
	{
	case 1:	CLCSBP_AVX2_INTR_Impl<1>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 2:	CLCSBP_AVX2_INTR_Impl<2>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 3:	CLCSBP_AVX2_INTR_Impl<3>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 4:	CLCSBP_AVX2_INTR_Impl<4>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 5:	CLCSBP_AVX2_INTR_Impl<5>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 6:	CLCSBP_AVX2_INTR_Impl<6>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 7:	CLCSBP_AVX2_INTR_Impl<7>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 8:	CLCSBP_AVX2_INTR_Impl<8>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 9:	CLCSBP_AVX2_INTR_Impl<9>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 10: CLCSBP_AVX2_INTR_Impl<10>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 11: CLCSBP_AVX2_INTR_Impl<11>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 12: CLCSBP_AVX2_INTR_Impl<12>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 13: CLCSBP_AVX2_INTR_Impl<13>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 14: CLCSBP_AVX2_INTR_Impl<14>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 15: CLCSBP_AVX2_INTR_Impl<15>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 16: CLCSBP_AVX2_INTR_Impl<16>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 17: CLCSBP_AVX2_INTR_Impl<17>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 18: CLCSBP_AVX2_INTR_Impl<18>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 19: CLCSBP_AVX2_INTR_Impl<19>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 20: CLCSBP_AVX2_INTR_Impl<20>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 21: CLCSBP_AVX2_INTR_Impl<21>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 22: CLCSBP_AVX2_INTR_Impl<22>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 23: CLCSBP_AVX2_INTR_Impl<23>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 24: CLCSBP_AVX2_INTR_Impl<24>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 25: CLCSBP_AVX2_INTR_Impl<25>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 26: CLCSBP_AVX2_INTR_Impl<26>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 27: CLCSBP_AVX2_INTR_Impl<27>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 28: CLCSBP_AVX2_INTR_Impl<28>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 29: CLCSBP_AVX2_INTR_Impl<29>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 30: CLCSBP_AVX2_INTR_Impl<30>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 31: CLCSBP_AVX2_INTR_Impl<31>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	case 32: CLCSBP_AVX2_INTR_Impl<32>::Calculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, res, max_len, X);					break;
	default: calculate(seq0, seq1, seq2, seq3, seq4, res, bv_len, max_len);		break;
	}

	//	calculate(seq0, seq1, seq2, seq3, seq4, res, bv_len, max_len);

	dist1 = res[0];
	dist2 = res[1];
	dist3 = res[2];
	dist4 = res[3];
}
