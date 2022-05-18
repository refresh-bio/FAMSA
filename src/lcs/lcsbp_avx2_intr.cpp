/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#include "lcsbp_avx2_intr.h"
#include "../core/defs.h"

#include <algorithm>
#include <memory>
#include <immintrin.h>

using namespace std;


// *******************************************************************
// Prepares (if necessary) sufficient amount of memory for LCS calculation
void CLCSBP_AVX2_INTR::prepare_X(uint32_t bv_len)
{
	uint32_t new_X_size = bv_len * sizeof(__m256i);

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
void CLCSBP_AVX2_INTR::prepare_mask_pairs(uint32_t bv_len, CSequence* seq0)
{
	if (seq0->sequence_no == prev_sequence_no)
		return;

	size_t new_size_precomp_masks = bv_len * sizeof(__m128i) * NO_SYMBOLS * NO_SYMBOLS;

	bool need_alloc = !raw_precomp_masks;

	if (raw_precomp_masks)
		need_alloc = new_size_precomp_masks > size_precomp_masks;

	if (need_alloc)
	{
		if (raw_precomp_masks)
			free(raw_precomp_masks);

		size_precomp_masks = new_size_precomp_masks;
		size_t raw_size_precomp_masks = size_precomp_masks + 64;

		auto p = raw_precomp_masks = malloc(raw_size_precomp_masks);
		precomp_masks = (__m128i*) my_align(64, size_precomp_masks, p, raw_size_precomp_masks);
	}

	prev_sequence_no = seq0->sequence_no;

	bit_vec_t * bit_masks = seq0->p_bit_masks;
	uint64_t* bm = (uint64_t*)bit_masks;
	uint64_t bm_len = seq0->p_bv_len;

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
void CLCSBP_AVX2_INTR::Calculate(CSequence* seq0, CSequence* seq1, CSequence* seq2, CSequence* seq3, CSequence* seq4,
	uint32_t* dist)
{
	uint32_t max_len = max4(seq1->length, seq2->length, seq3->length, seq4->length);
	uint32_t bv_len = (seq0->length + bv_size256 - 1) / bv_size256;

	prepare_X(bv_len);
	prepare_mask_pairs(bv_len, seq0);

	dist[0] = dist[1] = dist[2] = dist[3] = 0;

	switch (bv_len)
	{
	case 1:	CLCSBP_AVX2_INTR_Impl<1, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 2:	CLCSBP_AVX2_INTR_Impl<2, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 3:	CLCSBP_AVX2_INTR_Impl<3, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 4:	CLCSBP_AVX2_INTR_Impl<4, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 5:	CLCSBP_AVX2_INTR_Impl<5, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 6:	CLCSBP_AVX2_INTR_Impl<6, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 7:	CLCSBP_AVX2_INTR_Impl<7, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 8:	CLCSBP_AVX2_INTR_Impl<8, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 9:	CLCSBP_AVX2_INTR_Impl<9, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 10: CLCSBP_AVX2_INTR_Impl<10, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 11: CLCSBP_AVX2_INTR_Impl<11, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 12: CLCSBP_AVX2_INTR_Impl<12, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 13: CLCSBP_AVX2_INTR_Impl<13, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 14: CLCSBP_AVX2_INTR_Impl<14, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 15: CLCSBP_AVX2_INTR_Impl<15, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 16: CLCSBP_AVX2_INTR_Impl<16, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 17: CLCSBP_AVX2_INTR_Impl<17, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 18: CLCSBP_AVX2_INTR_Impl<18, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 19: CLCSBP_AVX2_INTR_Impl<19, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 20: CLCSBP_AVX2_INTR_Impl<20, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 21: CLCSBP_AVX2_INTR_Impl<21, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 22: CLCSBP_AVX2_INTR_Impl<22, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 23: CLCSBP_AVX2_INTR_Impl<23, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 24: CLCSBP_AVX2_INTR_Impl<24, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 25: CLCSBP_AVX2_INTR_Impl<25, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 26: CLCSBP_AVX2_INTR_Impl<26, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 27: CLCSBP_AVX2_INTR_Impl<27, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 28: CLCSBP_AVX2_INTR_Impl<28, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 29: CLCSBP_AVX2_INTR_Impl<29, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 30: CLCSBP_AVX2_INTR_Impl<30, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 31: CLCSBP_AVX2_INTR_Impl<31, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 32: CLCSBP_AVX2_INTR_Impl<32, CSequence>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	default: CLCSBP_AVX2_INTR_Impl<1, CSequence>::LoopCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, bv_len, max_len, X);		break;
	}
}

// *******************************************************************
void CLCSBP_AVX2_INTR::Calculate(CSequence* seq0, CSequenceView* seq1, CSequenceView* seq2, CSequenceView* seq3, CSequenceView* seq4,
	uint32_t* dist)
{
	uint32_t max_len = max4(seq1->length, seq2->length, seq3->length, seq4->length);
	uint32_t bv_len = (seq0->length + bv_size256 - 1) / bv_size256;

	prepare_X(bv_len);
	prepare_mask_pairs(bv_len, seq0);

	dist[0] = dist[1] = dist[2] = dist[3] = 0;

	switch (bv_len)
	{
	case 1:	CLCSBP_AVX2_INTR_Impl<1, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 2:	CLCSBP_AVX2_INTR_Impl<2, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 3:	CLCSBP_AVX2_INTR_Impl<3, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 4:	CLCSBP_AVX2_INTR_Impl<4, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 5:	CLCSBP_AVX2_INTR_Impl<5, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 6:	CLCSBP_AVX2_INTR_Impl<6, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 7:	CLCSBP_AVX2_INTR_Impl<7, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 8:	CLCSBP_AVX2_INTR_Impl<8, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 9:	CLCSBP_AVX2_INTR_Impl<9, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 10: CLCSBP_AVX2_INTR_Impl<10, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 11: CLCSBP_AVX2_INTR_Impl<11, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 12: CLCSBP_AVX2_INTR_Impl<12, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 13: CLCSBP_AVX2_INTR_Impl<13, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 14: CLCSBP_AVX2_INTR_Impl<14, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 15: CLCSBP_AVX2_INTR_Impl<15, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 16: CLCSBP_AVX2_INTR_Impl<16, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 17: CLCSBP_AVX2_INTR_Impl<17, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 18: CLCSBP_AVX2_INTR_Impl<18, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 19: CLCSBP_AVX2_INTR_Impl<19, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 20: CLCSBP_AVX2_INTR_Impl<20, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 21: CLCSBP_AVX2_INTR_Impl<21, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 22: CLCSBP_AVX2_INTR_Impl<22, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 23: CLCSBP_AVX2_INTR_Impl<23, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 24: CLCSBP_AVX2_INTR_Impl<24, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 25: CLCSBP_AVX2_INTR_Impl<25, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 26: CLCSBP_AVX2_INTR_Impl<26, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 27: CLCSBP_AVX2_INTR_Impl<27, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 28: CLCSBP_AVX2_INTR_Impl<28, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 29: CLCSBP_AVX2_INTR_Impl<29, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 30: CLCSBP_AVX2_INTR_Impl<30, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 31: CLCSBP_AVX2_INTR_Impl<31, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	case 32: CLCSBP_AVX2_INTR_Impl<32, CSequenceView>::UnrolledCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, max_len, X);					break;
	default: CLCSBP_AVX2_INTR_Impl<1, CSequenceView>::LoopCalculate(precomp_masks, seq0, seq1, seq2, seq3, seq4, dist, bv_len, max_len, X);		break;
	}
}

