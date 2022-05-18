/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#include "lcsbp_avx_intr.h"
#include "../core/defs.h"

#include "algorithm"
#include <memory>

using namespace std;

// *******************************************************************
// Prepares (if necessary sufficient amount of memory for LCS calculation
void CLCSBP_AVX_INTR::prepare_X(uint32_t bv_len)
{
	uint32_t new_X_size = bv_len * sizeof(__m128i);

	if (new_X_size <= X_size)
		return;

	if (orig_X)
		free(orig_X);

	X_size = new_X_size;
	raw_X_size = X_size + 64;
	raw_X = malloc(raw_X_size);
	orig_X = raw_X;

	X = (__m128i*)my_align(64, X_size, raw_X, raw_X_size);
}

// *******************************************************************
// AVX variant of the bit-parallel LCS len calculation (processes 2 pairs of sequences in parallel)
void CLCSBP_AVX_INTR::Calculate(CSequence *seq0, CSequence *seq1, CSequence *seq2,
	uint32_t *dist)
{
	uint32_t max_len = max(seq1->length, seq2->length);
//	uint32_t bv_len = (seq0->length + bv_size128 - 1) / bv_size128;
	uint32_t bv_len = seq0->p_bv_len;

	prepare_X(bv_len);

	dist[0] = dist[1] = 0;

	switch (bv_len)
	{
	case 1:	CLCSBP_AVX_INTR_Impl<1, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 2: CLCSBP_AVX_INTR_Impl<2, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 3: CLCSBP_AVX_INTR_Impl<3, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 4: CLCSBP_AVX_INTR_Impl<4, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 5: CLCSBP_AVX_INTR_Impl<5, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 6: CLCSBP_AVX_INTR_Impl<6, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 7: CLCSBP_AVX_INTR_Impl<7, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 8: CLCSBP_AVX_INTR_Impl<8, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 9: CLCSBP_AVX_INTR_Impl<9, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 10: CLCSBP_AVX_INTR_Impl<10, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 11: CLCSBP_AVX_INTR_Impl<11, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 12: CLCSBP_AVX_INTR_Impl<12, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 13: CLCSBP_AVX_INTR_Impl<13, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 14: CLCSBP_AVX_INTR_Impl<14, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 15: CLCSBP_AVX_INTR_Impl<15, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 16: CLCSBP_AVX_INTR_Impl<16, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 17: CLCSBP_AVX_INTR_Impl<17, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 18: CLCSBP_AVX_INTR_Impl<18, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 19: CLCSBP_AVX_INTR_Impl<19, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 20: CLCSBP_AVX_INTR_Impl<20, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 21: CLCSBP_AVX_INTR_Impl<21, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 22: CLCSBP_AVX_INTR_Impl<22, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 23: CLCSBP_AVX_INTR_Impl<23, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 24: CLCSBP_AVX_INTR_Impl<24, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 25: CLCSBP_AVX_INTR_Impl<25, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 26: CLCSBP_AVX_INTR_Impl<26, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 27: CLCSBP_AVX_INTR_Impl<27, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 28: CLCSBP_AVX_INTR_Impl<28, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 29: CLCSBP_AVX_INTR_Impl<29, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 30: CLCSBP_AVX_INTR_Impl<30, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 31: CLCSBP_AVX_INTR_Impl<31, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 32: CLCSBP_AVX_INTR_Impl<32, CSequence>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	default: CLCSBP_AVX_INTR_Impl<1, CSequence>::LoopCalculate(seq0, seq1, seq2, dist, bv_len, max_len, X);		break;
	}
}

// *******************************************************************
// AVX variant of the bit-parallel LCS len calculation (processes 2 pairs of sequences in parallel)
void CLCSBP_AVX_INTR::Calculate(CSequence *seq0, CSequenceView *seq1, CSequenceView *seq2,
	uint32_t *dist)
{
	uint32_t max_len = max(seq1->length, seq2->length);
//	uint32_t bv_len = (seq0->length + bv_size128 - 1) / bv_size128;
	uint32_t bv_len = seq0->p_bv_len;

	prepare_X(bv_len);

	dist[0] = dist[1] = 0;

	switch (bv_len)
	{
	case 1:	CLCSBP_AVX_INTR_Impl<1, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 2: CLCSBP_AVX_INTR_Impl<2, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 3: CLCSBP_AVX_INTR_Impl<3, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 4: CLCSBP_AVX_INTR_Impl<4, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 5: CLCSBP_AVX_INTR_Impl<5, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 6: CLCSBP_AVX_INTR_Impl<6, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 7: CLCSBP_AVX_INTR_Impl<7, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 8: CLCSBP_AVX_INTR_Impl<8, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 9: CLCSBP_AVX_INTR_Impl<9, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 10: CLCSBP_AVX_INTR_Impl<10, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 11: CLCSBP_AVX_INTR_Impl<11, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 12: CLCSBP_AVX_INTR_Impl<12, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 13: CLCSBP_AVX_INTR_Impl<13, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 14: CLCSBP_AVX_INTR_Impl<14, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 15: CLCSBP_AVX_INTR_Impl<15, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 16: CLCSBP_AVX_INTR_Impl<16, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 17: CLCSBP_AVX_INTR_Impl<17, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 18: CLCSBP_AVX_INTR_Impl<18, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 19: CLCSBP_AVX_INTR_Impl<19, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 20: CLCSBP_AVX_INTR_Impl<20, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 21: CLCSBP_AVX_INTR_Impl<21, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 22: CLCSBP_AVX_INTR_Impl<22, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 23: CLCSBP_AVX_INTR_Impl<23, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 24: CLCSBP_AVX_INTR_Impl<24, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 25: CLCSBP_AVX_INTR_Impl<25, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 26: CLCSBP_AVX_INTR_Impl<26, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 27: CLCSBP_AVX_INTR_Impl<27, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 28: CLCSBP_AVX_INTR_Impl<28, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 29: CLCSBP_AVX_INTR_Impl<29, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 30: CLCSBP_AVX_INTR_Impl<30, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 31: CLCSBP_AVX_INTR_Impl<31, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	case 32: CLCSBP_AVX_INTR_Impl<32, CSequenceView>::UnrolledCalculate(seq0, seq1, seq2, dist, max_len, X);					break;
	default: CLCSBP_AVX_INTR_Impl<1, CSequenceView>::LoopCalculate(seq0, seq1, seq2, dist, bv_len, max_len, X);		break;
	}
}

