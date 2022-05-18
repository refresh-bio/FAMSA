/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#include "../lcs/lcsbp_classic.h"

// *******************************************************************
// Prepares (if necessary sufficient amount of memory for LCS calculation
void CLCSBP_Classic::prepare_X(uint32_t bv_len)
{
	if (bv_len > X_size)
	{
		if (X)
			delete[] X;

		X_size = bv_len;
		X = new bit_vec_t[X_size];
	}
}

// *******************************************************************
void CLCSBP_Classic::prefetch_bitmasks(CSequence *seq0)
{
	if (seq0 == pf_seq0)
		return;

	pf_seq0 = seq0;

	for (int i = 0; i < (int) NO_SYMBOLS; ++i) {
		s0bm[i] = seq0->p_bit_masks + i * seq0->p_bv_len;
	}
}

// *******************************************************************
void CLCSBP_Classic::Calculate(CSequence *seq0, CSequence *seq1,
	uint32_t *dist)
{
	uint32_t bv_len = (seq0->length + bv_size - 1) / bv_size;

	prepare_X(bv_len);
	prefetch_bitmasks(seq0);

	dist[0] = 0;

	switch (bv_len)
	{
	case 1:	CLCSBP_Classic_Impl<1, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 2:	CLCSBP_Classic_Impl<2, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 3:	CLCSBP_Classic_Impl<3, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 4:	CLCSBP_Classic_Impl<4, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 5:	CLCSBP_Classic_Impl<5, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 6:	CLCSBP_Classic_Impl<6, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 7:	CLCSBP_Classic_Impl<7, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 8:	CLCSBP_Classic_Impl<8, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 9:	CLCSBP_Classic_Impl<9, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 10:	CLCSBP_Classic_Impl<10, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 11:	CLCSBP_Classic_Impl<11, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 12:	CLCSBP_Classic_Impl<12, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 13:	CLCSBP_Classic_Impl<13, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 14:	CLCSBP_Classic_Impl<14, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 15:	CLCSBP_Classic_Impl<15, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 16:	CLCSBP_Classic_Impl<16, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 17:	CLCSBP_Classic_Impl<17, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 18:	CLCSBP_Classic_Impl<18, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 19:	CLCSBP_Classic_Impl<19, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 20:	CLCSBP_Classic_Impl<20, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 21:	CLCSBP_Classic_Impl<21, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 22:	CLCSBP_Classic_Impl<22, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 23:	CLCSBP_Classic_Impl<23, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 24:	CLCSBP_Classic_Impl<24, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 25:	CLCSBP_Classic_Impl<25, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 26:	CLCSBP_Classic_Impl<26, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 27:	CLCSBP_Classic_Impl<27, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 28:	CLCSBP_Classic_Impl<28, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 29:	CLCSBP_Classic_Impl<29, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 30:	CLCSBP_Classic_Impl<30, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 31:	CLCSBP_Classic_Impl<31, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 32:	CLCSBP_Classic_Impl<32, CSequence>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	default: 
		CLCSBP_Classic_Impl<1, CSequence>::LoopCalculate(seq0, seq1, dist, bv_len, X, s0bm);
	}
}

// *******************************************************************
void CLCSBP_Classic::Calculate(CSequence *seq0, CSequenceView *seq1,
	uint32_t *dist)
{
	uint32_t bv_len = (seq0->length + bv_size - 1) / bv_size;

	prepare_X(bv_len);
	prefetch_bitmasks(seq0);

	dist[0] = 0;

	switch (bv_len)
	{
	case 1:	CLCSBP_Classic_Impl<1, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 2:	CLCSBP_Classic_Impl<2, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 3:	CLCSBP_Classic_Impl<3, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 4:	CLCSBP_Classic_Impl<4, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 5:	CLCSBP_Classic_Impl<5, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 6:	CLCSBP_Classic_Impl<6, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 7:	CLCSBP_Classic_Impl<7, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 8:	CLCSBP_Classic_Impl<8, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 9:	CLCSBP_Classic_Impl<9, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 10:	CLCSBP_Classic_Impl<10, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 11:	CLCSBP_Classic_Impl<11, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 12:	CLCSBP_Classic_Impl<12, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 13:	CLCSBP_Classic_Impl<13, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 14:	CLCSBP_Classic_Impl<14, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 15:	CLCSBP_Classic_Impl<15, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 16:	CLCSBP_Classic_Impl<16, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 17:	CLCSBP_Classic_Impl<17, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 18:	CLCSBP_Classic_Impl<18, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 19:	CLCSBP_Classic_Impl<19, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 20:	CLCSBP_Classic_Impl<20, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 21:	CLCSBP_Classic_Impl<21, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 22:	CLCSBP_Classic_Impl<22, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 23:	CLCSBP_Classic_Impl<23, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 24:	CLCSBP_Classic_Impl<24, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 25:	CLCSBP_Classic_Impl<25, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 26:	CLCSBP_Classic_Impl<26, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 27:	CLCSBP_Classic_Impl<27, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 28:	CLCSBP_Classic_Impl<28, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 29:	CLCSBP_Classic_Impl<29, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 30:	CLCSBP_Classic_Impl<30, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 31:	CLCSBP_Classic_Impl<31, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	case 32:	CLCSBP_Classic_Impl<32, CSequenceView>::UnrolledCalculate(seq0, seq1, dist, X, s0bm);					break;
	default: 
		CLCSBP_Classic_Impl<1, CSequenceView>::LoopCalculate(seq0, seq1, dist, bv_len, X, s0bm);
	}
}
