/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#include "../core/lcsbp_avx.h"
#include "../core/defs.h"

#include "algorithm"
#include <memory>

using namespace std;

// *******************************************************************
// Prepares (if necessary sufficient amount of memory for LCS calculation
void CLCSBP_AVX::prepare_X(size_t bv_len)
{
	size_t new_X_size = bv_len * sizeof(Vec2uq);

	if (new_X_size <= X_size)
		return;

	if (orig_X)
		free(orig_X);

	X_size = new_X_size;
	raw_X_size = X_size + 64;
	raw_X = malloc(raw_X_size);
	orig_X = raw_X;
	
	X = (Vec2uq*)my_align(64, X_size, raw_X, raw_X_size);
}

// *******************************************************************
void CLCSBP_AVX::calculate(CSequence *seq0, CSequence *seq1, CSequence *seq2, uint32_t *res, uint32_t bv_len, uint32_t max_len)
{
	Vec2uq V, tB, V2, sB;

	for (size_t i = 0; i < bv_len; ++i)
		X[i] = Vec2uq(~(uint64_t)0);

	for (size_t i = 0; i < max_len; ++i)
	{
		sB = 0;
		symbol_t c1 = seq1->data[i];
		symbol_t c2 = seq2->data[i];

		for (size_t j = 0; j < bv_len; ++j)
		{
			V = X[j];
			tB = V & Vec2uq(seq0->bit_masks[c1][j], seq0->bit_masks[c2][j]);
			V2 = V + tB - sB;			// Important: true in boolean vectors is represented as -1, so here we use "-sB" instead of "+sB" in the classic code
			sB = V2 < V;
			X[j] = V2 | (V - tB);
		}
	}

#ifdef _MSC_VER					// Visual C++
	for (size_t i = 0; i < bv_len; ++i)
	{
		res[0] += __popcnt64(~X[i][0]);
		res[1] += __popcnt64(~X[i][1]);
	}
#else
#ifdef __GNUC__
	for (size_t i = 0; i < bv_len; ++i)
	{
		res[0] += __builtin_popcountll(~X[i][0]);
		res[1] += __builtin_popcountll(~X[i][1]);
	}
#else
	for (size_t i = 0; i < bv_len; ++i)
		for (int v = 0; v < 2; ++v)
			for (uint64_t T = ~X[i][v]; T; T &= T - 1)
				++res[v];
#endif
#endif
}

// *******************************************************************
// AVX variant of the bit-parallel LCS len calculation (processes 2 pairs of sequences in parallel)
void CLCSBP_AVX::Calculate(CSequence *seq0, CSequence *seq1, CSequence *seq2,
	uint32_t &dist1, uint32_t &dist2)
{
	size_t max_len;
	max_len = max(seq1->length, seq2->length);

	size_t bv_len = (seq0->length + bv_size128 - 1) / bv_size128;

	prepare_X(bv_len);

	uint32_t res[2] = { 0 };

	switch (bv_len)
	{
	case 1:	CLCSBP_AVX_Impl<1>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 2: CLCSBP_AVX_Impl<2>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 3: CLCSBP_AVX_Impl<3>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 4: CLCSBP_AVX_Impl<4>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 5: CLCSBP_AVX_Impl<5>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 6: CLCSBP_AVX_Impl<6>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 7: CLCSBP_AVX_Impl<7>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 8: CLCSBP_AVX_Impl<8>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 9: CLCSBP_AVX_Impl<9>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 10: CLCSBP_AVX_Impl<10>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 11: CLCSBP_AVX_Impl<11>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 12: CLCSBP_AVX_Impl<12>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 13: CLCSBP_AVX_Impl<13>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 14: CLCSBP_AVX_Impl<14>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 15: CLCSBP_AVX_Impl<15>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 16: CLCSBP_AVX_Impl<16>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 17: CLCSBP_AVX_Impl<17>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 18: CLCSBP_AVX_Impl<18>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 19: CLCSBP_AVX_Impl<19>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 20: CLCSBP_AVX_Impl<20>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 21: CLCSBP_AVX_Impl<21>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 22: CLCSBP_AVX_Impl<22>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 23: CLCSBP_AVX_Impl<23>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 24: CLCSBP_AVX_Impl<24>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 25: CLCSBP_AVX_Impl<25>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 26: CLCSBP_AVX_Impl<26>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 27: CLCSBP_AVX_Impl<27>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 28: CLCSBP_AVX_Impl<28>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 29: CLCSBP_AVX_Impl<29>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 30: CLCSBP_AVX_Impl<30>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 31: CLCSBP_AVX_Impl<31>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 32: CLCSBP_AVX_Impl<32>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	default: calculate(seq0, seq1, seq2, res, bv_len, max_len);		break;
	}
	
	dist1 = res[0];
	dist2 = res[1];
}