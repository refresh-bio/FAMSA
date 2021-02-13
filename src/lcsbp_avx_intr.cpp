/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#include "lcsbp_avx_intr.h"
#include "defs.h"

#include "algorithm"
#include <memory>

using namespace std;

// *******************************************************************
// Prepares (if necessary sufficient amount of memory for LCS calculation
void CLCSBP_AVX_INTR::prepare_X(size_t bv_len)
{
	size_t new_X_size = bv_len * sizeof(__m128i);

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
void CLCSBP_AVX_INTR::calculate(CSequence* seq0, CSequence* seq1, CSequence* seq2, uint32_t* res, uint32_t bv_len, uint32_t max_len)
{
	__m128i V, tB, V2, sB;
	__m128i sign64_bit = _mm_set1_epi64x(1ull << 63);
	__m128i ones = _mm_set1_epi64x(~0ull);

	const Array<bit_vec_t>& bit_masks = *(seq0->bit_masks);

	for (size_t i = 0; i < bv_len; ++i)
		X[i] = ones;

	for (size_t i = 0; i < max_len; ++i)
	{
		sB = _mm_setzero_si128();
		symbol_t c1 = seq1->data[i];
		symbol_t c2 = seq2->data[i];

		for (size_t j = 0; j < bv_len; ++j)
		{
			V = X[j];
			tB = _mm_and_si128(V, _mm_set_epi64x(bit_masks[c2][j], bit_masks[c1][j]));
			V2 = _mm_sub_epi64(_mm_add_epi64(V, tB), sB);
			sB = _mm_cmpgt_epi64(_mm_xor_si128(V, sign64_bit), _mm_xor_si128(V2, sign64_bit));
			X[j] = _mm_or_si128(V2, _mm_sub_epi64(V, tB));
		}
	}

	alignas(32) uint64_t p[2];
#ifdef _MSC_VER					// Visual C++
	for (size_t i = 0; i < bv_len; ++i)
	{
		_mm_storeu_si128((__m128i*) p, X[i]);

		res[0] += POPCNT(~p[0]);
		res[1] += POPCNT(~p[1]);
	}
#else
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
#endif
}

// *******************************************************************
// AVX variant of the bit-parallel LCS len calculation (processes 2 pairs of sequences in parallel)
void CLCSBP_AVX_INTR::Calculate(CSequence *seq0, CSequence *seq1, CSequence *seq2,
	uint32_t &dist1, uint32_t &dist2)
{
	size_t max_len;
	max_len = max(seq1->length, seq2->length);

	size_t bv_len = (seq0->length + bv_size128 - 1) / bv_size128;

	prepare_X(bv_len);

	uint32_t res[2] = { 0 };

	switch (bv_len)
	{
	case 1:	CLCSBP_AVX_INTR_Impl<1>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 2: CLCSBP_AVX_INTR_Impl<2>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 3: CLCSBP_AVX_INTR_Impl<3>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 4: CLCSBP_AVX_INTR_Impl<4>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 5: CLCSBP_AVX_INTR_Impl<5>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 6: CLCSBP_AVX_INTR_Impl<6>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 7: CLCSBP_AVX_INTR_Impl<7>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 8: CLCSBP_AVX_INTR_Impl<8>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 9: CLCSBP_AVX_INTR_Impl<9>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 10: CLCSBP_AVX_INTR_Impl<10>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 11: CLCSBP_AVX_INTR_Impl<11>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 12: CLCSBP_AVX_INTR_Impl<12>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 13: CLCSBP_AVX_INTR_Impl<13>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 14: CLCSBP_AVX_INTR_Impl<14>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 15: CLCSBP_AVX_INTR_Impl<15>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 16: CLCSBP_AVX_INTR_Impl<16>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 17: CLCSBP_AVX_INTR_Impl<17>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 18: CLCSBP_AVX_INTR_Impl<18>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 19: CLCSBP_AVX_INTR_Impl<19>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 20: CLCSBP_AVX_INTR_Impl<20>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 21: CLCSBP_AVX_INTR_Impl<21>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 22: CLCSBP_AVX_INTR_Impl<22>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 23: CLCSBP_AVX_INTR_Impl<23>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 24: CLCSBP_AVX_INTR_Impl<24>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 25: CLCSBP_AVX_INTR_Impl<25>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 26: CLCSBP_AVX_INTR_Impl<26>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 27: CLCSBP_AVX_INTR_Impl<27>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 28: CLCSBP_AVX_INTR_Impl<28>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 29: CLCSBP_AVX_INTR_Impl<29>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 30: CLCSBP_AVX_INTR_Impl<30>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 31: CLCSBP_AVX_INTR_Impl<31>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	case 32: CLCSBP_AVX_INTR_Impl<32>::Calculate(seq0, seq1, seq2, res, max_len, X);					break;
	default: calculate(seq0, seq1, seq2, res, bv_len, max_len);		break;
	}

	dist1 = res[0];
	dist2 = res[1];
}
