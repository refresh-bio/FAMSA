/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#ifndef _LCSBP_AVX_INTR_H
#define _LCSBP_AVX_INTR_H

#include "../core/sequence.h"
#include "../utils/meta_oper.h"
#include <nmmintrin.h>

template <unsigned BV_LEN, typename SeqType> class CLCSBP_AVX_INTR_Impl;

// *******************************************************************
//
class CLCSBP_AVX_INTR
{
public:
	void* raw_X;
	void* orig_X;
	uint32_t X_size;
	size_t raw_X_size;
	__m128i* X;

	inline void prepare_X(uint32_t bv_len);

public:
	CLCSBP_AVX_INTR() {
		X = nullptr;
		raw_X = nullptr;
		orig_X = nullptr;
		X_size = 0;
		raw_X_size = 0;
	};

	~CLCSBP_AVX_INTR()
	{
		if (orig_X)
			free(orig_X);
	}
	
	void Calculate(CSequence* seq0, CSequence* seq1, CSequence* seq2, uint32_t* dist);
	void Calculate(CSequence* seq0, CSequenceView* seq1, CSequenceView* seq2, uint32_t* dist);
};

// *******************************************************************
//
template <unsigned BV_LEN, typename SeqType> class CLCSBP_AVX_INTR_Impl {
public:
#define AVX_INTR_POP_CNT_LOOP					\
	{										\
		_mm_storeu_si128((__m128i*) p, *pX);		\
		res[0] += POPCNT(~p[0]);				\
		res[1] += POPCNT(~p[1]);				\
		++pX;								\
	}

	// Important: true in boolean vectors is represented as -1, so here we use "-sB" instead of "+sB" in the classic code
#define AVX_INTR_LCS_INNER_LOOP_FIRST					\
	{										\
		V = *pX;					\
		tB = _mm_and_si128(V, _mm_set_epi64x(*pbm2++, *pbm1++));				\
		V2 = _mm_add_epi64(V, tB);												\
		sB = _mm_cmpgt_epi64(_mm_xor_si128(V, sign64_bit), _mm_xor_si128(V2, sign64_bit));		\
		*pX++ = _mm_or_si128(V2, _mm_xor_si128(V, tB));											\
	}	
#define AVX_INTR_LCS_INNER_LOOP_SINGLE					\
	{										\
		V = *pX;					\
		tB = _mm_and_si128(V, _mm_set_epi64x(*pbm2++, *pbm1++));				\
		V2 = _mm_add_epi64(V, tB);												\
		*pX++ = _mm_or_si128(V2, _mm_xor_si128(V, tB));											\
	}	
#define AVX_INTR_LCS_INNER_LOOP					\
	{										\
		V = *pX;					\
		tB = _mm_and_si128(V, _mm_set_epi64x(*pbm2++, *pbm1++));				\
		V2 = _mm_sub_epi64(_mm_add_epi64(V, tB), sB);												\
		sB = _mm_cmpgt_epi64(_mm_xor_si128(V, sign64_bit), _mm_xor_si128(V2, sign64_bit));		\
		*pX++ = _mm_or_si128(V2, _mm_xor_si128(V, tB));											\
	}	
#define AVX_INTR_LCS_INNER_LOOP_LAST					\
	{										\
		V = *pX;					\
		tB = _mm_and_si128(V, _mm_set_epi64x(*pbm2++, *pbm1++));				\
		V2 = _mm_sub_epi64(_mm_add_epi64(V, tB), sB);												\
		*pX++ = _mm_or_si128(V2, _mm_xor_si128(V, tB));											\
	}	

	// *******************************************************************
	static void LoopCalculate(CSequence* seq0, SeqType* seq1, SeqType* seq2, uint32_t* res, uint32_t bv_len, uint32_t max_len, __m128i* X)
	{
		__m128i V, tB, V2, sB;
		__m128i sign64_bit = _mm_set1_epi64x(1ull << 63);
		__m128i ones = _mm_set1_epi64x(~0ull);

		bit_vec_t* bit_masks = seq0->p_bit_masks;
		uint32_t loc_bv_len = seq0->p_bv_len;

		auto pc1 = seq1->data;
		auto pc2 = seq2->data;

		for (size_t i = 0; i < loc_bv_len; ++i)
			X[i] = ones;

		for (size_t i = 0; i < max_len; ++i)
		{
			sB = _mm_setzero_si128();

			for (size_t j = 0; j < loc_bv_len; ++j)
			{
				V = X[j];
				tB = _mm_and_si128(V, _mm_set_epi64x(bit_masks[*pc2 * loc_bv_len + j], bit_masks[*pc1 * loc_bv_len + j]));
				V2 = _mm_sub_epi64(_mm_add_epi64(V, tB), sB);
				sB = _mm_cmpgt_epi64(_mm_xor_si128(V, sign64_bit), _mm_xor_si128(V2, sign64_bit));
				X[j] = _mm_or_si128(V2, _mm_sub_epi64(V, tB));
			}

			++pc1; ++pc2;
		}

		alignas(32) uint64_t p[2];
#ifdef _MSC_VER					// Visual C++
		for (size_t i = 0; i < loc_bv_len; ++i)
		{
			_mm_storeu_si128((__m128i*) p, X[i]);

			res[0] += POPCNT(~p[0]);
			res[1] += POPCNT(~p[1]);
		}
#else
#ifdef __GNUC__
		for (size_t i = 0; i < loc_bv_len; ++i)
		{
			_mm_storeu_si128((__m128i*) p, X[i]);

			res[0] += POPCNT(~p[0]);
			res[1] += POPCNT(~p[1]);
		}
#else
		for (size_t i = 0; i < loc_bv_len; ++i)
		{
			_mm_storeu_si128((__m128i*) p, X[i]);

			for (int v = 0; v < 2; ++v)
				for (uint64_t T = ~p[v]; T; T &= T - 1)
					++res[v];
		}
#endif
#endif
	}

	// *******************************************************************
	static void UnrolledCalculate(CSequence *seq0, SeqType* seq1, SeqType* seq2, uint32_t *res, uint32_t max_len, __m128i *X)
	{
		__m128i V, tB, V2, sB;
		__m128i sign64_bit = _mm_set1_epi64x(1ull << 63);
		__m128i ones = _mm_set1_epi64x(~0ull);

		auto pX0 = X;
		const bit_vec_t* bit_masks = seq0->p_bit_masks;

		if (BV_LEN > 0)				*pX0++ = ones;
		if (BV_LEN > 1)				*pX0++ = ones;
		if (BV_LEN > 2)				*pX0++ = ones;
		if (BV_LEN > 3)				*pX0++ = ones;
		if (BV_LEN > 4)				*pX0++ = ones;
		if (BV_LEN > 5)				*pX0++ = ones;
		if (BV_LEN > 6)				*pX0++ = ones;
		if (BV_LEN > 7)				*pX0++ = ones;
		if (BV_LEN > 8)				*pX0++ = ones;
		if (BV_LEN > 9)				*pX0++ = ones;
		if (BV_LEN > 10)			*pX0++ = ones;
		if (BV_LEN > 11)			*pX0++ = ones;
		if (BV_LEN > 12)			*pX0++ = ones;
		if (BV_LEN > 13)			*pX0++ = ones;
		if (BV_LEN > 14)			*pX0++ = ones;
		if (BV_LEN > 15)			*pX0++ = ones;
		if (BV_LEN > 16)			*pX0++ = ones;
		if (BV_LEN > 17)			*pX0++ = ones;
		if (BV_LEN > 18)			*pX0++ = ones;
		if (BV_LEN > 19)			*pX0++ = ones;
		if (BV_LEN > 20)			*pX0++ = ones;
		if (BV_LEN > 21)			*pX0++ = ones;
		if (BV_LEN > 22)			*pX0++ = ones;
		if (BV_LEN > 23)			*pX0++ = ones;
		if (BV_LEN > 24)			*pX0++ = ones;
		if (BV_LEN > 25)			*pX0++ = ones;
		if (BV_LEN > 26)			*pX0++ = ones;
		if (BV_LEN > 27)			*pX0++ = ones;
		if (BV_LEN > 28)			*pX0++ = ones;
		if (BV_LEN > 29)			*pX0++ = ones;
		if (BV_LEN > 30)			*pX0++ = ones;
		if (BV_LEN > 31)			*pX0++ = ones;

		uint64_t* bm = (uint64_t*)bit_masks;
		uint64_t bm_len = seq0->p_bv_len;
//		uint32_t bm_len = (seq0->length + bv_size128 - 1) / bv_size128;

		auto pc1 = seq1->data;
		auto pc2 = seq2->data;


		for (size_t i = 0; i < max_len; ++i)
		{
			sB = _mm_setzero_si128();

			auto pX = X;

			uint64_t* pbm1 = bm + *pc1 * bm_len;
			uint64_t* pbm2 = bm + *pc2 * bm_len;

			if (BV_LEN == 1)		AVX_INTR_LCS_INNER_LOOP_SINGLE
			else if (BV_LEN > 0)	AVX_INTR_LCS_INNER_LOOP_FIRST;
			if (BV_LEN == 2)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 1)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 3)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 2)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 4)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 3)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 5)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 4)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 6)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 5)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 7)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 6)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 8)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 7)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 9)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 8)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 10)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 9)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 11)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 10)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 12)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 11)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 13)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 12)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 14)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 13)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 15)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 14)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 16)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 15)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 17)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 16)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 18)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 17)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 19)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 18)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 20)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 19)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 21)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 20)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 22)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 21)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 23)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 22)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 24)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 23)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 25)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 24)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 26)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 25)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 27)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 26)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 28)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 27)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 29)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 28)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 30)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 29)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 31)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 30)	AVX_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 32)		AVX_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 31)	AVX_INTR_LCS_INNER_LOOP;

			++pc1; ++pc2;
		}

		alignas(32) uint64_t p[2];
		auto pX = X;
		if (BV_LEN > 0)		AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 1)		AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 2)		AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 3)		AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 4)		AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 5)		AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 6)		AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 7)		AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 8)		AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 9)		AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 10)	AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 11)	AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 12)	AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 13)	AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 14)	AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 15)	AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 16)	AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 17)	AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 18)	AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 19)	AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 20)	AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 21)	AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 22)	AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 23)	AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 24)	AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 25)	AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 26)	AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 27)	AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 28)	AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 29)	AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 30)	AVX_INTR_POP_CNT_LOOP;
		if (BV_LEN > 31)	AVX_INTR_POP_CNT_LOOP;
	}
};

#endif
