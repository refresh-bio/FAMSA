/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#ifndef _LCSBP_AVX2_INTR_H
#define _LCSBP_AVX2_INTR_H

#include "../core/sequence.h"
#include "../utils/meta_oper.h"
#include "../utils/utils.h"
#include "../core/defs.h"

#include <algorithm>
#include <memory>
#include <immintrin.h>

template <unsigned BV_LEN, typename SeqType> class CLCSBP_AVX2_INTR_Impl;

// *******************************************************************
//
class CLCSBP_AVX2_INTR
{
	void* raw_X;
	void* orig_X;
	uint32_t X_size;
	size_t raw_X_size;
	__m256i* X;
	int prev_sequence_no;

	__m128i* precomp_masks;
	void* raw_precomp_masks;
	size_t size_precomp_masks;

	inline void prepare_X(uint32_t bv_len);
	inline void prepare_mask_pairs(uint32_t bv_len, CSequence* seq0);

public:
	CLCSBP_AVX2_INTR() {
		X = nullptr;
		raw_X = nullptr;
		orig_X = nullptr;
		X_size = 0;
		raw_X_size = 0;
		precomp_masks = nullptr;
		raw_precomp_masks = nullptr;
		size_precomp_masks = 0;

		prev_sequence_no = -1;
	};

	~CLCSBP_AVX2_INTR()
	{
		if (orig_X)
			free(orig_X);

		if (raw_precomp_masks)
			free(raw_precomp_masks);
	}

	void Calculate(CSequence* seq0, CSequence* seq1, CSequence* seq2, CSequence* seq3, CSequence* seq4,
		uint32_t* dist);
	void Calculate(CSequence* seq0, CSequenceView* seq1, CSequenceView* seq2, CSequenceView* seq3, CSequenceView* seq4,
		uint32_t* dist);
};

//#undef _mm256_set_m128i
#ifndef _mm256_set_m128i
#define _mm256_set_m128i(v0, v1)  _mm256_insertf128_si256(_mm256_castsi128_si256(v1), (v0), 1)
#endif

// *******************************************************************
//
template <unsigned BV_LEN, typename SeqType> class CLCSBP_AVX2_INTR_Impl {
public:
#define AVX2_INTR_POP_CNT_LOOP						\
	{												\
		_mm256_storeu_si256((__m256i*) p, *pX);		\
		res[0] += POPCNT(~p[0]);				\
		res[1] += POPCNT(~p[1]);				\
		res[2] += POPCNT(~p[2]);				\
		res[3] += POPCNT(~p[3]);				\
		++pX;										\
	}

	// Important: true in boolean vectors is represented as -1, so here we use "-sB" instead of "+sB" in the classic code
/*#define AVX2_INTR_LCS_INNER_LOOP_FIRST										\
	{															\
		V = *pX;		\
		tB = _mm256_and_si256(V, _mm256_set_epi64x(*pbm4++, *pbm3++, *pbm2++, *pbm1++));				\
		V2 = _mm256_add_epi64(V, tB);												\
		sB = _mm256_cmpgt_epi64(_mm256_xor_si256(V, sign64_bit), _mm256_xor_si256(V2, sign64_bit));		\
		*pX++ = _mm256_or_si256(V2, _mm256_xor_si256(V, tB));											\
	}
#define AVX2_INTR_LCS_INNER_LOOP_SINGLE										\
	{															\
		V = *pX;		\
		tB = _mm256_and_si256(V, _mm256_set_epi64x(*pbm4++, *pbm3++, *pbm2++, *pbm1++));				\
		V2 = _mm256_add_epi64(V, tB);												\
		*pX++ = _mm256_or_si256(V2, _mm256_xor_si256(V, tB));											\
	}
#define AVX2_INTR_LCS_INNER_LOOP										\
	{															\
		V = *pX;		\
		tB = _mm256_and_si256(V, _mm256_set_epi64x(*pbm4++, *pbm3++, *pbm2++, *pbm1++));				\
		V2 = _mm256_sub_epi64(_mm256_add_epi64(V, tB), sB);												\
		sB = _mm256_cmpgt_epi64(_mm256_xor_si256(V, sign64_bit), _mm256_xor_si256(V2, sign64_bit));		\
		*pX++ = _mm256_or_si256(V2, _mm256_xor_si256(V, tB));											\
	}
#define AVX2_INTR_LCS_INNER_LOOP_LAST										\
	{															\
		V = *pX;		\
		tB = _mm256_and_si256(V, _mm256_set_epi64x(*pbm4++, *pbm3++, *pbm2++, *pbm1++));				\
		V2 = _mm256_sub_epi64(_mm256_add_epi64(V, tB), sB);												\
		*pX++ = _mm256_or_si256(V2, _mm256_xor_si256(V, tB));											\
	}
	*/
#define AVX2_INTR_LCS_INNER_LOOP_FIRST										\
	{															\
		V = *pX;		\
		tB = _mm256_and_si256(V, _mm256_set_m128i(*pp34++, *pp12++));				\
		V2 = _mm256_add_epi64(V, tB);												\
		sB = _mm256_cmpgt_epi64(_mm256_xor_si256(V, sign64_bit), _mm256_xor_si256(V2, sign64_bit));		\
		*pX++ = _mm256_or_si256(V2, _mm256_xor_si256(V, tB));											\
	}	
#define AVX2_INTR_LCS_INNER_LOOP_SINGLE										\
	{															\
		V = *pX;		\
		tB = _mm256_and_si256(V, _mm256_set_m128i(*pp34++, *pp12++));				\
		V2 = _mm256_add_epi64(V, tB);												\
		*pX++ = _mm256_or_si256(V2, _mm256_xor_si256(V, tB));											\
	}	
#define AVX2_INTR_LCS_INNER_LOOP										\
	{															\
		V = *pX;		\
		tB = _mm256_and_si256(V, _mm256_set_m128i(*pp34++, *pp12++));				\
		V2 = _mm256_sub_epi64(_mm256_add_epi64(V, tB), sB);												\
		sB = _mm256_cmpgt_epi64(_mm256_xor_si256(V, sign64_bit), _mm256_xor_si256(V2, sign64_bit));		\
		*pX++ = _mm256_or_si256(V2, _mm256_xor_si256(V, tB));											\
	}	
#define AVX2_INTR_LCS_INNER_LOOP_LAST										\
	{															\
		V = *pX;		\
		tB = _mm256_and_si256(V, _mm256_set_m128i(*pp34++, *pp12++));				\
		V2 = _mm256_sub_epi64(_mm256_add_epi64(V, tB), sB);												\
		*pX++ = _mm256_or_si256(V2, _mm256_xor_si256(V, tB));											\
	}	

	// *******************************************************************
	static void LoopCalculate(__m128i* precomp_masks, CSequence* seq0, SeqType* seq1, SeqType* seq2, SeqType* seq3, SeqType* seq4, uint32_t* res, uint32_t bv_len, uint32_t max_len, __m256i* X)
	{
		__m256i V, tB, V2, sB, U2;
		__m256i sign64_bit = _mm256_set1_epi64x(1ull << 63);
		__m256i ones = _mm256_set1_epi64x(~0ull);

		for (size_t i = 0; i < bv_len; ++i)
			X[i] = ones;

		auto pc1 = seq1->data;
		auto pc2 = seq2->data;
		auto pc3 = seq3->data;
		auto pc4 = seq4->data;

		for (size_t i = 0; i < max_len; ++i)
		{
			sB = _mm256_setzero_si256();

			for (size_t j = 0; j < bv_len; ++j)
			{
				V = X[j];
				U2 = _mm256_set_m128i(precomp_masks[(size_t)*pc3 * NO_SYMBOLS * bv_len + *pc4 * bv_len + j], precomp_masks[(size_t)*pc1 * NO_SYMBOLS * bv_len + *pc2 * bv_len + j]);
				tB = _mm256_and_si256(V, U2);
				V2 = _mm256_sub_epi64(_mm256_add_epi64(V, tB), sB);
				sB = _mm256_cmpgt_epi64(_mm256_xor_si256(V, sign64_bit), _mm256_xor_si256(V2, sign64_bit));
				X[j] = _mm256_or_si256(V2, _mm256_sub_epi64(V, tB));
			}
			++pc1; ++pc2; ++pc3; ++pc4;
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
	static void UnrolledCalculate(__m128i* precomp_masks, CSequence* seq0, SeqType* seq1, SeqType* seq2, SeqType* seq3, SeqType* seq4, uint32_t* res, uint32_t max_len, __m256i* X)
	{
		__m256i V, tB, V2, sB;
		__m256i sign64_bit = _mm256_set1_epi64x(1ull << 63);
		__m256i ones = _mm256_set1_epi64x(~0ull);

		uint32_t bv_len = (seq0->length + bv_size256 - 1) / bv_size256;

		auto pX0 = X;

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

		auto pc1 = seq1->data;
		auto pc2 = seq2->data;
		auto pc3 = seq3->data;
		auto pc4 = seq4->data;


		for (size_t i = 0; i < max_len; ++i)
		{
			sB = _mm256_setzero_si256();
			auto pX = X;

			__m128i* pp34 = precomp_masks + *pc3 * NO_SYMBOLS * bv_len + *pc4 * bv_len;
			__m128i* pp12 = precomp_masks + *pc1 * NO_SYMBOLS * bv_len + *pc2 * bv_len;

			if (BV_LEN == 1)		AVX2_INTR_LCS_INNER_LOOP_SINGLE
			else if (BV_LEN > 0)	AVX2_INTR_LCS_INNER_LOOP_FIRST;
			if (BV_LEN == 2)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 1)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 3)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 2)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 4)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 3)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 5)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 4)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 6)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 5)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 7)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 6)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 8)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 7)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 9)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 8)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 10)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 9)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 11)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 10)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 12)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 11)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 13)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 12)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 14)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 13)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 15)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 14)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 16)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 15)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 17)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 16)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 18)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 17)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 19)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 18)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 20)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 19)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 21)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 20)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 22)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 21)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 23)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 22)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 24)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 23)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 25)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 24)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 26)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 25)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 27)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 26)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 28)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 27)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 29)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 28)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 30)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 29)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 31)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 30)	AVX2_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 32)		AVX2_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 31)	AVX2_INTR_LCS_INNER_LOOP;

			++pc1; ++pc2; ++pc3; ++pc4;
		}

		alignas(32) uint64_t p[4];
		auto pX = X;
		if (BV_LEN > 0)		AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 1)		AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 2)		AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 3)		AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 4)		AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 5)		AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 6)		AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 7)		AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 8)		AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 9)		AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 10)	AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 11)	AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 12)	AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 13)	AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 14)	AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 15)	AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 16)	AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 17)	AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 18)	AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 19)	AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 20)	AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 21)	AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 22)	AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 23)	AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 24)	AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 25)	AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 26)	AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 27)	AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 28)	AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 29)	AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 30)	AVX2_INTR_POP_CNT_LOOP;
		if (BV_LEN > 31)	AVX2_INTR_POP_CNT_LOOP;
	}
};

#endif
