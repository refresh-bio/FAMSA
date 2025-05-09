/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#ifndef _LCSBP_AVX512_INTR_H
#define _LCSBP_AVX512_INTR_H

#include "../core/sequence.h"
#include "../utils/meta_oper.h"
#include "../utils/utils.h"
#include "../core/defs.h"

#include <algorithm>
#include <memory>
#include <immintrin.h>

template <unsigned BV_LEN, typename SeqType> class CLCSBP_AVX512_INTR_Impl;

// *******************************************************************
//
class CLCSBP_AVX512_INTR
{
	void* raw_X;
	void* orig_X;
	uint32_t X_size;
	size_t raw_X_size;
	__m512i* X;
	int prev_sequence_no;

	__m128i* precomp_masks;
	void* raw_precomp_masks;
	size_t size_precomp_masks;

	inline void prepare_X(uint32_t bv_len);
	inline void prepare_mask_pairs(uint32_t bv_len, CSequence* seq0);

public:
	CLCSBP_AVX512_INTR() {
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

	~CLCSBP_AVX512_INTR()
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

	void Calculate(CSequence* seq0, CSequence* seq1, CSequence* seq2, CSequence* seq3, CSequence* seq4, CSequence* seq5, CSequence* seq6, CSequence* seq7, CSequence* seq8,
		uint32_t* dist);
	void Calculate(CSequence* seq0, CSequenceView* sv1, CSequenceView* sv2, CSequenceView* sv3, CSequenceView* sv4, CSequenceView* sv5, CSequenceView* sv6, CSequenceView* sv7, CSequenceView* sv8,
		uint32_t* dist);
};

// *******************************************************************
//
template <unsigned BV_LEN, typename SeqType> class CLCSBP_AVX512_INTR_Impl {
public:
#define AVX512_INTR_POP_CNT_LOOP						\
	{												\
		_mm512_storeu_si512((__m512i*) p, *pX);		\
		res[0] += POPCNT(~p[0]);				\
		res[1] += POPCNT(~p[1]);				\
		res[2] += POPCNT(~p[2]);				\
		res[3] += POPCNT(~p[3]);				\
		res[4] += POPCNT(~p[4]);				\
		res[5] += POPCNT(~p[5]);				\
		res[6] += POPCNT(~p[6]);				\
		res[7] += POPCNT(~p[7]);				\
		++pX;										\
	}

#define AVX512_LOADING_MODE 0

#if AVX512_LOADING_MODE==0
#define AVX512_CONSTRUCT_U2(pa, pb, pc, pd)		\
	    m4A = _mm256_set_m128i(*(pb), *(pa));	\
		m4B = _mm256_set_m128i(*(pd), *(pc));	\
		U2 = _mm512_inserti64x4(_mm512_zextsi256_si512(m4A), m4B, 1);
#elif AVX512_LOADING_MODE==1
#define AVX512_CONSTRUCT_U2(pa, pb, pc, pd)		\
	U2 = _mm512_castsi128_si512(*(pa));		\
	U2 = _mm512_inserti64x2(U2, *(pb), 1);	\
	U2 = _mm512_inserti64x2(U2, *(pc), 2);	\
	U2 = _mm512_inserti64x2(U2, *(pd), 3);
#elif AVX512_LOADING_MODE==2
#define AVX512_CONSTRUCT_U2(pa, pb, pc, pd)		\
	m2A = *(pa);								\
	m2B = *(pc);								\
	m4A = _mm256_inserti128_si256(_mm256_zextsi128_si256(m2A), *(pb), 1);	\
	m4B = _mm256_inserti128_si256(_mm256_zextsi128_si256(m2B), *(pd), 1);	\
	U2 = _mm512_inserti64x4(_mm512_zextsi256_si512(m4A), m4B, 1);
#endif

#define AVX512_INTR_LCS_INNER_LOOP_FIRST										\
	{															\
		V = *X;		\
		AVX512_CONSTRUCT_U2(pp12, pp34, pp56, pp78)	\
		tB = _mm512_and_si512(V, U2);				\
		V2 = _mm512_add_epi64(V, tB);												\
		cmp_mask = _mm512_cmpgt_epu64_mask(V, V2);\
		sB = _mm512_movm_epi64 (cmp_mask);		\
		*X = _mm512_or_si512(V2, _mm512_xor_si512(V, tB));											\
	}	

#define AVX512_INTR_LCS_INNER_LOOP_SINGLE										\
	{															\
		V = *X;		\
		AVX512_CONSTRUCT_U2(pp12, pp34, pp56, pp78)	\
		tB = _mm512_and_si512(V, U2);				\
		V2 = _mm512_add_epi64(V, tB);												\
		*X = _mm512_or_si512(V2, _mm512_xor_si512(V, tB));											\
	}

#define AVX512_INTR_LCS_INNER_LOOP(n)										\
	{															\
		V = X[n];		\
		AVX512_CONSTRUCT_U2(pp12+n, pp34+n, pp56+n, pp78+n)	\
		tB = _mm512_and_si512(V, U2);				\
		V2 = _mm512_sub_epi64(_mm512_add_epi64(V, tB), sB);												\
		cmp_mask = _mm512_cmpgt_epu64_mask(V, V2);\
		sB = _mm512_movm_epi64 (cmp_mask);		\
		X[n] = _mm512_or_si512(V2, _mm512_xor_si512(V, tB));											\
	}	

#define AVX512_INTR_LCS_INNER_LOOP_LAST(n)										\
	{															\
		V = X[n];		\
		AVX512_CONSTRUCT_U2(pp12+n, pp34+n, pp56+n, pp78+n)	\
		tB = _mm512_and_si512(V, U2);				\
		V2 = _mm512_sub_epi64(_mm512_add_epi64(V, tB), sB);					\
		X[n] = _mm512_or_si512(V2, _mm512_xor_si512(V, tB));											\
	}

	// *******************************************************************
	static void LoopCalculate(__m128i* precomp_masks, CSequence* seq0, 
		SeqType* seq1, SeqType* seq2, SeqType* seq3, SeqType* seq4, SeqType* seq5, SeqType* seq6, SeqType* seq7, SeqType* seq8,
		uint32_t* res, uint32_t bv_len, uint32_t max_len, __m512i* X)
	{
		__m512i V, tB, V2, sB, U2;
//		__mmask8 cmp_mask;
		const __m512i sign64_bit = _mm512_set1_epi64(1ull << 63);
		const __m512i ones = _mm512_set1_epi64(~0ull);

		for (size_t i = 0; i < bv_len; ++i)
			X[i] = ones;

		auto pc1 = seq1->data;
		auto pc2 = seq2->data;
		auto pc3 = seq3->data;
		auto pc4 = seq4->data;
		auto pc5 = seq5->data;
		auto pc6 = seq6->data;
		auto pc7 = seq7->data;
		auto pc8 = seq8->data;

#if 1
		for (size_t i = 0; i < max_len; ++i)
		{
			sB = _mm512_setzero_si512();

			__m128i* pp12 = precomp_masks + (*pc1 * NO_AMINOACIDS + *pc2) * bv_len;
			__m128i* pp34 = precomp_masks + (*pc3 * NO_AMINOACIDS + *pc4) * bv_len;
			__m128i* pp56 = precomp_masks + (*pc5 * NO_AMINOACIDS + *pc6) * bv_len;
			__m128i* pp78 = precomp_masks + (*pc7 * NO_AMINOACIDS + *pc8) * bv_len;

			for (size_t j = 0; j < bv_len; ++j)
			{
				V = X[j];
				U2 = _mm512_castsi128_si512(*pp12++);
				U2 = _mm512_inserti64x2(U2, *pp34++, 1);
				U2 = _mm512_inserti64x2(U2, *pp56++, 2);
				U2 = _mm512_inserti64x2(U2, *pp78++, 3);

				tB = _mm512_and_si512(V, U2);
				V2 = _mm512_add_epi64(_mm512_add_epi64(V, tB), sB);

				__mmask8 cmp_mask = _mm512_cmpgt_epi64_mask(_mm512_xor_si512(V, sign64_bit), _mm512_xor_si512(V2, sign64_bit));
				sB = _mm512_maskz_mov_epi64(cmp_mask, _mm512_set1_epi64(1));

				X[j] = _mm512_or_si512(V2, _mm512_xor_si512(V, tB));
			}

			++pc1; ++pc2; ++pc3; ++pc4; ++pc5; ++pc6; ++pc7; ++pc8;
		}
#else
		for (size_t i = 0; i < max_len; ++i)
		{
			sB = _mm256_setzero_si256();
			size_t j;

			__m128i* pp12 = precomp_masks + (*pc1 * NO_AMINOACIDS + *pc2) * bv_len;
			__m128i* pp34 = precomp_masks + (*pc3 * NO_AMINOACIDS + *pc4) * bv_len;

			for (j = 0; j + 1 < bv_len; ++j)
			{
				V = X[j];
				U2 = _mm256_set_m128i(*pp34++, *pp12++);
				tB = _mm256_and_si256(V, U2);
				V2 = _mm256_sub_epi64(_mm256_add_epi64(V, tB), sB);
				sB = _mm256_cmpgt_epi64(_mm256_xor_si256(V, sign64_bit), _mm256_xor_si256(V2, sign64_bit));
				X[j] = _mm256_or_si256(V2, _mm256_sub_epi64(V, tB));

				++j;

				V = X[j];
				U2 = _mm256_set_m128i(*pp34++, *pp12++);
				tB = _mm256_and_si256(V, U2);
				V2 = _mm256_sub_epi64(_mm256_add_epi64(V, tB), sB);
				sB = _mm256_cmpgt_epi64(_mm256_xor_si256(V, sign64_bit), _mm256_xor_si256(V2, sign64_bit));
				X[j] = _mm256_or_si256(V2, _mm256_sub_epi64(V, tB));
			}

			if (j < bv_len)
			{
				V = X[j];
				U2 = _mm256_set_m128i(*pp34, *pp12);
				tB = _mm256_and_si256(V, U2);
				V2 = _mm256_sub_epi64(_mm256_add_epi64(V, tB), sB);
				sB = _mm256_cmpgt_epi64(_mm256_xor_si256(V, sign64_bit), _mm256_xor_si256(V2, sign64_bit));
				X[j] = _mm256_or_si256(V2, _mm256_sub_epi64(V, tB));
			}

			++pc1; ++pc2; ++pc3; ++pc4;
		}
#endif

		alignas(64) uint64_t p[8];
#ifdef _MSC_VER					// Visual C++

		for (size_t i = 0; i < bv_len; ++i)
		{
			_mm512_storeu_si512((__m512i*) p, X[i]);

			res[0] += POPCNT(~p[0]);
			res[1] += POPCNT(~p[1]);
			res[2] += POPCNT(~p[2]);
			res[3] += POPCNT(~p[3]);
			res[4] += POPCNT(~p[4]);
			res[5] += POPCNT(~p[5]);
			res[6] += POPCNT(~p[6]);
			res[7] += POPCNT(~p[7]);
		}
#else
#ifdef __GNUC__
		for (size_t i = 0; i < bv_len; ++i)
		{
			_mm512_storeu_si512((__m512i*) p, X[i]);

			res[0] += POPCNT(~p[0]);
			res[1] += POPCNT(~p[1]);
			res[2] += POPCNT(~p[2]);
			res[3] += POPCNT(~p[3]);
			res[4] += POPCNT(~p[4]);
			res[5] += POPCNT(~p[5]);
			res[6] += POPCNT(~p[6]);
			res[7] += POPCNT(~p[7]);
		}
#else
		for (size_t i = 0; i < bv_len; ++i)
		{
			_mm512_storeu_si512((__m512i*) p, X[i]);
			for (int v = 0; v < 8; ++v)
				for (uint64_t T = ~p[v]; T; T &= T - 1)
					++res[v];
		}
#endif
#endif

	}

	// *******************************************************************
	static void UnrolledCalculate(__m128i* precomp_masks, CSequence* seq0, 
		SeqType* seq1, SeqType* seq2, SeqType* seq3, SeqType* seq4, SeqType* seq5, SeqType* seq6, SeqType* seq7, SeqType* seq8,
		uint32_t* res, uint32_t max_len, __m512i* X1)
	{
		__m512i alignas(64) V, tB, V2, sB, U2;
		__mmask8 alignas(16) cmp_mask;
		const __m512i alignas(64) sign64_bit = _mm512_set1_epi64(1ull << 63);
		const __m512i alignas(64) ones = _mm512_set1_epi64(~0ull);

		const uint32_t bv_len = (seq0->length + bv_size512 - 1) / bv_size512;

		__m128i alignas(16) m2A, m2B;
		__m256i alignas(32) m4A, m4B;

		__m512i alignas(64) X[BV_LEN];

		if constexpr (BV_LEN > 0)			X[0] = ones;
		if constexpr (BV_LEN > 1)			X[1] = ones;
		if constexpr (BV_LEN > 2)			X[2] = ones;
		if constexpr (BV_LEN > 3)			X[3] = ones;
		if constexpr (BV_LEN > 4)			X[4] = ones;
		if constexpr (BV_LEN > 5)			X[5] = ones;
		if constexpr (BV_LEN > 6)			X[6] = ones;
		if constexpr (BV_LEN > 7)			X[7] = ones;
		if constexpr (BV_LEN > 8)			X[8] = ones;
		if constexpr (BV_LEN > 9)			X[9] = ones;
		if constexpr (BV_LEN > 10)			X[10] = ones;
		if constexpr (BV_LEN > 11)			X[11] = ones;
		if constexpr (BV_LEN > 12)			X[12] = ones;
		if constexpr (BV_LEN > 13)			X[13] = ones;
		if constexpr (BV_LEN > 14)			X[14] = ones;
		if constexpr (BV_LEN > 15)			X[15] = ones;
		if constexpr (BV_LEN > 16)			X[16] = ones;
		if constexpr (BV_LEN > 17)			X[17] = ones;
		if constexpr (BV_LEN > 18)			X[18] = ones;
		if constexpr (BV_LEN > 19)			X[19] = ones;
		if constexpr (BV_LEN > 20)			X[20] = ones;
		if constexpr (BV_LEN > 21)			X[21] = ones;
		if constexpr (BV_LEN > 22)			X[22] = ones;
		if constexpr (BV_LEN > 23)			X[23] = ones;
		if constexpr (BV_LEN > 24)			X[24] = ones;
		if constexpr (BV_LEN > 25)			X[25] = ones;
		if constexpr (BV_LEN > 26)			X[26] = ones;
		if constexpr (BV_LEN > 27)			X[27] = ones;
		if constexpr (BV_LEN > 28)			X[28] = ones;
		if constexpr (BV_LEN > 29)			X[29] = ones;
		if constexpr (BV_LEN > 30)			X[30] = ones;
		if constexpr (BV_LEN > 31)			X[31] = ones;

		auto pc1 = seq1->data;
		auto pc2 = seq2->data;
		auto pc3 = seq3->data;
		auto pc4 = seq4->data;
		auto pc5 = seq5->data;
		auto pc6 = seq6->data;
		auto pc7 = seq7->data;
		auto pc8 = seq8->data;

		size_t i;

#if 0
		// For some reason, this code is not faster than the one below
		const __m256i v_bv_len = _mm256_set1_epi32(bv_len);
		const __m256i v_NO_AA = _mm256_set1_epi32(NO_AMINOACIDS);
		int alignas(32) idx12[8], idx34[8], idx56[8], idx78[8];

		for (i = 0; i + 8 < max_len; i += 8)
		{
			__m128i vec_p1 = _mm_loadl_epi64((__m128i const*)(pc1));
			__m128i vec_p2 = _mm_loadl_epi64((__m128i const*)(pc2));
			__m128i vec_p3 = _mm_loadl_epi64((__m128i const*)(pc3));
			__m128i vec_p4 = _mm_loadl_epi64((__m128i const*)(pc4));
			__m128i vec_p5 = _mm_loadl_epi64((__m128i const*)(pc5));
			__m128i vec_p6 = _mm_loadl_epi64((__m128i const*)(pc6));
			__m128i vec_p7 = _mm_loadl_epi64((__m128i const*)(pc7));
			__m128i vec_p8 = _mm_loadl_epi64((__m128i const*)(pc8));

			__m256i v_p1 = _mm256_cvtepu8_epi32(vec_p1);
			__m256i v_p2 = _mm256_cvtepu8_epi32(vec_p2);
			__m256i v_p3 = _mm256_cvtepu8_epi32(vec_p3);
			__m256i v_p4 = _mm256_cvtepu8_epi32(vec_p4);
			__m256i v_p5 = _mm256_cvtepu8_epi32(vec_p5);
			__m256i v_p6 = _mm256_cvtepu8_epi32(vec_p6);
			__m256i v_p7 = _mm256_cvtepu8_epi32(vec_p7);
			__m256i v_p8 = _mm256_cvtepu8_epi32(vec_p8);

			__m256i v_idx12 = _mm256_mullo_epi32(_mm256_add_epi32(_mm256_mullo_epi32(v_p1, v_NO_AA), v_p2), v_bv_len);
			__m256i v_idx34 = _mm256_mullo_epi32(_mm256_add_epi32(_mm256_mullo_epi32(v_p3, v_NO_AA), v_p4), v_bv_len);
			__m256i v_idx56 = _mm256_mullo_epi32(_mm256_add_epi32(_mm256_mullo_epi32(v_p5, v_NO_AA), v_p6), v_bv_len);
			__m256i v_idx78 = _mm256_mullo_epi32(_mm256_add_epi32(_mm256_mullo_epi32(v_p7, v_NO_AA), v_p8), v_bv_len);

			_mm256_storeu_si256((__m256i*)idx12, v_idx12);
			_mm256_storeu_si256((__m256i*)idx34, v_idx34);
			_mm256_storeu_si256((__m256i*)idx56, v_idx56);
			_mm256_storeu_si256((__m256i*)idx78, v_idx78);

			for (int j = 0; j < 8; ++j)
			{
				__m128i* pp12 = precomp_masks + idx12[j];
				__m128i* pp34 = precomp_masks + idx34[j];
				__m128i* pp56 = precomp_masks + idx56[j];
				__m128i* pp78 = precomp_masks + idx78[j];

				if constexpr (BV_LEN == 1)		AVX512_INTR_LCS_INNER_LOOP_SINGLE
				else if constexpr (BV_LEN > 0)	AVX512_INTR_LCS_INNER_LOOP_FIRST;
				if constexpr (BV_LEN == 2)		AVX512_INTR_LCS_INNER_LOOP_LAST(1)
				else if constexpr (BV_LEN > 1)	AVX512_INTR_LCS_INNER_LOOP(1);
				if constexpr (BV_LEN == 3)		AVX512_INTR_LCS_INNER_LOOP_LAST(2)
				else if constexpr (BV_LEN > 2)	AVX512_INTR_LCS_INNER_LOOP(2);
				if constexpr (BV_LEN == 4)		AVX512_INTR_LCS_INNER_LOOP_LAST(3)
				else if constexpr (BV_LEN > 3)	AVX512_INTR_LCS_INNER_LOOP(3);
				if constexpr (BV_LEN == 5)		AVX512_INTR_LCS_INNER_LOOP_LAST(4)
				else if constexpr (BV_LEN > 4)	AVX512_INTR_LCS_INNER_LOOP(4);
				if constexpr (BV_LEN == 6)		AVX512_INTR_LCS_INNER_LOOP_LAST(5)
				else if constexpr (BV_LEN > 5)	AVX512_INTR_LCS_INNER_LOOP(5);
				if constexpr (BV_LEN == 7)		AVX512_INTR_LCS_INNER_LOOP_LAST(6)
				else if constexpr (BV_LEN > 6)	AVX512_INTR_LCS_INNER_LOOP(6);
				if constexpr (BV_LEN == 8)		AVX512_INTR_LCS_INNER_LOOP_LAST(7)
				else if constexpr (BV_LEN > 7)	AVX512_INTR_LCS_INNER_LOOP(7);
				if constexpr (BV_LEN == 9)		AVX512_INTR_LCS_INNER_LOOP_LAST(8)
				else if constexpr (BV_LEN > 8)	AVX512_INTR_LCS_INNER_LOOP(8);
				if constexpr (BV_LEN == 10)		AVX512_INTR_LCS_INNER_LOOP_LAST(9)
				else if constexpr (BV_LEN > 9)	AVX512_INTR_LCS_INNER_LOOP(9);
				if constexpr (BV_LEN == 11)		AVX512_INTR_LCS_INNER_LOOP_LAST(10)
				else if constexpr (BV_LEN > 10)	AVX512_INTR_LCS_INNER_LOOP(10);
				if constexpr (BV_LEN == 12)		AVX512_INTR_LCS_INNER_LOOP_LAST(11)
				else if constexpr (BV_LEN > 11)	AVX512_INTR_LCS_INNER_LOOP(11);
				if constexpr (BV_LEN == 13)		AVX512_INTR_LCS_INNER_LOOP_LAST(12)
				else if constexpr (BV_LEN > 12)	AVX512_INTR_LCS_INNER_LOOP(12);
				if constexpr (BV_LEN == 14)		AVX512_INTR_LCS_INNER_LOOP_LAST(13)
				else if constexpr (BV_LEN > 13)	AVX512_INTR_LCS_INNER_LOOP(13);
				if constexpr (BV_LEN == 15)		AVX512_INTR_LCS_INNER_LOOP_LAST(14)
				else if constexpr (BV_LEN > 14)	AVX512_INTR_LCS_INNER_LOOP(14);
				if constexpr (BV_LEN == 16)		AVX512_INTR_LCS_INNER_LOOP_LAST(15)
				else if constexpr (BV_LEN > 15)	AVX512_INTR_LCS_INNER_LOOP(15);
				if constexpr (BV_LEN == 17)		AVX512_INTR_LCS_INNER_LOOP_LAST(16)
				else if constexpr (BV_LEN > 16)	AVX512_INTR_LCS_INNER_LOOP(16);
				if constexpr (BV_LEN == 18)		AVX512_INTR_LCS_INNER_LOOP_LAST(17)
				else if constexpr (BV_LEN > 17)	AVX512_INTR_LCS_INNER_LOOP(17);
				if constexpr (BV_LEN == 19)		AVX512_INTR_LCS_INNER_LOOP_LAST(18)
				else if constexpr (BV_LEN > 18)	AVX512_INTR_LCS_INNER_LOOP(18);
				if constexpr (BV_LEN == 20)		AVX512_INTR_LCS_INNER_LOOP_LAST(19)
				else if constexpr (BV_LEN > 19)	AVX512_INTR_LCS_INNER_LOOP(19);
				if constexpr (BV_LEN == 21)		AVX512_INTR_LCS_INNER_LOOP_LAST(20)
				else if constexpr (BV_LEN > 20)	AVX512_INTR_LCS_INNER_LOOP(20);
				if constexpr (BV_LEN == 22)		AVX512_INTR_LCS_INNER_LOOP_LAST(21)
				else if constexpr (BV_LEN > 21)	AVX512_INTR_LCS_INNER_LOOP(21);
				if constexpr (BV_LEN == 23)		AVX512_INTR_LCS_INNER_LOOP_LAST(22)
				else if constexpr (BV_LEN > 22)	AVX512_INTR_LCS_INNER_LOOP(22);
				if constexpr (BV_LEN == 24)		AVX512_INTR_LCS_INNER_LOOP_LAST(23)
				else if constexpr (BV_LEN > 23)	AVX512_INTR_LCS_INNER_LOOP(23);
				if constexpr (BV_LEN == 25)		AVX512_INTR_LCS_INNER_LOOP_LAST(24)
				else if constexpr (BV_LEN > 24)	AVX512_INTR_LCS_INNER_LOOP(24);
				if constexpr (BV_LEN == 26)		AVX512_INTR_LCS_INNER_LOOP_LAST(25)
				else if constexpr (BV_LEN > 25)	AVX512_INTR_LCS_INNER_LOOP(25);
				if constexpr (BV_LEN == 27)		AVX512_INTR_LCS_INNER_LOOP_LAST(26)
				else if constexpr (BV_LEN > 26)	AVX512_INTR_LCS_INNER_LOOP(26);
				if constexpr (BV_LEN == 28)		AVX512_INTR_LCS_INNER_LOOP_LAST(27)
				else if constexpr (BV_LEN > 27)	AVX512_INTR_LCS_INNER_LOOP(27);
				if constexpr (BV_LEN == 29)		AVX512_INTR_LCS_INNER_LOOP_LAST(28)
				else if constexpr (BV_LEN > 28)	AVX512_INTR_LCS_INNER_LOOP(28);
				if constexpr (BV_LEN == 30)		AVX512_INTR_LCS_INNER_LOOP_LAST(29)
				else if constexpr (BV_LEN > 29)	AVX512_INTR_LCS_INNER_LOOP(29);
				if constexpr (BV_LEN == 31)		AVX512_INTR_LCS_INNER_LOOP_LAST(30)
				else if constexpr (BV_LEN > 30)	AVX512_INTR_LCS_INNER_LOOP(30);
				if constexpr (BV_LEN == 32)		AVX512_INTR_LCS_INNER_LOOP_LAST(31)
				else if constexpr (BV_LEN > 31)	AVX512_INTR_LCS_INNER_LOOP(31);
			}

			pc1 += 8; pc2 += 8; pc3 += 8; pc4 += 8; pc5 += 8; pc6 += 8; pc7 += 8; pc8 += 8;
		}
#endif

		for (; i < max_len; ++i)
		{
			__m128i* pp12 = precomp_masks + (*pc1 * NO_AMINOACIDS + *pc2) * bv_len;
			__m128i* pp34 = precomp_masks + (*pc3 * NO_AMINOACIDS + *pc4) * bv_len;
			__m128i* pp56 = precomp_masks + (*pc5 * NO_AMINOACIDS + *pc6) * bv_len;
			__m128i* pp78 = precomp_masks + (*pc7 * NO_AMINOACIDS + *pc8) * bv_len;

			if constexpr (BV_LEN == 1)		AVX512_INTR_LCS_INNER_LOOP_SINGLE
			else if constexpr (BV_LEN > 0)	AVX512_INTR_LCS_INNER_LOOP_FIRST;
			if constexpr (BV_LEN == 2)		AVX512_INTR_LCS_INNER_LOOP_LAST(1)
			else if constexpr (BV_LEN > 1)	AVX512_INTR_LCS_INNER_LOOP(1);
			if constexpr (BV_LEN == 3)		AVX512_INTR_LCS_INNER_LOOP_LAST(2)
			else if constexpr (BV_LEN > 2)	AVX512_INTR_LCS_INNER_LOOP(2);
			if constexpr (BV_LEN == 4)		AVX512_INTR_LCS_INNER_LOOP_LAST(3)
			else if constexpr (BV_LEN > 3)	AVX512_INTR_LCS_INNER_LOOP(3);
			if constexpr (BV_LEN == 5)		AVX512_INTR_LCS_INNER_LOOP_LAST(4)
			else if constexpr (BV_LEN > 4)	AVX512_INTR_LCS_INNER_LOOP(4);
			if constexpr (BV_LEN == 6)		AVX512_INTR_LCS_INNER_LOOP_LAST(5)
			else if constexpr (BV_LEN > 5)	AVX512_INTR_LCS_INNER_LOOP(5);
			if constexpr (BV_LEN == 7)		AVX512_INTR_LCS_INNER_LOOP_LAST(6)
			else if constexpr (BV_LEN > 6)	AVX512_INTR_LCS_INNER_LOOP(6);
			if constexpr (BV_LEN == 8)		AVX512_INTR_LCS_INNER_LOOP_LAST(7)
			else if constexpr (BV_LEN > 7)	AVX512_INTR_LCS_INNER_LOOP(7);
			if constexpr (BV_LEN == 9)		AVX512_INTR_LCS_INNER_LOOP_LAST(8)
			else if constexpr (BV_LEN > 8)	AVX512_INTR_LCS_INNER_LOOP(8);
			if constexpr (BV_LEN == 10)		AVX512_INTR_LCS_INNER_LOOP_LAST(9)
			else if constexpr (BV_LEN > 9)	AVX512_INTR_LCS_INNER_LOOP(9);
			if constexpr (BV_LEN == 11)		AVX512_INTR_LCS_INNER_LOOP_LAST(10)
			else if constexpr (BV_LEN > 10)	AVX512_INTR_LCS_INNER_LOOP(10);
			if constexpr (BV_LEN == 12)		AVX512_INTR_LCS_INNER_LOOP_LAST(11)
			else if constexpr (BV_LEN > 11)	AVX512_INTR_LCS_INNER_LOOP(11);
			if constexpr (BV_LEN == 13)		AVX512_INTR_LCS_INNER_LOOP_LAST(12)
			else if constexpr (BV_LEN > 12)	AVX512_INTR_LCS_INNER_LOOP(12);
			if constexpr (BV_LEN == 14)		AVX512_INTR_LCS_INNER_LOOP_LAST(13)
			else if constexpr (BV_LEN > 13)	AVX512_INTR_LCS_INNER_LOOP(13);
			if constexpr (BV_LEN == 15)		AVX512_INTR_LCS_INNER_LOOP_LAST(14)
			else if constexpr (BV_LEN > 14)	AVX512_INTR_LCS_INNER_LOOP(14);
			if constexpr (BV_LEN == 16)		AVX512_INTR_LCS_INNER_LOOP_LAST(15)
			else if constexpr (BV_LEN > 15)	AVX512_INTR_LCS_INNER_LOOP(15);
			if constexpr (BV_LEN == 17)		AVX512_INTR_LCS_INNER_LOOP_LAST(16)
			else if constexpr (BV_LEN > 16)	AVX512_INTR_LCS_INNER_LOOP(16);
			if constexpr (BV_LEN == 18)		AVX512_INTR_LCS_INNER_LOOP_LAST(17)
			else if constexpr (BV_LEN > 17)	AVX512_INTR_LCS_INNER_LOOP(17);
			if constexpr (BV_LEN == 19)		AVX512_INTR_LCS_INNER_LOOP_LAST(18)
			else if constexpr (BV_LEN > 18)	AVX512_INTR_LCS_INNER_LOOP(18);
			if constexpr (BV_LEN == 20)		AVX512_INTR_LCS_INNER_LOOP_LAST(19)
			else if constexpr (BV_LEN > 19)	AVX512_INTR_LCS_INNER_LOOP(19);
			if constexpr (BV_LEN == 21)		AVX512_INTR_LCS_INNER_LOOP_LAST(20)
			else if constexpr (BV_LEN > 20)	AVX512_INTR_LCS_INNER_LOOP(20);
			if constexpr (BV_LEN == 22)		AVX512_INTR_LCS_INNER_LOOP_LAST(21)
			else if constexpr (BV_LEN > 21)	AVX512_INTR_LCS_INNER_LOOP(21);
			if constexpr (BV_LEN == 23)		AVX512_INTR_LCS_INNER_LOOP_LAST(22)
			else if constexpr (BV_LEN > 22)	AVX512_INTR_LCS_INNER_LOOP(22);
			if constexpr (BV_LEN == 24)		AVX512_INTR_LCS_INNER_LOOP_LAST(23)
			else if constexpr (BV_LEN > 23)	AVX512_INTR_LCS_INNER_LOOP(23);
			if constexpr (BV_LEN == 25)		AVX512_INTR_LCS_INNER_LOOP_LAST(24)
			else if constexpr (BV_LEN > 24)	AVX512_INTR_LCS_INNER_LOOP(24);
			if constexpr (BV_LEN == 26)		AVX512_INTR_LCS_INNER_LOOP_LAST(25)
			else if constexpr (BV_LEN > 25)	AVX512_INTR_LCS_INNER_LOOP(25);
			if constexpr (BV_LEN == 27)		AVX512_INTR_LCS_INNER_LOOP_LAST(26)
			else if constexpr (BV_LEN > 26)	AVX512_INTR_LCS_INNER_LOOP(26);
			if constexpr (BV_LEN == 28)		AVX512_INTR_LCS_INNER_LOOP_LAST(27)
			else if constexpr (BV_LEN > 27)	AVX512_INTR_LCS_INNER_LOOP(27);
			if constexpr (BV_LEN == 29)		AVX512_INTR_LCS_INNER_LOOP_LAST(28)
			else if constexpr (BV_LEN > 28)	AVX512_INTR_LCS_INNER_LOOP(28);
			if constexpr (BV_LEN == 30)		AVX512_INTR_LCS_INNER_LOOP_LAST(29)
			else if constexpr (BV_LEN > 29)	AVX512_INTR_LCS_INNER_LOOP(29);
			if constexpr (BV_LEN == 31)		AVX512_INTR_LCS_INNER_LOOP_LAST(30)
			else if constexpr (BV_LEN > 30)	AVX512_INTR_LCS_INNER_LOOP(30);
			if constexpr (BV_LEN == 32)		AVX512_INTR_LCS_INNER_LOOP_LAST(31)
			else if constexpr (BV_LEN > 31)	AVX512_INTR_LCS_INNER_LOOP(31);

			++pc1; ++pc2; ++pc3; ++pc4; ++pc5; ++pc6; ++pc7; ++pc8;
		}

		alignas(64) uint64_t p[8];
		auto pX = X;
		if constexpr (BV_LEN > 0)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 1)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 2)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 3)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 4)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 5)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 6)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 7)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 8)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 9)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 10)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 11)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 12)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 13)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 14)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 15)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 16)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 17)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 18)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 19)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 20)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 21)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 22)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 23)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 24)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 25)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 26)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 27)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 28)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 29)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 30)	AVX512_INTR_POP_CNT_LOOP;
		if constexpr (BV_LEN > 31)	AVX512_INTR_POP_CNT_LOOP;
	}
};

#endif
