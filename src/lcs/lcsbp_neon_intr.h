/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#ifndef _LCSBP_NEON_INTR_H
#define _LCSBP_NEON_INTR_H

#ifdef __ARM_NEON

#include "../core/sequence.h"
#include "../utils/meta_oper.h"
#include <arm_neon.h>

template <unsigned BV_LEN, typename SeqType> class CLCSBP_NEON_INTR_Impl;

// *******************************************************************
//
class CLCSBP_NEON_INTR
{
public:
	void *raw_X;
	void *orig_X;
	uint32_t X_size;
	size_t raw_X_size;
	int64x2_t*X;

	inline void prepare_X(uint32_t bv_len);

public:
	CLCSBP_NEON_INTR() {
		X = nullptr;
		raw_X = nullptr;
		orig_X = nullptr;
		X_size = 0;
		raw_X_size = 0;
	};

	~CLCSBP_NEON_INTR()
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
#define NEON_INTR_POP_CNT_LOOP					\
	{										\
		vst1q_s64((int64_t*) p, *pX);		\
		res[0] += POPCNT(~p[0]);				\
		res[1] += POPCNT(~p[1]);				\
		++pX;								\
	}

	// Important: true in boolean vectors is represented as -1, so here we use "-sB" instead of "+sB" in the classic code
#define NEON_INTR_LCS_INNER_LOOP_FIRST					\
	{										\
		V = *pX;					\
		tB = vandq_s64(V, vcombine_s64(vcreate_s64(*pbm1++), vcreate_s64(*pbm2++)));				\
		V2 = vaddq_s64(V, tB);												\
		sB = (int64x2_t) vcgtq_s64(veorq_s64(V, sign64_bit), veorq_s64(V2, sign64_bit));		\
		*pX++ = vorrq_s64(V2, veorq_s64(V, tB));											\
	}	
#define NEON_INTR_LCS_INNER_LOOP_SINGLE					\
	{										\
		V = *pX;					\
		tB = vandq_s64(V, vcombine_s64(vcreate_s64(*pbm1++), vcreate_s64(*pbm2++)));				\
		V2 = vaddq_s64(V, tB);												\
		*pX++ = vorrq_s64(V2, veorq_s64(V, tB));											\
	}	
#define NEON_INTR_LCS_INNER_LOOP					\
	{										\
		V = *pX;					\
		tB = vandq_s64(V, vcombine_s64(vcreate_s64(*pbm1++), vcreate_s64(*pbm2++)));				\
		V2 = vsubq_s64(vaddq_s64(V, tB), sB);												\
		sB = (int64x2_t) vcgtq_s64(veorq_s64(V, sign64_bit), veorq_s64(V2, sign64_bit));		\
		*pX++ = vorrq_s64(V2, veorq_s64(V, tB));											\
	}	
#define NEON_INTR_LCS_INNER_LOOP_LAST					\
	{										\
		V = *pX;					\
		tB = vandq_s64(V, vcombine_s64(vcreate_s64(*pbm1++), vcreate_s64(*pbm2++)));				\
		V2 = vsubq_s64(vaddq_s64(V, tB), sB);												\
		*pX++ = vorrq_s64(V2, veorq_s64(V, tB));											\
	}	

	// *******************************************************************
	static void LoopCalculate(CSequence* seq0, SeqType* seq1, SeqType* seq2, uint32_t* res, uint32_t bv_len, uint32_t max_len, int64x2_t* X)
	{
		int64x2_t V, tB, V2, sB;
		int64x2_t sign64_bit = vdupq_n_s64(1ull << 63);
		int64x2_t ones = vdupq_n_s64(~0ull);

		bit_vec_t* bit_masks = seq0->p_bit_masks;
		uint32_t loc_bv_len = seq0->p_bv_len;

		auto pc1 = seq1->data;
		auto pc2 = seq2->data;

		for (size_t i = 0; i < loc_bv_len; ++i)
			X[i] = ones;

		for (size_t i = 0; i < max_len; ++i)
		{
			sB = vdupq_n_s64(0);

			for (size_t j = 0; j < loc_bv_len; ++j)
			{
				V = X[j];
				tB = vandq_s64(V, vcombine_s64(vcreate_s64(bit_masks[*pc1 * bv_len + j]), vcreate_s64(bit_masks[*pc2 * bv_len + j])));
				V2 = vsubq_s64(vaddq_s64(V, tB), sB);
				sB = (int64x2_t)vcgtq_s64(veorq_s64(V, sign64_bit), veorq_s64(V2, sign64_bit));
				X[j] = vorrq_s64(V2, vsubq_s64(V, tB));
			}

			++pc1; ++pc2;
		}

		alignas(32) uint64_t p[2];
		for (size_t i = 0; i < loc_bv_len; ++i)
		{
			vst1q_s64((int64_t*)p, X[i]);

			res[0] += POPCNT(~p[0]);
			res[1] += POPCNT(~p[1]);
		}
	}

	// *******************************************************************
	static void UnrolledCalculate(CSequence *seq0, SeqType*seq1, SeqType*seq2, uint32_t *res, uint32_t max_len, int64x2_t *X)
	{
		int64x2_t V, tB, V2, sB;
		int64x2_t sign64_bit = vdupq_n_s64(1ull << 63);
		int64x2_t ones = vdupq_n_s64(~0ull);

		auto pX0 = X;
		bit_vec_t* bit_masks = seq0->p_bit_masks;

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
			sB = vdupq_n_s64(0);

			auto pX = X;

			uint64_t* pbm1 = bm + *pc1 * bm_len;
			uint64_t* pbm2 = bm + *pc2 * bm_len;

			if (BV_LEN == 1)		NEON_INTR_LCS_INNER_LOOP_SINGLE
			else if (BV_LEN > 0)	NEON_INTR_LCS_INNER_LOOP_FIRST;
			if (BV_LEN == 2)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 1)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 3)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 2)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 4)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 3)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 5)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 4)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 6)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 5)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 7)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 6)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 8)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 7)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 9)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 8)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 10)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 9)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 11)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 10)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 12)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 11)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 13)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 12)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 14)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 13)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 15)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 14)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 16)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 15)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 17)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 16)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 18)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 17)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 19)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 18)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 20)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 19)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 21)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 20)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 22)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 21)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 23)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 22)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 24)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 23)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 25)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 24)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 26)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 25)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 27)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 26)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 28)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 27)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 29)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 28)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 30)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 29)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 31)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 30)	NEON_INTR_LCS_INNER_LOOP;
			if (BV_LEN == 32)		NEON_INTR_LCS_INNER_LOOP_LAST
			else if (BV_LEN > 31)	NEON_INTR_LCS_INNER_LOOP;

			++pc1; ++pc2;
		}

		alignas(32) uint64_t p[2];
		auto pX = X;
		if (BV_LEN > 0)		NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 1)		NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 2)		NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 3)		NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 4)		NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 5)		NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 6)		NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 7)		NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 8)		NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 9)		NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 10)	NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 11)	NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 12)	NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 13)	NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 14)	NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 15)	NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 16)	NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 17)	NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 18)	NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 19)	NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 20)	NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 21)	NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 22)	NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 23)	NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 24)	NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 25)	NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 26)	NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 27)	NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 28)	NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 29)	NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 30)	NEON_INTR_POP_CNT_LOOP;
		if (BV_LEN > 31)	NEON_INTR_POP_CNT_LOOP;
	}
};

#endif
#endif
