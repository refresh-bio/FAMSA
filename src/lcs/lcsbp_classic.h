/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#ifndef _LCSBP_CLASSIC_H
#define _LCSBP_CLASSIC_H

#include "../core/sequence.h"
#include "../utils/meta_oper.h"

template <unsigned BV_LEN, typename SeqType> class CLCSBP_Classic_Impl;

using bit_vec_iterator_t = bit_vec_t*;

class CLCSBP_Classic
{
	bit_vec_t *X;
	uint32_t X_size;

	CSequence *pf_seq0;
	bit_vec_iterator_t s0bm[NO_SYMBOLS];
	
	inline void prepare_X(uint32_t bv_len);
	void prefetch_bitmasks(CSequence* seq0);

public:
	CLCSBP_Classic() {
		X_size = 0;
		X = nullptr;

		pf_seq0 = nullptr;
	};
	
	~CLCSBP_Classic() {
		if (X)
			delete[] X;
	};
		
	void Calculate(CSequence* seq0, CSequence* seq1,
		uint32_t* dist);
	void Calculate(CSequence* seq0, CSequenceView* seq1,
		uint32_t* dist);
};

template <unsigned BV_LEN, typename SeqType> class CLCSBP_Classic_Impl {
public:
#define CLASSIC_LCS_INNER_LOOP		\
	{								\
		V = *pX;					\
		tB = V & *s0b++;			\
		V2 = V + tB + sB;			\
		sB = V2 < V;				\
		*pX++ = V2 | (V - tB);		\
	}

#define CLASSIC_POP_CNT_LOOP					\
	{											\
		for (V = ~*pX; V; V &= V - 1)			\
			++res[0];							\
		pX++;									\
	}

	static void LoopCalculate(CSequence* seq0, SeqType* seq1, uint32_t* res, uint32_t bv_len, bit_vec_t* X, bit_vec_iterator_t* s0bm)
	{
		bit_vec_t V, tB, V2, sB;

		for (size_t i = 0; i < bv_len; ++i)
			X[i] = ~(uint64_t)0;

		auto pc = seq1->data;

		for (size_t i = 0; i < seq1->length; ++i)
		{
			sB = (bit_vec_t)0;
			auto s0b = s0bm[*pc];
			auto pX = X;

			if (*pc++ == UNKNOWN_SYMBOL)				// Unknown aminoacid
				continue;

			for (size_t j = 0; j < bv_len; ++j)
			{
				V = *pX;
				tB = V & *s0b++;
				V2 = V + tB + sB;
				sB = V2 < V;
				*pX++ = V2 | (V - tB);
			}
		}

		for (size_t i = 0; i < bv_len; ++i)
			for (V = ~X[i]; V; V &= V - 1)
				++res[0];
	}



	static void UnrolledCalculate(CSequence *seq0, SeqType *seq1, uint32_t *res, bit_vec_t *X, bit_vec_iterator_t *s0bm)
	{
		bit_vec_t V, tB, V2, sB;

		auto pc = seq1->data;
		auto pX0 = X;

		if (BV_LEN > 0)				*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 1)				*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 2)				*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 3)				*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 4)				*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 5)				*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 6)				*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 7)				*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 8)				*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 9)				*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 10)			*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 11)			*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 12)			*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 13)			*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 14)			*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 15)			*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 16)			*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 17)			*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 18)			*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 19)			*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 20)			*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 21)			*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 22)			*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 23)			*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 24)			*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 25)			*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 26)			*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 27)			*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 28)			*pX0++ = ~(uint64_t)0; 
		if (BV_LEN > 29)			*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 30)			*pX0++ = ~(uint64_t)0;
		if (BV_LEN > 31)			*pX0++ = ~(uint64_t)0;

		for (size_t i = 0; i < seq1->length; ++i)
		{
			sB = (bit_vec_t)0;

			auto pX = X;
			auto s0b = s0bm[*pc];

			if (*pc++ == UNKNOWN_SYMBOL)				// Unknown aminoacid
				continue;

			if (BV_LEN > 0)			CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 1)			CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 2)			CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 3)			CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 4)			CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 5)			CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 6)			CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 7)			CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 8)			CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 9)			CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 10)		CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 11)		CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 12)		CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 13)		CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 14)		CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 15)		CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 16)		CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 17)		CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 18)		CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 19)		CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 20)		CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 21)		CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 22)		CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 23)		CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 24)		CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 25)		CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 26)		CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 27)		CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 28)		CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 29)		CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 30)		CLASSIC_LCS_INNER_LOOP;
			if (BV_LEN > 31)		CLASSIC_LCS_INNER_LOOP;
		}


		auto pX = X;

		if (BV_LEN > 0)		CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 1)		CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 2)		CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 3)		CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 4)		CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 5)		CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 6)		CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 7)		CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 8)		CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 9)		CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 10)	CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 11)	CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 12)	CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 13)	CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 14)	CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 15)	CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 16)	CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 17)	CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 18)	CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 19)	CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 20)	CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 21)	CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 22)	CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 23)	CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 24)	CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 25)	CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 26)	CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 27)	CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 28)	CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 29)	CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 30)	CLASSIC_POP_CNT_LOOP;
		if (BV_LEN > 31)	CLASSIC_POP_CNT_LOOP;
	}
};

#endif
