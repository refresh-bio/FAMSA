/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#ifndef _LCSBP_AVX_H
#define _LCSBP_AVX_H

#include "../libs/vectorclass.h"
#include "../core/sequence.h"
#include "../core/meta_oper.h"

template <unsigned BV_LEN> class CLCSBP_AVX_Impl;

class CLCSBP_AVX
{
public:
	void *raw_X;
	void *orig_X;
	size_t X_size;
	size_t raw_X_size;
	Vec2uq *X;

	inline void prepare_X(size_t bv_len);
	void calculate(CSequence *seq0, CSequence *seq1, CSequence *seq2, uint32_t *res, uint32_t bv_len, uint32_t max_len);

public:
	CLCSBP_AVX() {
		X = nullptr;
		raw_X = nullptr;
		orig_X = nullptr;
		X_size = 0;
		raw_X_size = 0;
	};

	~CLCSBP_AVX()
	{
		if (orig_X)
			free(orig_X);
	}
	
	void Calculate(CSequence *seq0, CSequence *seq1, CSequence *seq2, uint32_t &dist1, uint32_t &dist2);
};

template <unsigned BV_LEN> class CLCSBP_AVX_Impl {
public:
	static void CalculateOld(CSequence *seq0, CSequence *seq1, CSequence *seq2, uint32_t *res, uint32_t max_len, Vec2uq *X)
	{
		Vec2uq V, tB, V2, sB;
		auto pX0 = X;

		for (int i = 0; i < BV_LEN; ++i)
			*pX0++ = Vec2uq(~(uint64_t)0);

		for (size_t i = 0; i < max_len; ++i)
		{
			sB = 0;
			symbol_t c1 = seq1->data[i];
			symbol_t c2 = seq2->data[i];

			auto pX = X;
			auto pbm1 = seq0->bit_masks[c1].begin();
			auto pbm2 = seq0->bit_masks[c2].begin();

			for (int j = 0; j < BV_LEN; ++j)
			{
				V = *pX;
				tB = V & Vec2uq(*pbm1++, *pbm2++);
				V2 = V + tB - sB;			// Important: true in boolean vectors is represented as -1, so here we use "-sB" instead of "+sB" in the classic code
				sB = V2 < V;
				*pX++ = V2 | (V - tB);
			}
		}


#ifdef _MSC_VER					// Visual C++
		auto pX = X;
		IterFwd([&](const int i){
			res[0] += __popcnt64(~(*pX)[0]);
			res[1] += __popcnt64(~(*pX)[1]);
			++pX;
		}, uint_<BV_LEN - 1>());
#else
#ifdef __GNUC__
		for(int i = 0; i < BV_LEN; ++i)
		{
			res[0] += __builtin_popcountll(~X[i][0]);
			res[1] += __builtin_popcountll(~X[i][1]);
		}
#else
		for (size_t i = 0; i < BV_LEN; ++i)
			for (int v = 0; v < 2; ++v)
				for (uint64_t T = ~X[i][v]; T; T &= T - 1)
					++res[v];
#endif
#endif
	}

#define AVX_POP_CNT_LOOP					\
	{										\
		res[0] += POPCNT(~(*pX)[0]);	\
		res[1] += POPCNT(~(*pX)[1]);	\
		++pX;								\
	}

	// Important: true in boolean vectors is represented as -1, so here we use "-sB" instead of "+sB" in the classic code
#define AVX_LCS_INNER_LOOP					\
	{										\
		V = *pX;							\
		tB = V & Vec2uq(*pbm1++, *pbm2++);	\
		V2 = V + tB - sB;					\
		sB = V2 < V;						\
		*pX++ = V2 | (V - tB);				\
	}	


	static void Calculate(CSequence *seq0, CSequence *seq1, CSequence *seq2, uint32_t *res, uint32_t max_len, Vec2uq *X)
	{
		Vec2uq V, tB, V2, sB;
		auto pX0 = X;

		if (BV_LEN > 0)				*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 1)				*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 2)				*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 3)				*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 4)				*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 5)				*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 6)				*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 7)				*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 8)				*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 9)				*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 10)			*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 11)			*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 12)			*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 13)			*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 14)			*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 15)			*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 16)			*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 17)			*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 18)			*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 19)			*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 20)			*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 21)			*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 22)			*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 23)			*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 24)			*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 25)			*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 26)			*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 27)			*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 28)			*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 29)			*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 30)			*pX0++ = Vec2uq(~(uint64_t)0);
		if (BV_LEN > 31)			*pX0++ = Vec2uq(~(uint64_t)0);

		for (size_t i = 0; i < max_len; ++i)
		{
			sB = 0;
			symbol_t c1 = seq1->data[i];
			symbol_t c2 = seq2->data[i];

			auto pX = X;
			auto pbm1 = seq0->bit_masks[c1].begin();
			auto pbm2 = seq0->bit_masks[c2].begin();

			if (BV_LEN > 0)			AVX_LCS_INNER_LOOP;
			if (BV_LEN > 1)			AVX_LCS_INNER_LOOP;
			if (BV_LEN > 2)			AVX_LCS_INNER_LOOP;
			if (BV_LEN > 3)			AVX_LCS_INNER_LOOP;
			if (BV_LEN > 4)			AVX_LCS_INNER_LOOP;
			if (BV_LEN > 5)			AVX_LCS_INNER_LOOP;
			if (BV_LEN > 6)			AVX_LCS_INNER_LOOP;
			if (BV_LEN > 7)			AVX_LCS_INNER_LOOP;
			if (BV_LEN > 8)			AVX_LCS_INNER_LOOP;
			if (BV_LEN > 9)			AVX_LCS_INNER_LOOP;
			if (BV_LEN > 10)		AVX_LCS_INNER_LOOP;
			if (BV_LEN > 11)		AVX_LCS_INNER_LOOP;
			if (BV_LEN > 12)		AVX_LCS_INNER_LOOP;
			if (BV_LEN > 13)		AVX_LCS_INNER_LOOP;
			if (BV_LEN > 14)		AVX_LCS_INNER_LOOP;
			if (BV_LEN > 15)		AVX_LCS_INNER_LOOP;
			if (BV_LEN > 16)		AVX_LCS_INNER_LOOP;
			if (BV_LEN > 17)		AVX_LCS_INNER_LOOP;
			if (BV_LEN > 18)		AVX_LCS_INNER_LOOP;
			if (BV_LEN > 19)		AVX_LCS_INNER_LOOP;
			if (BV_LEN > 20)		AVX_LCS_INNER_LOOP;
			if (BV_LEN > 21)		AVX_LCS_INNER_LOOP;
			if (BV_LEN > 22)		AVX_LCS_INNER_LOOP;
			if (BV_LEN > 23)		AVX_LCS_INNER_LOOP;
			if (BV_LEN > 24)		AVX_LCS_INNER_LOOP;
			if (BV_LEN > 25)		AVX_LCS_INNER_LOOP;
			if (BV_LEN > 26)		AVX_LCS_INNER_LOOP;
			if (BV_LEN > 27)		AVX_LCS_INNER_LOOP;
			if (BV_LEN > 28)		AVX_LCS_INNER_LOOP;
			if (BV_LEN > 29)		AVX_LCS_INNER_LOOP;
			if (BV_LEN > 30)		AVX_LCS_INNER_LOOP;
			if (BV_LEN > 31)		AVX_LCS_INNER_LOOP;
		}

		auto pX = X;
		if (BV_LEN > 0)		AVX_POP_CNT_LOOP;
		if (BV_LEN > 1)		AVX_POP_CNT_LOOP;
		if (BV_LEN > 2)		AVX_POP_CNT_LOOP;
		if (BV_LEN > 3)		AVX_POP_CNT_LOOP;
		if (BV_LEN > 4)		AVX_POP_CNT_LOOP;
		if (BV_LEN > 5)		AVX_POP_CNT_LOOP;
		if (BV_LEN > 6)		AVX_POP_CNT_LOOP;
		if (BV_LEN > 7)		AVX_POP_CNT_LOOP;
		if (BV_LEN > 8)		AVX_POP_CNT_LOOP;
		if (BV_LEN > 9)		AVX_POP_CNT_LOOP;
		if (BV_LEN > 10)	AVX_POP_CNT_LOOP;
		if (BV_LEN > 11)	AVX_POP_CNT_LOOP;
		if (BV_LEN > 12)	AVX_POP_CNT_LOOP;
		if (BV_LEN > 13)	AVX_POP_CNT_LOOP;
		if (BV_LEN > 14)	AVX_POP_CNT_LOOP;
		if (BV_LEN > 15)	AVX_POP_CNT_LOOP;
		if (BV_LEN > 16)	AVX_POP_CNT_LOOP;
		if (BV_LEN > 17)	AVX_POP_CNT_LOOP;
		if (BV_LEN > 18)	AVX_POP_CNT_LOOP;
		if (BV_LEN > 19)	AVX_POP_CNT_LOOP;
		if (BV_LEN > 20)	AVX_POP_CNT_LOOP;
		if (BV_LEN > 21)	AVX_POP_CNT_LOOP;
		if (BV_LEN > 22)	AVX_POP_CNT_LOOP;
		if (BV_LEN > 23)	AVX_POP_CNT_LOOP;
		if (BV_LEN > 24)	AVX_POP_CNT_LOOP;
		if (BV_LEN > 25)	AVX_POP_CNT_LOOP;
		if (BV_LEN > 26)	AVX_POP_CNT_LOOP;
		if (BV_LEN > 27)	AVX_POP_CNT_LOOP;
		if (BV_LEN > 28)	AVX_POP_CNT_LOOP;
		if (BV_LEN > 29)	AVX_POP_CNT_LOOP;
		if (BV_LEN > 30)	AVX_POP_CNT_LOOP;
		if (BV_LEN > 31)	AVX_POP_CNT_LOOP;
	}
};

#endif
