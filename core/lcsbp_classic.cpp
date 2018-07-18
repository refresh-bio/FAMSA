/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#include "../core/lcsbp_classic.h"
#include "../libs/vectorclass.h"

// *******************************************************************
// Prepares (if necessary sufficient amount of memory for LCS calculation
void CLCSBP_Classic::prepare_X(size_t bv_len)
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

	for (int i = 0; i < NO_SYMBOLS; ++i)
		s0bm[i] = seq0->bit_masks[i].begin();
}

// *******************************************************************
void CLCSBP_Classic::calculate(CSequence *seq0, CSequence *seq1, uint32_t *res, uint32_t bv_len)
{
	bit_vec_t V, tB, V2, sB;

	for (size_t i = 0; i < bv_len; ++i)
		X[i] = ~(uint64_t)0;

	auto pc = seq1->data.begin();

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

// *******************************************************************
void CLCSBP_Classic::Calculate(CSequence *seq0, CSequence *seq1,
	uint32_t &dist1)
{
	size_t bv_len = (seq0->length + bv_size - 1) / bv_size;

	prepare_X(bv_len);
	prefetch_bitmasks(seq0);

	uint32_t res[1] = { 0 };

	switch (bv_len)
	{
	case 1:	CLCSBP_Classic_Impl<1>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 2:	CLCSBP_Classic_Impl<2>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 3:	CLCSBP_Classic_Impl<3>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 4:	CLCSBP_Classic_Impl<4>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 5:	CLCSBP_Classic_Impl<5>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 6:	CLCSBP_Classic_Impl<6>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 7:	CLCSBP_Classic_Impl<7>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 8:	CLCSBP_Classic_Impl<8>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 9:	CLCSBP_Classic_Impl<9>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 10:	CLCSBP_Classic_Impl<10>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 11:	CLCSBP_Classic_Impl<11>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 12:	CLCSBP_Classic_Impl<12>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 13:	CLCSBP_Classic_Impl<13>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 14:	CLCSBP_Classic_Impl<14>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 15:	CLCSBP_Classic_Impl<15>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 16:	CLCSBP_Classic_Impl<16>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 17:	CLCSBP_Classic_Impl<17>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 18:	CLCSBP_Classic_Impl<18>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 19:	CLCSBP_Classic_Impl<19>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 20:	CLCSBP_Classic_Impl<20>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 21:	CLCSBP_Classic_Impl<21>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 22:	CLCSBP_Classic_Impl<22>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 23:	CLCSBP_Classic_Impl<23>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 24:	CLCSBP_Classic_Impl<24>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 25:	CLCSBP_Classic_Impl<25>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 26:	CLCSBP_Classic_Impl<26>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 27:	CLCSBP_Classic_Impl<27>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 28:	CLCSBP_Classic_Impl<28>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 29:	CLCSBP_Classic_Impl<29>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 30:	CLCSBP_Classic_Impl<30>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 31:	CLCSBP_Classic_Impl<31>::Calculate(seq0, seq1, res, X, s0bm);					break;
	case 32:	CLCSBP_Classic_Impl<32>::Calculate(seq0, seq1, res, X, s0bm);					break;
	default: 
		calculate(seq0, seq1, res, bv_len);	
	}
	
	dist1 = res[0];
}

