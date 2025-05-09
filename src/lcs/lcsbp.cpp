#include "../core/sequence.h"
#include "lcsbp.h"
#include "lcsbp_classic.h"

#if defined(SIMD_AVX) || defined(SIMD_AVX2) || defined(SIMD_AVX512)
#include "../simd/lcsbp_avx_intr.h"
#endif

#if defined(SIMD_AVX2) || defined(SIMD_AVX512)
#include "../simd/lcsbp_avx2_intr.h"
#endif

#if defined(SIMD_AVX512)
#include "../simd/lcsbp_avx512_intr.h"
#endif

#if defined(SIMD_NEON)
#include "../simd/lcsbp_neon_intr.h"
#endif

#include <algorithm>

// *******************************************************************
CLCSBP::CLCSBP(instruction_set_t _instruction_set)
{
	instruction_set = _instruction_set;

	lcsbp_classic = std::shared_ptr<CLCSBP_Classic>(new CLCSBP_Classic());	

#if defined(SIMD_AVX) || defined(SIMD_AVX2) || defined(SIMD_AVX512)
	lcsbp_avx_intr = std::shared_ptr<CLCSBP_AVX_INTR>(new CLCSBP_AVX_INTR());
#endif

#if defined(SIMD_AVX2) || defined(SIMD_AVX512)
	lcsbp_avx2_intr = std::shared_ptr<CLCSBP_AVX2_INTR>(new CLCSBP_AVX2_INTR());
#endif

#if defined(SIMD_AVX512)
	lcsbp_avx512_intr = std::shared_ptr<CLCSBP_AVX512_INTR>(new CLCSBP_AVX512_INTR());
#endif

#if defined(SIMD_NEON)
	lcsbp_neon_intr = std::shared_ptr<CLCSBP_NEON_INTR>(new CLCSBP_NEON_INTR());
#endif
}

// *******************************************************************
void CLCSBP::GetLCSBP(CSequence *seq0, CSequence *seq1, CSequence *seq2, CSequence *seq3, CSequence *seq4,
	uint32_t *dist)
{
	if (seq4 == nullptr)
	{
		if (seq1 != nullptr)
			lcsbp_classic->Calculate(seq0, seq1, dist + 0);
		if (seq2 != nullptr)
			lcsbp_classic->Calculate(seq0, seq2, dist + 1);
		if (seq3 != nullptr)
			lcsbp_classic->Calculate(seq0, seq3, dist + 2);
		if (seq4 != nullptr)
			lcsbp_classic->Calculate(seq0, seq4, dist + 3);
	}
	else {
#if defined(SIMD_AVX2) || defined(SIMD_AVX512)
		if (instruction_set < instruction_set_t::avx)				// In theory SSE2 will suffice, but the SSE2-compiled code is too slow
		{
			lcsbp_classic->Calculate(seq0, seq1, dist + 0);
			lcsbp_classic->Calculate(seq0, seq2, dist + 1);
			lcsbp_classic->Calculate(seq0, seq3, dist + 2);
			lcsbp_classic->Calculate(seq0, seq4, dist + 3);
		}
		else if (instruction_set < instruction_set_t::avx2) {
			lcsbp_avx_intr->Calculate(seq0, seq1, seq2, dist + 0);
			lcsbp_avx_intr->Calculate(seq0, seq3, seq4, dist + 2);
		}
		else {
			lcsbp_avx2_intr->Calculate(seq0, seq1, seq2, seq3, seq4, dist);
		}
#elif defined(SIMD_AVX)
		if (instruction_set < instruction_set_t::avx)
		{
			lcsbp_classic->Calculate(seq0, seq1, dist + 0);
			lcsbp_classic->Calculate(seq0, seq2, dist + 1);
			lcsbp_classic->Calculate(seq0, seq3, dist + 2);
			lcsbp_classic->Calculate(seq0, seq4, dist + 3);
		}
		else 
		{
			lcsbp_avx_intr->Calculate(seq0, seq1, seq2, dist + 0);
			lcsbp_avx_intr->Calculate(seq0, seq3, seq4, dist + 2);
		}
#elif defined(SIMD_NEON)
		lcsbp_neon_intr->Calculate(seq0, seq1, seq2, dist + 0);
		lcsbp_neon_intr->Calculate(seq0, seq3, seq4, dist + 2);
#else
		lcsbp_classic->Calculate(seq0, seq1, dist + 0);
		lcsbp_classic->Calculate(seq0, seq2, dist + 1);
		lcsbp_classic->Calculate(seq0, seq3, dist + 2);
		lcsbp_classic->Calculate(seq0, seq4, dist + 3);
#endif
	}
}

// *******************************************************************
void CLCSBP::GetLCSBP(CSequence *seq0, CSequenceView *sv1, CSequenceView *sv2, CSequenceView *sv3, CSequenceView *sv4,
	uint32_t *dist)
{
	if (sv4 == nullptr)
	{
		if (sv1 != nullptr)
			lcsbp_classic->Calculate(seq0, sv1, dist + 0);
		if (sv2 != nullptr)
			lcsbp_classic->Calculate(seq0, sv2, dist + 1);
		if (sv3 != nullptr)
			lcsbp_classic->Calculate(seq0, sv3, dist + 2);
		if (sv4 != nullptr)
			lcsbp_classic->Calculate(seq0, sv4, dist + 3);
	}
	else 
	{
#if defined(SIMD_AVX2) || defined(SIMD_AVX512)
		if (instruction_set < instruction_set_t::avx)
		{
			lcsbp_classic->Calculate(seq0, sv1, dist + 0);
			lcsbp_classic->Calculate(seq0, sv2, dist + 1);
			lcsbp_classic->Calculate(seq0, sv3, dist + 2);
			lcsbp_classic->Calculate(seq0, sv4, dist + 3);
		}
		else if (instruction_set < instruction_set_t::avx2) 
		{
			lcsbp_avx_intr->Calculate(seq0, sv1, sv2, dist + 0);
			lcsbp_avx_intr->Calculate(seq0, sv3, sv4, dist + 2);
		}
		else 
		{
			lcsbp_avx2_intr->Calculate(seq0, sv1, sv2, sv3, sv4, dist);
		}
#elif defined(SIMD_AVX)
		if (instruction_set < instruction_set_t::avx)
		{
			lcsbp_classic->Calculate(seq0, sv1, dist + 0);
			lcsbp_classic->Calculate(seq0, sv2, dist + 1);
			lcsbp_classic->Calculate(seq0, sv3, dist + 2);
			lcsbp_classic->Calculate(seq0, sv4, dist + 3);
		}
		else 
		{
			lcsbp_avx_intr->Calculate(seq0, sv1, sv2, dist + 0);
			lcsbp_avx_intr->Calculate(seq0, sv3, sv4, dist + 2);
		}
#elif defined(SIMD_NEON)
		lcsbp_neon_intr->Calculate(seq0, sv1, sv2, dist + 0);
		lcsbp_neon_intr->Calculate(seq0, sv3, sv4, dist + 2);
#else
		lcsbp_classic->Calculate(seq0, sv1, dist + 0);
		lcsbp_classic->Calculate(seq0, sv2, dist + 1);
		lcsbp_classic->Calculate(seq0, sv3, dist + 2);
		lcsbp_classic->Calculate(seq0, sv4, dist + 3);
#endif
	}
}

// *******************************************************************
void CLCSBP::GetLCSBP(CSequence *seq0, CSequence *seq1, CSequence *seq2, CSequence *seq3, CSequence *seq4, CSequence* seq5, CSequence* seq6, CSequence* seq7, CSequence* seq8,
	uint32_t *dist)
{
	if (seq8 == nullptr)
	{
		if (seq1 != nullptr)			lcsbp_classic->Calculate(seq0, seq1, dist + 0);
		if (seq2 != nullptr)			lcsbp_classic->Calculate(seq0, seq2, dist + 1);
		if (seq3 != nullptr)			lcsbp_classic->Calculate(seq0, seq3, dist + 2);
		if (seq4 != nullptr)			lcsbp_classic->Calculate(seq0, seq4, dist + 3);
		if (seq5 != nullptr)			lcsbp_classic->Calculate(seq0, seq5, dist + 4);
		if (seq6 != nullptr)			lcsbp_classic->Calculate(seq0, seq6, dist + 5);
		if (seq7 != nullptr)			lcsbp_classic->Calculate(seq0, seq7, dist + 6);
		if (seq8 != nullptr)			lcsbp_classic->Calculate(seq0, seq8, dist + 7);
	}
	else 
	{
#if defined(SIMD_AVX512)
		if (instruction_set < instruction_set_t::avx)
		{
			lcsbp_classic->Calculate(seq0, seq1, dist + 0);
			lcsbp_classic->Calculate(seq0, seq2, dist + 1);
			lcsbp_classic->Calculate(seq0, seq3, dist + 2);
			lcsbp_classic->Calculate(seq0, seq4, dist + 3);
			lcsbp_classic->Calculate(seq0, seq5, dist + 4);
			lcsbp_classic->Calculate(seq0, seq6, dist + 5);
			lcsbp_classic->Calculate(seq0, seq7, dist + 6);
			lcsbp_classic->Calculate(seq0, seq8, dist + 7);
		}
		else if (instruction_set < instruction_set_t::avx2) 
		{
			lcsbp_avx_intr->Calculate(seq0, seq1, seq2, dist + 0);
			lcsbp_avx_intr->Calculate(seq0, seq3, seq4, dist + 2);
			lcsbp_avx_intr->Calculate(seq0, seq5, seq6, dist + 4);
			lcsbp_avx_intr->Calculate(seq0, seq7, seq8, dist + 6);
		}
		else if (instruction_set < instruction_set_t::avx512) 
		{
			lcsbp_avx2_intr->Calculate(seq0, seq1, seq2, seq3, seq4, dist);
			lcsbp_avx2_intr->Calculate(seq0, seq5, seq6, seq7, seq8, dist + 4);
		}
		else
			lcsbp_avx512_intr->Calculate(seq0, seq1, seq2, seq3, seq4, seq5, seq6, seq7, seq8, dist);
#elif defined(SIMD_AVX2)
		if (instruction_set < instruction_set_t::avx)
		{
			lcsbp_classic->Calculate(seq0, seq1, dist + 0);
			lcsbp_classic->Calculate(seq0, seq2, dist + 1);
			lcsbp_classic->Calculate(seq0, seq3, dist + 2);
			lcsbp_classic->Calculate(seq0, seq4, dist + 3);
			lcsbp_classic->Calculate(seq0, seq5, dist + 4);
			lcsbp_classic->Calculate(seq0, seq6, dist + 5);
			lcsbp_classic->Calculate(seq0, seq7, dist + 6);
			lcsbp_classic->Calculate(seq0, seq8, dist + 7);
		}
		else if (instruction_set < instruction_set_t::avx2) 
		{
			lcsbp_avx_intr->Calculate(seq0, seq1, seq2, dist + 0);
			lcsbp_avx_intr->Calculate(seq0, seq3, seq4, dist + 2);
			lcsbp_avx_intr->Calculate(seq0, seq5, seq6, dist + 4);
			lcsbp_avx_intr->Calculate(seq0, seq7, seq8, dist + 6);
		}
		else 
		{
			lcsbp_avx2_intr->Calculate(seq0, seq1, seq2, seq3, seq4, dist);
			lcsbp_avx2_intr->Calculate(seq0, seq5, seq6, seq7, seq8, dist + 4);
		}
#elif defined(SIMD_AVX)
		if (instruction_set < instruction_set_t::avx)
		{
			lcsbp_classic->Calculate(seq0, seq1, dist + 0);
			lcsbp_classic->Calculate(seq0, seq2, dist + 1);
			lcsbp_classic->Calculate(seq0, seq3, dist + 2);
			lcsbp_classic->Calculate(seq0, seq4, dist + 3);
			lcsbp_classic->Calculate(seq0, seq5, dist + 4);
			lcsbp_classic->Calculate(seq0, seq6, dist + 5);
			lcsbp_classic->Calculate(seq0, seq7, dist + 6);
			lcsbp_classic->Calculate(seq0, seq8, dist + 7);
		}
		else 
		{
			lcsbp_avx_intr->Calculate(seq0, seq1, seq2, dist + 0);
			lcsbp_avx_intr->Calculate(seq0, seq3, seq4, dist + 2);
			lcsbp_avx_intr->Calculate(seq0, seq5, seq6, dist + 4);
			lcsbp_avx_intr->Calculate(seq0, seq7, seq8, dist + 6);
		}
#elif defined(SIMD_NEON)
		lcsbp_neon_intr->Calculate(seq0, seq1, seq2, dist + 0);
		lcsbp_neon_intr->Calculate(seq0, seq3, seq4, dist + 2);
		lcsbp_neon_intr->Calculate(seq0, seq5, seq6, dist + 4);
		lcsbp_neon_intr->Calculate(seq0, seq7, seq8, dist + 6);
#else
		lcsbp_classic->Calculate(seq0, seq1, dist + 0);
		lcsbp_classic->Calculate(seq0, seq2, dist + 1);
		lcsbp_classic->Calculate(seq0, seq3, dist + 2);
		lcsbp_classic->Calculate(seq0, seq4, dist + 3);
		lcsbp_classic->Calculate(seq0, seq5, dist + 4);
		lcsbp_classic->Calculate(seq0, seq6, dist + 5);
		lcsbp_classic->Calculate(seq0, seq7, dist + 6);
		lcsbp_classic->Calculate(seq0, seq8, dist + 7);
#endif
	}
}

// *******************************************************************
void CLCSBP::GetLCSBP(CSequence *seq0, CSequenceView *sv1, CSequenceView *sv2, CSequenceView *sv3, CSequenceView *sv4, CSequenceView* sv5, CSequenceView* sv6, CSequenceView* sv7, CSequenceView* sv8,
	uint32_t *dist)
{
	if (sv8 == nullptr)
	{
		if (sv1 != nullptr)			lcsbp_classic->Calculate(seq0, sv1, dist + 0);
		if (sv2 != nullptr)			lcsbp_classic->Calculate(seq0, sv2, dist + 1);
		if (sv3 != nullptr)			lcsbp_classic->Calculate(seq0, sv3, dist + 2);
		if (sv4 != nullptr)			lcsbp_classic->Calculate(seq0, sv4, dist + 3);
		if (sv5 != nullptr)			lcsbp_classic->Calculate(seq0, sv5, dist + 4);
		if (sv6 != nullptr)			lcsbp_classic->Calculate(seq0, sv6, dist + 5);
		if (sv7 != nullptr)			lcsbp_classic->Calculate(seq0, sv7, dist + 6);
		if (sv8 != nullptr)			lcsbp_classic->Calculate(seq0, sv8, dist + 7);
	}
	else
	{
#if defined(SIMD_AVX512)
		if (instruction_set < instruction_set_t::avx)		
		{
			lcsbp_classic->Calculate(seq0, sv1, dist + 0);
			lcsbp_classic->Calculate(seq0, sv2, dist + 1);
			lcsbp_classic->Calculate(seq0, sv3, dist + 2);
			lcsbp_classic->Calculate(seq0, sv4, dist + 3);
			lcsbp_classic->Calculate(seq0, sv5, dist + 4);
			lcsbp_classic->Calculate(seq0, sv6, dist + 5);
			lcsbp_classic->Calculate(seq0, sv7, dist + 6);
			lcsbp_classic->Calculate(seq0, sv8, dist + 7);
		}
		else if (instruction_set < instruction_set_t::avx2) 
		{
			lcsbp_avx_intr->Calculate(seq0, sv1, sv2, dist + 0);
			lcsbp_avx_intr->Calculate(seq0, sv3, sv4, dist + 2);
			lcsbp_avx_intr->Calculate(seq0, sv5, sv6, dist + 4);
			lcsbp_avx_intr->Calculate(seq0, sv7, sv8, dist + 6);
		}
		else if (instruction_set < instruction_set_t::avx512) 
		{
			lcsbp_avx2_intr->Calculate(seq0, sv1, sv2, sv3, sv4, dist);
			lcsbp_avx2_intr->Calculate(seq0, sv5, sv6, sv7, sv8, dist + 4);
		}
		else
			lcsbp_avx512_intr->Calculate(seq0, sv1, sv2, sv3, sv4, sv5, sv6, sv7, sv8, dist);
#elif defined(SIMD_AVX2)
		if (instruction_set < instruction_set_t::avx)		
		{
			lcsbp_classic->Calculate(seq0, sv1, dist + 0);
			lcsbp_classic->Calculate(seq0, sv2, dist + 1);
			lcsbp_classic->Calculate(seq0, sv3, dist + 2);
			lcsbp_classic->Calculate(seq0, sv4, dist + 3);
			lcsbp_classic->Calculate(seq0, sv5, dist + 4);
			lcsbp_classic->Calculate(seq0, sv6, dist + 5);
			lcsbp_classic->Calculate(seq0, sv7, dist + 6);
			lcsbp_classic->Calculate(seq0, sv8, dist + 7);
		}
		else if (instruction_set < instruction_set_t::avx2) 
		{
			lcsbp_avx_intr->Calculate(seq0, sv1, sv2, dist + 0);
			lcsbp_avx_intr->Calculate(seq0, sv3, sv4, dist + 2);
			lcsbp_avx_intr->Calculate(seq0, sv5, sv6, dist + 4);
			lcsbp_avx_intr->Calculate(seq0, sv7, sv8, dist + 6);
		}
		else 
		{
			lcsbp_avx2_intr->Calculate(seq0, sv1, sv2, sv3, sv4, dist);
			lcsbp_avx2_intr->Calculate(seq0, sv5, sv6, sv7, sv8, dist + 4);
		}
#elif defined(SIMD_AVX)
		if (instruction_set < instruction_set_t::avx)		
		{
			lcsbp_classic->Calculate(seq0, sv1, dist + 0);
			lcsbp_classic->Calculate(seq0, sv2, dist + 1);
			lcsbp_classic->Calculate(seq0, sv3, dist + 2);
			lcsbp_classic->Calculate(seq0, sv4, dist + 3);
			lcsbp_classic->Calculate(seq0, sv5, dist + 4);
			lcsbp_classic->Calculate(seq0, sv6, dist + 5);
			lcsbp_classic->Calculate(seq0, sv7, dist + 6);
			lcsbp_classic->Calculate(seq0, sv8, dist + 7);
		}
		else 
		{
			lcsbp_avx_intr->Calculate(seq0, sv1, sv2, dist + 0);
			lcsbp_avx_intr->Calculate(seq0, sv3, sv4, dist + 2);
			lcsbp_avx_intr->Calculate(seq0, sv5, sv6, dist + 4);
			lcsbp_avx_intr->Calculate(seq0, sv7, sv8, dist + 6);
		}
#elif defined(SIMD_NEON)
		lcsbp_neon_intr->Calculate(seq0, sv1, sv2, dist + 0);
		lcsbp_neon_intr->Calculate(seq0, sv3, sv4, dist + 2);
		lcsbp_neon_intr->Calculate(seq0, sv5, sv6, dist + 4);
		lcsbp_neon_intr->Calculate(seq0, sv7, sv8, dist + 6);
#else
		lcsbp_classic->Calculate(seq0, sv1, dist + 0);
		lcsbp_classic->Calculate(seq0, sv2, dist + 1);
		lcsbp_classic->Calculate(seq0, sv3, dist + 2);
		lcsbp_classic->Calculate(seq0, sv4, dist + 3);
		lcsbp_classic->Calculate(seq0, sv5, dist + 4);
		lcsbp_classic->Calculate(seq0, sv6, dist + 5);
		lcsbp_classic->Calculate(seq0, sv7, dist + 6);
		lcsbp_classic->Calculate(seq0, sv8, dist + 7);
#endif
	}
}

#ifdef DEVELOPER_MODE
// *******************************************************************
// Compute LCS length for two sequences in the classical way - just for development
double CLCSBP::GetLCS(CSequence &seq1, CSequence &seq2)
{
	int **dp_row = new int*[2];

	for (int i = 0; i < 2; ++i)
		dp_row[i] = new int[seq2.length + 1];

	fill(dp_row[0], dp_row[0] + seq2.length + 1, 0);

	for (int i = 1; i <= (int)seq1.length; ++i)
	{
		int ii = i % 2;
		dp_row[ii][0] = 0;
		for (int j = 1; j <= (int)seq2.length; ++j)
			if (seq1.data[i - 1] == seq2.data[j - 1])
				dp_row[ii][j] = dp_row[!ii][j - 1] + 1;
			else
				dp_row[ii][j] = max(dp_row[ii][j - 1], dp_row[!ii][j]);
	}

	return dp_row[seq1.length % 2][seq2.length];
}
#endif
