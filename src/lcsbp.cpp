#include "sequence.h"
#include "lcsbp.h"
#include "lcsbp_classic.h"

#ifndef NO_AVX
#include "lcsbp_avx_intr.h"
#ifndef NO_AVX2
#include "lcsbp_avx2_intr.h"
#endif
#endif


CLCSBP::CLCSBP(instruction_set_t _instruction_set)
{
	instruction_set = _instruction_set;

	lcsbp_classic = std::shared_ptr<CLCSBP_Classic>(new CLCSBP_Classic());
#ifndef NO_AVX
	lcsbp_avx_intr = std::shared_ptr<CLCSBP_AVX_INTR>(new CLCSBP_AVX_INTR());
#ifndef NO_AVX2
	lcsbp_avx2_intr = std::shared_ptr<CLCSBP_AVX2_INTR>(new CLCSBP_AVX2_INTR());
#endif
#endif


}

void CLCSBP::GetLCSBP(CSequence *seq0, CSequence *seq1, CSequence *seq2, CSequence *seq3, CSequence *seq4,
	uint32_t &dist1, uint32_t &dist2, uint32_t &dist3, uint32_t &dist4)
{
	if (seq4 == nullptr)
	{
		if (seq1 != nullptr)
			lcsbp_classic->Calculate(seq0, seq1, dist1);
		if (seq2 != nullptr)
			lcsbp_classic->Calculate(seq0, seq2, dist2);
		if (seq3 != nullptr)
			lcsbp_classic->Calculate(seq0, seq3, dist3);
		if (seq4 != nullptr)
			lcsbp_classic->Calculate(seq0, seq4, dist4);
	}
	else {

		if (instruction_set < instruction_set_t::avx)				// In theory SSE2 will suffice, but the SSE2-compiled code is too slow
		{
			lcsbp_classic->Calculate(seq0, seq1, dist1);
			lcsbp_classic->Calculate(seq0, seq2, dist2);
			lcsbp_classic->Calculate(seq0, seq3, dist3);
			lcsbp_classic->Calculate(seq0, seq4, dist4);
		}
		else {
#ifndef NO_AVX
			if (instruction_set < instruction_set_t::avx2)
			{
//				lcsbp_avx->Calculate(seq0, seq1, seq2, dist1, dist2);
//				lcsbp_avx->Calculate(seq0, seq3, seq4, dist3, dist4);
				lcsbp_avx_intr->Calculate(seq0, seq1, seq2, dist1, dist2);
				lcsbp_avx_intr->Calculate(seq0, seq3, seq4, dist3, dist4);
			}
			else
			{
#ifndef NO_AVX2
//				lcsbp_avx2->Calculate(seq0, seq1, seq2, seq3, seq4, dist1, dist2, dist3, dist4);
				lcsbp_avx2_intr->Calculate(seq0, seq1, seq2, seq3, seq4, dist1, dist2, dist3, dist4);
#else
				lcsbp_avx->Calculate(seq0, seq1, seq2, dist1, dist2);
				lcsbp_avx->Calculate(seq0, seq3, seq4, dist3, dist4);
#endif
			}

#else
			lcsbp_classic->Calculate(seq0, seq1, dist1);
			lcsbp_classic->Calculate(seq0, seq2, dist2);
			lcsbp_classic->Calculate(seq0, seq3, dist3);
			lcsbp_classic->Calculate(seq0, seq4, dist4);
#endif
		}
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
