
#include "../core/lcsbp.h"
#include "../core/lcsbp_classic.h"

#ifndef NO_AVX
#include "../core/lcsbp_avx.h"
#ifndef NO_AVX2
#include "../core/lcsbp_avx2.h"
#endif
#endif


CLCSBP::CLCSBP(instruction_set_t _instruction_set)
{
	instruction_set = _instruction_set;

	lcsbp_classic = std::shared_ptr<CLCSBP_Classic>(new CLCSBP_Classic());	
#ifndef NO_AVX
	lcsbp_avx = std::shared_ptr<CLCSBP_AVX>(new CLCSBP_AVX());
#ifndef NO_AVX2
	lcsbp_avx2 = std::shared_ptr<CLCSBP_AVX2>(new CLCSBP_AVX2());
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
				lcsbp_avx->Calculate(seq0, seq1, seq2, dist1, dist2);
				lcsbp_avx->Calculate(seq0, seq3, seq4, dist3, dist4);
			}
			else
			{
#ifndef NO_AVX2
				lcsbp_avx2->Calculate(seq0, seq1, seq2, seq3, seq4, dist1, dist2, dist3, dist4);
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