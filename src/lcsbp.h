#ifndef _LSCBP_H
#define _LSCBP_H

#include <memory>
#include "defs.h"

class CLCSBP_Classic;
class CLCSBP_AVX_INTR;
class CLCSBP_AVX2_INTR;

class CSequence;


class CLCSBP
{
	instruction_set_t instruction_set;

	std::shared_ptr<CLCSBP_Classic> lcsbp_classic;
	std::shared_ptr<CLCSBP_AVX_INTR> lcsbp_avx_intr;
	std::shared_ptr<CLCSBP_AVX2_INTR> lcsbp_avx2_intr;

public:
	CLCSBP(instruction_set_t _instruction_set = instruction_set_t::none);

	void GetLCSBP(CSequence *seq0, CSequence *seq1, CSequence *seq2, CSequence *seq3, CSequence *seq4,
		uint32_t &dist1, uint32_t &dist2, uint32_t &dist3, uint32_t &dist4);

#ifdef DEVELOPER_MODE
	double GetLCS(CSequence &seq1, CSequence &seq2);
#endif
};


#endif
