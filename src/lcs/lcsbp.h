#ifndef _LSCBP_H
#define _LSCBP_H
	
#include <memory>
#include "../core/defs.h"

class CLCSBP_Classic;
class CLCSBP_AVX_INTR;
class CLCSBP_AVX2_INTR;
class CLCSBP_NEON_INTR;

class CSequence;
struct CSequenceView;


class CLCSBP
{
	instruction_set_t instruction_set;

	std::shared_ptr<CLCSBP_Classic> lcsbp_classic;
	std::shared_ptr<CLCSBP_AVX_INTR> lcsbp_avx_intr;
	std::shared_ptr<CLCSBP_AVX2_INTR> lcsbp_avx2_intr;
	std::shared_ptr<CLCSBP_NEON_INTR> lcsbp_neon_intr;

public:
	CLCSBP(instruction_set_t _instruction_set = instruction_set_t::none);
	
	void GetLCSBP(CSequence *seq0, CSequence *seq1, CSequence *seq2, CSequence *seq3, CSequence *seq4,
		uint32_t *dist);
	void GetLCSBP(CSequence *seq0, CSequenceView *sv1, CSequenceView *sv2, CSequenceView *sv3, CSequenceView *sv4,
		uint32_t *dist);

#ifdef DEVELOPER_MODE
	double GetLCS(CSequence &seq1, CSequence &seq2);
#endif
};

#endif