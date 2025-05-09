#ifndef _LSCBP_H
#define _LSCBP_H
	
#include <memory>
#include "../core/defs.h"

class CLCSBP_Classic;
class CLCSBP_AVX_INTR;
class CLCSBP_AVX2_INTR;
class CLCSBP_AVX512_INTR;
class CLCSBP_NEON_INTR;

class CSequence;
struct CSequenceView;

class CLCSBP
{
	instruction_set_t instruction_set;

	std::shared_ptr<CLCSBP_Classic> lcsbp_classic;

#if defined(SIMD_AVX) || defined(SIMD_AVX2) || defined(SIMD_AVX512)
	std::shared_ptr<CLCSBP_AVX_INTR> lcsbp_avx_intr;
#endif
#if defined(SIMD_AVX2) || defined(SIMD_AVX512)
	std::shared_ptr<CLCSBP_AVX2_INTR> lcsbp_avx2_intr;
#endif
#if defined(SIMD_AVX512)
	std::shared_ptr<CLCSBP_AVX512_INTR> lcsbp_avx512_intr;
#endif
#if defined(SIMD_NEON)
	std::shared_ptr<CLCSBP_NEON_INTR> lcsbp_neon_intr;
#endif

public:
	CLCSBP(instruction_set_t _instruction_set = instruction_set_t::none);
	
	void GetLCSBP(CSequence *seq0, CSequence *seq1, CSequence *seq2, CSequence *seq3, CSequence *seq4,
		uint32_t *dist);
	void GetLCSBP(CSequence *seq0, CSequenceView *sv1, CSequenceView *sv2, CSequenceView *sv3, CSequenceView *sv4,
		uint32_t *dist);

	void GetLCSBP(CSequence *seq0, CSequence *seq1, CSequence *seq2, CSequence *seq3, CSequence *seq4, CSequence* seq5, CSequence* seq6, CSequence* seq7, CSequence* seq8,
		uint32_t *dist);
	void GetLCSBP(CSequence *seq0, CSequenceView *sv1, CSequenceView *sv2, CSequenceView *sv3, CSequenceView *sv4, CSequenceView* sv5, CSequenceView* sv6, CSequenceView* sv7, CSequenceView* sv8,
		uint32_t *dist);

#ifdef DEVELOPER_MODE
	double GetLCS(CSequence &seq1, CSequence &seq2);
#endif
};

#endif