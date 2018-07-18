/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#ifndef _GPUMSA_H
#define _GPUMSA_H

#include "../core/msa.h"
#include "../opencl_utils/hardware/OpenCl.h"
#include "../opencl_utils/hardware/Kernel.h"

#include <memory>
#include <vector>

#define SEQS_AS_CHARS

#ifdef SEQS_AS_CHARS
	typedef cl_char type_seq;	
#else
	typedef cl_int type_seq;
#endif

typedef struct {
	cl_int x;
	cl_int y;
} Task;

class CGpuFAMSA : public CFAMSA 
{
public:
	CGpuFAMSA(int platform, int device, int pairsPerWave, int pairsPerTask, int threadsPerPair);

protected:
	const int pairsPerWave;
	
	const int threadsPerPair;
	const int pairsPerTask;
	bool localBitVectors;
	
	std::shared_ptr<clex::OpenCL> cl;
	std::shared_ptr<clex::Kernel> kernelLCS;

	std::shared_ptr<clex::Buffer> buffer_rawSequences;
	std::shared_ptr<clex::Buffer> buffer_lengths;
	std::shared_ptr<clex::Buffer> buffer_offsets;
	std::shared_ptr<clex::Buffer> buffer_bitVectors;
	std::shared_ptr<clex::Buffer> buffer_tasks;
	std::shared_ptr<clex::Buffer> buffer_lcs;
	
	std::vector<cl_int> lengths;
	std::vector<cl_int> offsets;
	std::vector<type_seq> rawSequences;
	std::vector<cl_uint> bitVectors;
	
	virtual void UPGMA_CalculateDistances(UPGMA_dist_t *dist_matrix);

	virtual void SingleLinkage();

	virtual void preprocess();
	virtual void calculateLCSs(
		::size_t verticalRangeLo, 
		::size_t verticalRangeHi, 
		std::vector<int>& lcs, 
		std::thread& postprocessingThread,
		cl::Event& finishedEvent);

	void calculateLCSs(
		::size_t verticalRangeLo,
		::size_t verticalRangeHi,
		std::vector<int>& lcs);

	void calculateReferenceLCSs(
		::size_t verticalRangeLo,
		::size_t verticalRangeHi,
		std::vector<int>& lcs);

	virtual void partialLinkage(
		const std::vector<int> &lcs,
		::size_t verticalRangeLo,
		::size_t verticalRangeHi,
		vector<int>& pi,
		vector<double>& lambda,
		vector<double>& sim_vector) const;
	
	inline ::size_t calculateRangeHi(::size_t rangeLo, ::size_t elements) const;

	inline ::size_t calculateRangeElements(::size_t rangeLo, ::size_t rangeHi) const;
	
};

#endif
