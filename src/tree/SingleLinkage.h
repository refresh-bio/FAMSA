/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#pragma once

#include "AbstractTreeGenerator.h"
#include "IPartialGenerator.h"

#include <vector>
#include <utility>


#define SLINK_HANDLE_TIES

#ifdef SLINK_HANDLE_TIES
using slink_similarity_t = pair<double, uint64_t>;
#else
using slink_similarity_t = double;
#endif


class SingleLinkage : public AbstractTreeGenerator, public IPartialGenerator {
	uint64_t ids_to_uint64(int id1, int id2)
	{
		if (id1 < 0 || id2 < 0)
			return 0u;
		if (id1 > id2)
			return (((uint64_t)id2) << 32) + (uint64_t)id1;
		return (((uint64_t)id1) << 32) + (uint64_t)id2;
	}

public:

	SingleLinkage(double indel_exp, size_t n_threads) : AbstractTreeGenerator(indel_exp, n_threads) {}

	void run(std::vector<CSequence>& sequences, tree_structure& tree) override;

	void runPartial(std::vector<CSequence*>& sequences, tree_structure& tree) override;
};

