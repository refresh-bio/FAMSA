/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#pragma once

#include "AbstractTreeGenerator.h"
#include "IPartialGenerator.h"

#include <vector>
#include <utility>

#define SLINK_HANDLE_TIES

#ifdef SLINK_HANDLE_TIES
struct slink_dist_t {
	double first;
	uint64_t second;

	// this is to preserve consistency with smilarity variant
	// - increasingly by distance
	// - decreasingly by id
	bool operator<(const slink_dist_t& rhs) const {
		return (this->first == rhs.first)
			? (this->second > rhs.second)
			: (this->first < rhs.first);
	}

	bool operator<=(const slink_dist_t& rhs) const {
		return (this->first == rhs.first)
			? (this->second >= rhs.second)
			: (this->first <= rhs.first);
	}

};

//using slink_similarity_t = pair<double, uint64_t>;
#else
using slink_dist_t = double;
#endif



template <Distance measure>
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

	SingleLinkage(int n_threads, instruction_set_t instruction_set) 
		: AbstractTreeGenerator(n_threads, instruction_set) {}

	void run(std::vector<CSequence*>& sequences, tree_structure& tree) override;

	void runPartial(std::vector<CSequence*>& sequences, tree_structure& tree) override;
};

