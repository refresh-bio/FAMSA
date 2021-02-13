#pragma once
#include "IPartialGenerator.h"
#include "AbstractTreeGenerator.h"

class NeighborJoining : public AbstractTreeGenerator, public IPartialGenerator {
public:
	
	NeighborJoining(double indel_exp, size_t n_threads) : AbstractTreeGenerator(indel_exp, n_threads) {}

	void run(std::vector<CSequence>& sequences, tree_structure& tree) override;
	
	void runPartial(std::vector<CSequence*>& sequences, tree_structure& tree) override;

protected:
	void computeTree(float* distances, size_t n_seq, tree_structure& tree);
};
