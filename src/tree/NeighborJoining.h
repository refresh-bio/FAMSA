#pragma once
#include "IPartialGenerator.h"
#include "AbstractTreeGenerator.h"

template <Distance _distance>
class NeighborJoining : public AbstractTreeGenerator, public IPartialGenerator {
public:
	
	NeighborJoining(size_t n_threads) : AbstractTreeGenerator(n_threads) {}

	void run(std::vector<CSequence>& sequences, tree_structure& tree) override;
	
	void runPartial(std::vector<CSequence*>& sequences, tree_structure& tree) override;

protected:
	void computeTree(float* distances, size_t n_seq, tree_structure& tree);
};
