#pragma once
#include "IPartialGenerator.h"
#include "AbstractTreeGenerator.h"

template <Distance _distance>
class NeighborJoining : public AbstractTreeGenerator, public IPartialGenerator {
public:
	
	NeighborJoining(int n_threads, instruction_set_t instruction_set) 
		: AbstractTreeGenerator(n_threads, instruction_set) {}

	void run(std::vector<CSequence*>& sequences, tree_structure& tree) override;
	
	void runPartial(std::vector<CSequence*>& sequences, tree_structure& tree) override;
	
protected:
	void computeTree(float* distances, int n_seq, tree_structure& tree);
};
