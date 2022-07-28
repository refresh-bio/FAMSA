/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#pragma once

#include "AbstractTreeGenerator.h"
#include "IPartialGenerator.h"
#include "Clustering.h"

#include <memory>

template <Distance _distance>
class FastTree : public AbstractTreeGenerator {
public:

	FastTree(
		int n_threads, 
		instruction_set_t instruction_set,
		std::shared_ptr<IPartialGenerator> partialGenerator, 
		int subtreeSize,
		std::shared_ptr<IClustering> clustering,
		int sampleSize);

	virtual void run(std::vector<CSequence*>& sequences, tree_structure& tree) override;

protected:
	std::shared_ptr<IPartialGenerator> partialGenerator;
	int subtreeSize;
	std::shared_ptr<IClustering> clustering;
	int sampleSize;
	int clusteringThreshold;

	void doStep(std::vector<CSequence*>& sequences, tree_structure& tree, int previousTop, bool parallel);

	int randomSeeds(
		std::vector<CSequence*>& sequences,
		int n_seeds,
		int * seed_ids,
		float * similarity_row);

	int clusterSeeds(
		std::vector<CSequence*>& sequences,
		int n_seeds,
		int n_samples,
		int * seed_ids,
		float * similarity_row);
};