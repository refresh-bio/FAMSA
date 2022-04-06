/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

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
		size_t n_threads, 
		std::shared_ptr<IPartialGenerator> partialGenerator, 
		int subtreeSize,
		std::shared_ptr<IClustering> clustering,
		int sampleSize);

	virtual void run(std::vector<CSequence>& sequences, tree_structure& tree) override;

protected:
	std::shared_ptr<IPartialGenerator> partialGenerator;
	int subtreeSize;
	std::shared_ptr<IClustering> clustering;
	int sampleSize;
	int clusteringThreshold;

	void doStep(std::vector<CSequence*>& sequences, tree_structure& tree, int previousTop, bool parallel);

	size_t randomSeeds(
		std::vector<CSequence*>& sequences,
		size_t n_seeds,
		int * seed_ids,
		float * similarity_row);

	size_t clusterSeeds(
		std::vector<CSequence*>& sequences,
		size_t n_seeds,
		size_t n_samples,
		int * seed_ids,
		float * similarity_row);
};