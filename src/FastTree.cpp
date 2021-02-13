/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#include "FastTree.h"
#include "AbstractTreeGenerator.hpp"
#include "lcsbp.h"
#include "deterministic_random.h"
#include "log.h"

#include <algorithm>
#include <random>
#include <numeric>

#define DIST_MEASURE Measure::DistanceReciprocal

// *******************************************************************
FastTree::FastTree(
	double indel_exp,
	size_t n_threads,
	std::shared_ptr<IPartialGenerator> partialGenerator,
	int subtreeSize,
	std::shared_ptr<IClustering> clustering,
	int sampleSize)
	:
		AbstractTreeGenerator(indel_exp, n_threads),
		partialGenerator(partialGenerator),
		subtreeSize(subtreeSize),
		clustering(clustering),
		sampleSize(sampleSize),
		clusteringThreshold(3 * subtreeSize)
{}


// *******************************************************************
void FastTree::run(std::vector<CSequence>& sequences, tree_structure& tree)
{
	// create vector of pointers to be passed to the recursion
	std::vector<CSequence*> sequencePtrs(sequences.size());
	std::transform(sequences.begin(), sequences.end(), sequencePtrs.begin(), [](CSequence& s)->CSequence* { return &s; });

	randomIds.resize(sequences.size());

	doStep(sequencePtrs, tree);
}


// *******************************************************************
void FastTree::doStep(std::vector<CSequence*>& sequences, tree_structure& tree)
{
	size_t n_seqs = sequences.size();
	CLCSBP lcsbp(instruction_set);

	if ((!clustering && n_seqs > subtreeSize) || (clustering && n_seqs > clusteringThreshold)) {

		float* dists = new float[n_seqs * 2]; // second row will be used later
		int* seed_ids = new int[subtreeSize];

		float* dist_row = dists;
		size_t n_seeds;

		if (clustering == nullptr) {
			n_seeds = randomSeeds(sequences, subtreeSize, seed_ids, dist_row);
		}
		else {
			n_seeds = clusterSeeds(sequences, subtreeSize, sampleSize, seed_ids, dist_row);
		}

		//
		// --- Clustering ---
		//

		// assume all sequences are clustered to 0'th seed at the beginning
		std::vector<CSequence*> seeds(n_seeds);
		int* assignments = new int[n_seqs];
		std::fill_n(assignments, n_seqs, 0);

		// make assignments of all sequences to seeds
		float* current_row = dists + n_seqs;
		seeds[0] = sequences[seed_ids[0]];
		for (int k = 1; k < seeds.size(); ++k) {
			seeds[k] = sequences[seed_ids[k]];
			calculateSimilarityVector<CSequence*, float, DIST_MEASURE>(seeds[k], sequences.data(), n_seqs, current_row, lcsbp);

			for (size_t j = 0; j < n_seqs; ++j) {
				if (current_row[j] < dist_row[j]) { 	// use dist_row for storing smallest distances
					dist_row[j] = current_row[j];
					assignments[j] = k;
				}
			}
		}

		// get histogram
		int* histogram = new int[seeds.size()];
		std::fill_n(histogram, seeds.size(), 0);
		for (size_t j = 0; j < n_seqs; ++j) {
			++histogram[assignments[j]];
		}

		// reserve memory for subgroups
		std::vector<std::vector<CSequence*>> subgroups(seeds.size());
		for (int k = 0; k < seeds.size(); ++k) {
			subgroups[k].reserve(histogram[k]);
			assignments[seed_ids[k]] = k; // mark seeds as assigned to themselves
		}

		// add sequences to subgroups
		for (size_t j = 0; j < n_seqs; ++j) {
			subgroups[assignments[j]].push_back(sequences[j]);
		}

		// delete everything before proceeding with the recursion
		delete[] histogram;
		delete[] assignments;
		delete[] seed_ids;
		delete[] dists;

		// process child nodes
		std::vector<int> subroots(seeds.size(), -1);
		for (int k = 0; k < seeds.size(); ++k) {
			auto &subgroup = subgroups[k];

			// create a subtree when more than 1 element
			if (subgroup.size() > 1) {
				doStep(subgroup, tree);
				subroots[k] = tree.size() - 1;
			}
		}

		size_t previousTop = tree.size();

		partialGenerator->runPartial(seeds, tree);

		// correct node identifiers in the guide tree
		for (size_t node_id = previousTop; node_id < tree.size(); ++node_id) {
			auto& node = tree[node_id];

			node.first = node.first < seeds.size()
				? (subgroups[node.first].size() > 1 ? subroots[node.first] : seeds[node.first]->sequence_no)	// case: seed id - change to subroot or seq id
				: node.first + previousTop - seeds.size();	 // case: intermediate node

			node.second = node.second < seeds.size()
				? (subgroups[node.second].size() > 1 ? subroots[node.second] : seeds[node.second]->sequence_no)	// case: seed id - change to subroot or seq id						 // case: seed id - change to subroot
				: node.second + previousTop - seeds.size();	 // case: intermediate node
		}
	}
	else {

		size_t previousTop = tree.size();

		partialGenerator->runPartial(sequences, tree);

		// correct node identifiers in the guide tree
		if (previousTop > sequences.size()) {
			for (size_t node_id = previousTop; node_id < tree.size(); ++node_id) {
				auto& node = tree[node_id];

				node.first = node.first < sequences.size()
					? sequences[node.first]->sequence_no				// case: sequence id
					: node.first + previousTop - sequences.size();	// case: intermediate node

				node.second = node.second < sequences.size()
					? sequences[node.second]->sequence_no				// case: sequence id
					: node.second + previousTop - sequences.size();	// case: intermediate node
			}
		}
	}
}


// *******************************************************************
size_t FastTree::randomSeeds(
	std::vector<CSequence*>& sequences,
	size_t n_seeds,
	int * seed_ids,
	float * dist_row)
{
	CLCSBP lcsbp(instruction_set);
	size_t n_seqs = sequences.size();

	// calculate distances 0'th (longest) vs all (first row)
	calculateSimilarityVector<CSequence*, float, DIST_MEASURE>(sequences[0], sequences.data(), n_seqs, dist_row, lcsbp);

	std::mt19937 mt;
	std::iota(randomIds.begin(), randomIds.begin() + n_seqs, 0);

	size_t furthestId = std::max_element(dist_row + 1, dist_row + n_seqs) - dist_row;
	std::swap(randomIds[1], randomIds[furthestId]);
	partial_shuffle(randomIds.begin() + 2, randomIds.begin() + n_seeds, randomIds.begin() + n_seqs, mt);

	std::copy(randomIds.begin(), randomIds.begin() + n_seeds, seed_ids);
	std::sort(seed_ids, seed_ids + n_seeds);

	return n_seeds;
}

// *******************************************************************
size_t FastTree::clusterSeeds(
	std::vector<CSequence*>& sequences,
	size_t n_seeds,
	size_t n_samples,
	int * seed_ids,
	float * dist_row)
{
	CLCSBP lcsbp(instruction_set);
	size_t n_seqs = sequences.size();

	CSequence** samples = nullptr;
	int * sample_ids = nullptr;

	// calculate distances 0'th (longest) vs all (first row)
	calculateSimilarityVector<CSequence*, float, DIST_MEASURE>(sequences[0], sequences.data(), n_seqs, dist_row, lcsbp);

	if (n_samples >= sequences.size()) {
		// use all sequences as a sample - don't need to sample anything
		n_samples = sequences.size();
		samples = sequences.data();
	}
	else {

		std::mt19937 mt;
		std::iota(randomIds.begin(), randomIds.begin() + n_seqs, 0);
		partial_shuffle(randomIds.begin() + 1, randomIds.begin() + n_samples, randomIds.begin() + n_seqs, mt);

		sample_ids = new int[n_samples];
		std::copy(randomIds.begin(), randomIds.begin() + n_samples, sample_ids);
		std::sort(sample_ids, sample_ids + n_samples);

		samples = new CSequence*[n_samples];
		for (size_t j = 0; j < n_samples; ++j) {
			samples[j] = sequences[sample_ids[j]];
		}
	}

	//LOG_DEBUG << "Seqs count: " << sequences.size() << ", samples count: " << n_samples << std::endl;


	// calculate distance matrix
	float* distances = TriangleMatrix::allocate<float>(n_samples);
	this->calculateSimilarityMatrix<CSequence*, float, DIST_MEASURE>(samples, n_samples, distances, lcsbp);

	// perform clustering
	(*clustering)(distances, n_samples, n_seeds, 1, seed_ids);

	// translate ids from sample space to sequence space
	if (sample_ids != nullptr) {
		for (size_t k = 0; k < n_seeds; ++k) {
			seed_ids[k] = sample_ids[seed_ids[k]];
		}

		//std::sort(seed_ids, seed_ids + n_seeds);

		delete[] samples;
	}

	// free memory
	delete[] distances;
	delete[] sample_ids;

	return n_seeds;
}
