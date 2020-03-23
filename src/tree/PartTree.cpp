/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#include "PartTree.h"
#include "AbstractTreeGenerator.hpp"
#include "../lcs/lcsbp.h"
#include "../utils/deterministic_random.h"

#include <algorithm>
#include <random>

// *******************************************************************
PartTree::PartTree(
	double indel_exp, 
	size_t n_threads, 
	std::shared_ptr<IPartialGenerator> partialGenerator, 
	int subtreeSize, 
	std::shared_ptr<IClustering> clustering,
	int sampleSize)
	: 
		AbstractTreeGeneator(indel_exp, n_threads),
		partialGenerator(partialGenerator),
		subtreeSize(subtreeSize),
		clustering(clustering),
		sampleSize(sampleSize),
		clusteringThreshold(3 * subtreeSize)

{}


// *******************************************************************
void PartTree::run(std::vector<CSequence>& sequences, tree_structure& tree)
{
	// create vector of pointers to be passed to the recursion
	std::vector<CSequence*> sequencePtrs(sequences.size());
	std::transform(sequences.begin(), sequences.end(), sequencePtrs.begin(), [](CSequence& s)->CSequence* { return &s; });

	doStep(sequencePtrs, tree);
}


// *******************************************************************
void PartTree::doStep(std::vector<CSequence*>& sequences, tree_structure& tree)
{
	size_t n_seqs = sequences.size();
	CLCSBP lcsbp(instruction_set);

	if ((!clustering && n_seqs > subtreeSize) || (clustering && n_seqs > clusteringThreshold)) {
		
		float* similarities = new float[n_seqs * 2]; // second row will be used later
		int* seed_ids = new int[subtreeSize];
		
		float* similarity_row = similarities;
		size_t n_seeds;

		if (clustering == nullptr) {
			n_seeds = randomSeeds(sequences, subtreeSize, seed_ids, similarity_row);
		}
		else {
			n_seeds = clusterSeeds(sequences, subtreeSize, sampleSize, seed_ids, similarity_row);
		}

		//
		// --- Clustering ---
		//

		// assume all sequences are clustered to 0'th seed at the beginning
		std::vector<CSequence*> seeds(n_seeds);
		int* assignments = new int[n_seqs];
		std::fill_n(assignments, n_seqs, 0);

		// make assignments of all sequences to seeds
		float* current_row = similarities + n_seqs;
		seeds[0] = sequences[seed_ids[0]];
		for (int k = 1; k < seeds.size(); ++k) {
			seeds[k] = sequences[seed_ids[k]];
			calculateSimilarityVector<CSequence*, float>(seeds[k], sequences.data(), n_seqs, current_row, lcsbp);

			for (size_t j = 0; j < n_seqs; ++j) {
				if (current_row[j] > similarity_row[j]) { 	// use similarity_row for storing maximum similarities 
					similarity_row[j] = current_row[j];
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
			assignments[seed_ids[k]] = -1; // mark seeds as not assigned to anything
		}

		// add sequences to subgroups
		for (size_t j = 0; j < n_seqs; ++j) {
			if (assignments[j] >= 0) { // do not assign seeds themselves
				subgroups[assignments[j]].push_back(sequences[j]);
			}
		}

		// delete everything before proceeding with the recursion
		delete[] histogram;
		delete[] assignments;
		delete[] seed_ids;
		delete[] similarities;

		// process child nodes
		std::vector<size_t> subroots(seeds.size());
		for (int k = 0; k < seeds.size(); ++k) {
			auto &subgroup = subgroups[k];

			if (subgroup.size() > 1) {
				// create a subtree
				doStep(subgroup, tree);
				// make an intermediate node to join a subtree with a seed
				tree.push_back(node_t(seeds[k]->sequence_no, tree.size() - 1));
				subroots[k] = tree.size() - 1;
			}
			else if (subgroup.size() == 1) {
				// merge seed with the only element in the group
				tree.push_back(node_t(seeds[k]->sequence_no, subgroup.front()->sequence_no));
				subroots[k] = tree.size() - 1;
			}
		}

		size_t previousTop = tree.size();

		partialGenerator->runPartial(seeds, tree);

		// correct node identifiers in the guide tree
		for (size_t node_id = previousTop; node_id < tree.size(); ++node_id) {
			auto& node = tree[node_id];

			node.first = node.first < seeds.size()
				? (subgroups[node.first].size() ? subroots[node.first] : seeds[node.first]->sequence_no)	// case: seed id - change to subroot or seq id 
				: node.first + previousTop - seeds.size();	 // case: intermediate node

			node.second = node.second < seeds.size()
				? (subgroups[node.second].size() ? subroots[node.second] : seeds[node.second]->sequence_no)	// case: seed id - change to subroot or seq id						 // case: seed id - change to subroot
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
size_t PartTree::randomSeeds(
	std::vector<CSequence*>& sequences,
	size_t n_seeds,
	int * seed_ids,
	float * similarity_row)
{
	CLCSBP lcsbp(instruction_set);
	size_t n_seqs = sequences.size();

	// calculate distances 0'th (longest) vs all (first row)
	calculateSimilarityVector<CSequence*, float>(sequences[0], sequences.data(), n_seqs, similarity_row, lcsbp);

	// two first seeds: the longest (0'th), the furthest 
	seed_ids[0] = 0;
	seed_ids[1] = std::min_element(similarity_row + 1, similarity_row + n_seqs) - similarity_row;

	// select randomly other seeds 
	det_uniform_int_distribution<size_t> dist(1, n_seqs - 1);
	std::mt19937 mt;
	std::generate(seed_ids + 2, seed_ids + n_seeds, [&dist, &mt]()->size_t { return dist(mt); });
	std::stable_sort(seed_ids, seed_ids + n_seeds);

	return std::unique(seed_ids, seed_ids + n_seeds) - seed_ids;
}

// *******************************************************************
size_t PartTree::clusterSeeds(
	std::vector<CSequence*>& sequences,
	size_t n_seeds,
	size_t n_samples,
	int * seed_ids,
	float * similarity_row)
{
	CLCSBP lcsbp(instruction_set);
	size_t n_seqs = sequences.size();

	// the longest sequences wll be used as a first seed
	seed_ids[0] = 0;
	++seed_ids;
	--n_seeds;

	// calculate distances 0'th (longest) vs all (first row)
	calculateSimilarityVector<CSequence*, float>(sequences[0], sequences.data(), n_seqs, similarity_row, lcsbp);
	
	// select randomly sample
	det_uniform_int_distribution<size_t> dist(1, n_seqs - 1);
	std::mt19937 mt;

	int * sample_ids = new int[n_samples];
	std::generate(sample_ids, sample_ids + n_samples, [&dist, &mt]()->size_t { return dist(mt); });
	std::sort(sample_ids, sample_ids + n_samples);
	n_samples = std::unique(sample_ids, sample_ids + n_samples) - sample_ids;

	CSequence** samples = new CSequence*[n_samples];
	for (size_t j = 0; j < n_samples; ++j) {
		samples[j] = sequences[sample_ids[j]];
	}

	// calculate distance matrix
	float* distances = TriangleMatrix::allocate<float>(n_samples);
	this->calculateSimilarityMatrix<CSequence*, float, Transformation::Reciprocal>(samples, n_samples, distances, lcsbp);
	
	// perform clustering
	(*clustering)(distances, n_samples, n_seeds, seed_ids);

	// translate ids from sample space to sequence space
	for (size_t k = 0; k < n_seeds; ++k) {
		seed_ids[k] = sample_ids[seed_ids[k]];
	}

	// restore pointer position
	seed_ids--;
	n_seeds++;
	
	// free memory
	delete[] distances;
	delete[] samples;
	delete[] sample_ids;

	return n_seeds;
}