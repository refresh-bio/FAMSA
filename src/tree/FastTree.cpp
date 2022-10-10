/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#include "FastTree.h"
#include "AbstractTreeGenerator.hpp"
#include "../lcs/lcsbp.h"
#include "../utils/deterministic_random.h"
#include "../utils/log.h"
#include "../core/queues.h"

#include <algorithm>
#include <random>
#include <numeric>
#include <thread>
#include <chrono>

// *******************************************************************
template <Distance _distance>
FastTree<_distance>::FastTree(
	int n_threads,
	instruction_set_t instruction_set,
	std::shared_ptr<IPartialGenerator> partialGenerator, 
	int subtreeSize, 
	std::shared_ptr<IClustering> clustering,
	int sampleSize)
	: 
		AbstractTreeGenerator(n_threads, instruction_set),
		partialGenerator(partialGenerator),
		subtreeSize(subtreeSize),
		clustering(clustering),
		sampleSize(sampleSize),
		clusteringThreshold(3 * subtreeSize)
{}


// *******************************************************************
template <Distance _distance>
void FastTree<_distance>::run(std::vector<CSequence*>& sequences, tree_structure& tree)
{
	tree_structure local_tree;
	doStep(sequences, local_tree, tree.size(), true);
	tree.insert(tree.end(), local_tree.begin(), local_tree.end());
	
}


// *******************************************************************
template <Distance _distance>
void FastTree<_distance>::doStep(std::vector<CSequence*>& sequences, tree_structure& tree, int previousTop, bool parallel)
{
	int n_seqs = (int)sequences.size();
	CLCSBP lcsbp(instruction_set);
	Transform<float, _distance> transform;

	if ((!clustering && n_seqs > subtreeSize) || (clustering && n_seqs > clusteringThreshold)) {

		float* dists = new float[sequences.size() * 2]; // second row will be used later
		int* seed_ids = new int[subtreeSize];

		float* dist_row = dists;
		int n_seeds;

		if (clustering == nullptr) {
			n_seeds = randomSeeds(sequences, subtreeSize, seed_ids, dist_row);
		}
		else {
			n_seeds = clusterSeeds(sequences, subtreeSize, sampleSize, seed_ids, dist_row);
		}

		//
		// Clustering 
		//
		// assume all sequences are clustered to 0'th seed at the beginning
		std::vector<CSequence*> seeds(n_seeds);
		int* assignments = new int[n_seqs];
		std::fill_n(assignments, n_seqs, 0);

		// make assignments of all sequences to seeds
		float* current_row = dists + n_seqs;
		seeds[0] = sequences[seed_ids[0]];
		for (int k = 1; k < n_seeds; ++k) {
			seeds[k] = sequences[seed_ids[k]];
			calculateDistanceVector<CSequence*, float, decltype(transform)>(transform, seeds[k], sequences.data(), n_seqs, current_row, lcsbp);

			for (int j = 0; j < n_seqs; ++j) {
				if (current_row[j] < dist_row[j]) { 	// use dist_row for storing smallest distances 
					dist_row[j] = current_row[j];
					assignments[j] = k;
				}
			}
		}

		// get histogram
		int* histogram = new int[seeds.size()];
		std::fill_n(histogram, seeds.size(), 0);
		for (int j = 0; j < n_seqs; ++j) {
			++histogram[assignments[j]];
		}

		// reserve memory for subgroups
		std::vector<std::vector<CSequence*>> subgroups(seeds.size());
		for (int k = 0; k < n_seeds; ++k) {
			subgroups[k].reserve(histogram[k]);
			assignments[seed_ids[k]] = k; // mark seeds as assigned to themselves
		}

		// add sequences to subgroups
		for (int j = 0; j < n_seqs; ++j) {
			subgroups[assignments[j]].push_back(sequences[j]);
		}

		// delete everything before proceeding with the recursion
		delete[] histogram;
		delete[] assignments;
		delete[] seed_ids;
		delete[] dists;

		// process child nodes
		std::vector<int> subroots(seeds.size(), -1);

		if (parallel) {
			//
			// Parallel subtree processing
			//
			std::vector<std::thread> workers(n_threads);
			std::vector<tree_structure> local_trees;
			local_trees.reserve(seeds.size()); // to avoid reallocations

			struct SubtreeTask {
				std::vector<CSequence*>* subgroup;
				tree_structure* localTree;
				int previousTop;
			};
			RegisteringQueue<SubtreeTask> queue(1);

			// schedule tasks
			for (int k = 0; k < n_seeds; ++k) {
				auto& subgroup = subgroups[k];

				// create a subtree when more than 1 element
				if (subgroup.size() > 1) {
					local_trees.push_back(tree_structure());
					SubtreeTask task{ &subgroup, &local_trees.back(), previousTop };
					queue.Push(task);

					previousTop += (int) subgroup.size() - 1;
					subroots[k] = previousTop - 1;
				}
			}

			queue.MarkCompleted();

			// process tasks in threads
			for (auto& w : workers) {
				w = std::thread([this, &queue] {
					SubtreeTask task;
					while (!queue.IsEmpty()) {
						if (queue.Pop(task)) {
							this->doStep(*task.subgroup, *task.localTree, task.previousTop, false);
						}
					}
				});
			}

			// wait for workers to finish
			for (auto& w : workers) {
				w.join();
			}

			// gather partial trees
			for (const auto& lt : local_trees) {
				tree.insert(tree.end(), lt.begin(), lt.end());
			}
		}
		else {
			//
			// Serial subtree processing
			//
			for (int k = 0; k < n_seeds; ++k) {
				auto& subgroup = subgroups[k];

				// create a subtree when more than 1 element
				if (subgroup.size() > 1) {
					tree_structure local_tree;
					doStep(subgroup, local_tree, previousTop, false);
					tree.insert(tree.end(), local_tree.begin(), local_tree.end());
					previousTop += (int) subgroup.size() - 1;
					subroots[k] = previousTop - 1;
				}
			}
		}

		//size_t previousTop = tree.size();
		tree_structure local_tree;
		partialGenerator->runPartial(seeds, local_tree);

		// correct node identifiers in the guide tree
		for (int node_id = 0; node_id < n_seeds - 1; ++node_id) {
			auto& node = local_tree[node_id];

			node.first = node.first < n_seeds
				? (subgroups[node.first].size() > 1 ? subroots[node.first] : seeds[node.first]->sequence_no)	// case: seed id - change to subroot or seq id 
				: node.first + previousTop - (int) n_seeds;	 // case: intermediate node

			node.second = node.second < n_seeds
				? (subgroups[node.second].size() > 1 ? subroots[node.second] : seeds[node.second]->sequence_no)	// case: seed id - change to subroot or seq id						 // case: seed id - change to subroot
				: node.second + previousTop - (int) n_seeds;	 // case: intermediate node
		}
		tree.insert(tree.end(), local_tree.begin(), local_tree.end());

	}
	else {

		//size_t previousTop = tree.size();

		partialGenerator->runPartial(sequences, tree);

		// correct node identifiers in the guide tree
		if (previousTop > n_seqs) {
			for (int node_id = 0; node_id < n_seqs - 1; ++node_id) {
				auto& node = tree[node_id];

				node.first = node.first < n_seqs
					? sequences[node.first]->sequence_no	// case: sequence id 
					: node.first + previousTop - n_seqs;	// case: intermediate node

				node.second = node.second < n_seqs
					? sequences[node.second]->sequence_no	 // case: sequence id 
					: node.second + previousTop - n_seqs;	// case: intermediate node
			}
		}
	}

	//LOG_VERBOSE << "Finished subtree of size: " << sequences.size() << endl;
}


// *******************************************************************
template <Distance _distance>
int FastTree<_distance>::randomSeeds(
	std::vector<CSequence*>& sequences,
	int n_seeds,
	int * seed_ids,
	float * dist_row)
{
	CLCSBP lcsbp(instruction_set);
	size_t n_seqs = sequences.size();

	Transform<float, _distance> transform;

	// calculate distances 0'th (longest) vs all (first row)
	calculateDistanceVector<CSequence*, float, decltype(transform)>(transform, sequences[0], sequences.data(), (int) n_seqs, dist_row, lcsbp);

	std::mt19937 mt;
	int* randomIds = new int[n_seqs];
	std::iota(randomIds, randomIds + n_seqs, 0);
	
	size_t furthestId = std::max_element(dist_row + 1, dist_row + n_seqs) - dist_row;
	std::swap(randomIds[1], randomIds[furthestId]);
	partial_shuffle(randomIds + 2, randomIds + n_seeds, randomIds + n_seqs, mt);

	std::copy(randomIds, randomIds + n_seeds, seed_ids);
	std::sort(seed_ids, seed_ids + n_seeds);
	delete[] randomIds;

	return n_seeds;
}

// *******************************************************************
template <Distance _distance>
int FastTree<_distance>::clusterSeeds(
	std::vector<CSequence*>& sequences,
	int n_seeds,
	int n_samples,
	int * seed_ids,
	float * dist_row)
{
	CLCSBP lcsbp(instruction_set);
	int n_seqs = (int)sequences.size();

	CSequence** samples = nullptr;
	int *sample_ids = nullptr;

	Transform<float, _distance> transform;

	// calculate distances 0'th (longest) vs all (first row)
	calculateDistanceVector<CSequence*, float, decltype(transform)>(transform, sequences[0], sequences.data(), (int) n_seqs, dist_row, lcsbp);

	if (n_samples >= n_seqs) {
		// use all sequences as a sample - don't need to sample anything
		n_samples = n_seqs;
		samples = sequences.data();
	}
	else {

		std::mt19937 mt;
		int* randomIds = new int[n_seqs];
		std::iota(randomIds, randomIds + n_seqs, 0);
		partial_shuffle(randomIds + 1, randomIds + n_samples, randomIds + n_seqs, mt);

		sample_ids = new int[n_samples];
		std::copy(randomIds, randomIds + n_samples, sample_ids);
		std::sort(sample_ids, sample_ids + n_samples);
	
		samples = new CSequence*[n_samples];
		for (int j = 0; j < n_samples; ++j) {
			samples[j] = sequences[sample_ids[j]];
		}
		delete[] randomIds;
	}

	//LOG_DEBUG << "Seqs count: " << sequences.size() << ", samples count: " << n_samples << std::endl;


	// calculate distance matrix
	float* distances = TriangleMatrix::allocate<float>(n_samples);
	this->calculateDistanceMatrix<CSequence*, float, decltype(transform)>(transform, samples, (int) n_samples, distances, lcsbp);
	
	// perform clustering
	(*clustering)(distances, n_samples, n_seeds, 1, seed_ids);

	// translate ids from sample space to sequence space
	if (sample_ids != nullptr) {
		for (int k = 0; k < n_seeds; ++k) {
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

// *******************************************************************
// Explicit template specializations for specified distance measures

template class FastTree<Distance::indel_div_lcs>;
template class FastTree<Distance::sqrt_indel_div_lcs>;

