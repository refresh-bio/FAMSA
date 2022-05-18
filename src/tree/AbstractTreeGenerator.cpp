/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#include "AbstractTreeGenerator.h"

#if SIMD==SIMD_AVX1 || SIMD==SIMD_AVX2 || SIMD==SIMD_AVX512
#include "../utils/cpuid.h"
#endif

#include <algorithm>

using namespace std;


// *******************************************************************
AbstractTreeGenerator::AbstractTreeGenerator(int n_threads, instruction_set_t instruction_set) 
	: n_threads(n_threads), instruction_set(instruction_set) {

}

// *******************************************************************
void AbstractTreeGenerator::operator()(std::vector<CSequence>& sequences, tree_structure& tree)
{
	tree.clear();
	tree.resize(sequences.size(), std::make_pair<int, int>(-1, -1));

	// Temporarily resize the sequences by adding unknown symbols (just to simplify the pairwise comparison)
	size_t max_seq_len =
		max_element(sequences.begin(), sequences.end(), [](const CSequence &x, CSequence &y) {return x.length < y.length; })->length;

	auto mma = sequences.front().get_mma();

	if(mma)
		mma->freeze();

	int n_seqs = (int)sequences.size();
	for (int i = 0; i < n_seqs; ++i) {
		sequences[i].DataResize(max_seq_len, UNKNOWN_SYMBOL);
	}

	if (mma)
		mma->release_freezed();

	// build the tree
	run(sequences, tree);

	if (mma)
		mma->freeze();

	// Bring the sequences to the valid length
	for (int i = 0; i < n_seqs; ++i)
		sequences[i].DataResize(sequences[i].length, UNKNOWN_SYMBOL);

	if (mma)
		mma->release_freezed();
}


// *******************************************************************
#ifdef DEVELOPER_MODE
size_t AbstractTreeGeneator::refSequencesSubTreeSize(
	const vector<CSequence>& sequences,
	const vector<CSequence>& ref_sequences,
	double *monte_carlo_subtree_size)
{
	const int monte_carlo_trials = 1000;

	set<int> ref_seq_ids;
	int n_seq = sequences.size();
	int r = 0;

	if (ref_sequences.size() == 1)
		return 1;

	// Find the ids of the referential sequences in the input file
	for (int i = 0; i < n_seq; ++i)
	{
		bool is_ref = false;
		for (auto &y : ref_sequences)
			if (sequences[i].id == y.id)
				is_ref = true;

		if (is_ref)
			ref_seq_ids.insert(i);
	}

	r = subTreeSize(sequences, ref_sequences, ref_seq_ids);

	if (monte_carlo_subtree_size)
	{
		mt19937 mt;
		double mc_r = 0;

		for (int i = 0; i < monte_carlo_trials; ++i)
		{
			set<int> mc_seq_ids;

			while (mc_seq_ids.size() < ref_seq_ids.size())
				mc_seq_ids.insert(mt() % n_seq);
			mc_r += subTreeSize(sequences, ref_sequences, mc_seq_ids);
		}

		*monte_carlo_subtree_size = mc_r / (double)monte_carlo_trials;
	}

	return r;
}

#endif
