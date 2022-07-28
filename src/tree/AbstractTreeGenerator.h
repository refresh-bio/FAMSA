/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#pragma once

#include "TreeDefs.h"

#include "../core/sequence.h"
#include "../core/defs.h"

#include <vector>

class CLCSBP;
#include "../lcs/lcsbp.h"


class AbstractTreeGenerator {
public:

	AbstractTreeGenerator(int n_threads, instruction_set_t instruction_set);

	virtual ~AbstractTreeGenerator() {}
	
	void operator()(std::vector<CSequence*>& sequences, tree_structure& tree);
		
	template <class seq_type, class distance_type, typename Transform>
	void calculateDistanceVector(
		Transform& transform,
		seq_type& ref,
		seq_type* sequences,
		int n_seqs,
		distance_type* out_vector,
		CLCSBP& lcsbp);

	template <class seq_type, class distance_type, typename Iter, typename Transform>
	void calculateDistanceRange(
		Transform& transform, 
		seq_type& ref,
		seq_type* sequences,
		pair<Iter, Iter> ids_range,
		distance_type* out_vector,
		CLCSBP& lcsbp);

	template <class seq_type, class sv_type, class distance_type, typename Iter, typename Transform>
	void calculateDistanceRangeSV(
		Transform& transform, 
		seq_type& ref,
		sv_type* sv,
		pair<Iter, Iter> ids_range,
		distance_type* out_vector,
		CLCSBP& lcsbp);

	template <class seq_type, class distance_type, typename Transform>
	void calculateDistanceMatrix(
		Transform& transform,
		seq_type* sequences,
		int n_seq,
		distance_type* out_matrix,
		CLCSBP& lcsbp);


#ifdef DEVELOPER_MODE
	size_t refSequencesSubTreeSize(
		const std::vector<CSequence>& sequences,
		const std::vector<CSequence>& ref_sequences,
		double *monte_carlo_subtree_size);

#endif

protected:
	int n_threads;
	instruction_set_t instruction_set;
	
	virtual void run(std::vector<CSequence*>& sequences, tree_structure& tree) = 0;

#ifdef DEVELOPER_MODE
	size_t subTreeSize(
		const std::vector<CSequence>& sequences,
		const std::vector<CSequence>& ref_sequences,
		const set<int> &seq_ids);
#endif

};
