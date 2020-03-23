/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#pragma once

#include "TreeDefs.h"

#include "../core/sequence.h"
#include "../core/defs.h"

#include <vector>

class CLCSBP;

enum class Transformation {
	None,
	Reciprocal,
	Inverse,
};


class AbstractTreeGeneator {
public:

	AbstractTreeGeneator(double indel_exp, size_t n_threads);
	
	void operator()(std::vector<CSequence>& sequences, tree_structure& tree);

	
	template <class seq_type, class similarity_type, Transformation transformation = Transformation::None>
	void calculateSimilarityVector(
		CSequence* ref,
		seq_type* sequences,
		size_t n_seqs,
		similarity_type* out_vector,
		CLCSBP& lcsbp);

	template <class seq_type, class similarity_type, Transformation transformation = Transformation::None>
	void calculateSimilarityMatrix(
		seq_type* sequences,
		size_t n_seq,
		similarity_type* out_matrix,
		CLCSBP& lcsbp);


#ifdef DEVELOPER_MODE
	size_t refSequencesSubTreeSize(
		const std::vector<CSequence>& sequences,
		const std::vector<CSequence>& ref_sequences,
		double *monte_carlo_subtree_size);

#endif

protected:
	double indel_exp;
	size_t n_threads;
	instruction_set_t instruction_set;
	
	virtual void run(std::vector<CSequence>& sequences, tree_structure& tree) = 0;

#ifdef DEVELOPER_MODE
	size_t subTreeSize(
		const std::vector<CSequence>& sequences,
		const std::vector<CSequence>& ref_sequences,
		const set<int> &seq_ids);
#endif

};
