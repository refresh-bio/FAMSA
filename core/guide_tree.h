/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#pragma once

#include "../core/guide_tree_defs.h"
#include "../core/defs.h"
#include "../core/params.h"
#include "../core/sequence.h"
#include "../core/lcsbp.h"


#include "../libs/vectorclass.h"


#include <vector>
#include <utility>
#include <thread>
#include <set>

using node_t = std::pair<int, int>;
using tree_structure = std::vector<node_t>;

class GuideTree {
public:
	
	tree_structure & getRaw() { return guide_tree; }

	GuideTree(double indel_exp, size_t n_threads, int seed, int parttree_size, bool are_sorted); 

	void compute(std::vector<CSequence>& sequences, GT_method::Value method);

	bool loadNewick(
		const std::string& file_name,
		std::vector<CSequence>& sequences);

	bool saveNewick(
		const std::string& file_name,
		const std::vector<CSequence>& sequences) const;

	bool saveDistances(
		const std::string& file_name,
		std::vector<CSequence>& sequences);

#ifdef DEVELOPER_MODE
	size_t refSequencesSubTreeSize(
		const std::vector<CSequence>& sequences,
		const std::vector<CSequence>& ref_sequences,
		double *monte_carlo_subtree_size);

#endif
	
protected:

	int parttree_size;
	double indel_exp;
	size_t n_threads;
	int seed;
	bool are_sorted;
	
	instruction_set_t instruction_set;
	tree_structure guide_tree;
	
	void computeSingleLinkage(std::vector<CSequence>& sequences, tree_structure& tree);
	void computeSingleLinkage_serial(std::vector<CSequence*>sequences, tree_structure& tree);

	void computeUPGMA(UPGMA_dist_t* distances, size_t n_seq, tree_structure& tree);
	void computePartTree(std::vector<CSequence>& sequences, GT_method::Value method, tree_structure& tree);

	void partTreeStep(std::vector<CSequence*>& sequences, GT_method::Value method, tree_structure& tree);

	void calculateDistances(std::vector<CSequence>& sequences, UPGMA_dist_t *dist_matrix);
	
	
	template <class seq_type, class similarity_type, bool invert = false>
	void calculateSimilarityVector(
		CSequence* ref, 
		seq_type* sequences, 
		size_t n_seqs, 
		similarity_type* out_vector, 
		CLCSBP& lcsbp);

	template <class seq_type, class similarity_type, bool invert = false>
	void calculateSimilarityMatrix(
		seq_type* sequences, 
		size_t n_seq, 
		similarity_type* out_matrix,
		CLCSBP& lcsbp);


#ifdef DEVELOPER_MODE
	void GuideTree::computeChained(std::vector<CSequence>& sequences);

	size_t subTreeSize(
		const std::vector<CSequence>& sequences,
		const std::vector<CSequence>& ref_sequences,
		const set<int> &seq_ids);
#endif

};


