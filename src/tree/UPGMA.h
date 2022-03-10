/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#pragma once

#include "AbstractTreeGenerator.h"
#include "IPartialGenerator.h"

#include <vector>
#include <utility>
#include <condition_variable>
#include <stack>
#include <mutex>

// UPGMA defines and consts
typedef float UPGMA_dist_t;

class CUPGMAQueue
{
	std::vector<CSequence> *sequences;

	uint32_t lowest_uncomputed_row;
	uint32_t n_rows;
	UPGMA_dist_t *dist_matrix;
	bool eoq_flag;
	std::mutex mtx;
	

public:
	CUPGMAQueue(std::vector<CSequence> *_sequences, uint32_t _n_rows, UPGMA_dist_t *_dist_matrix);
	~CUPGMAQueue();

	bool GetTask(int &row_id, std::vector<CSequence> *&_sequences, UPGMA_dist_t *&dist_row);
};

class UPGMA : public AbstractTreeGenerator, public IPartialGenerator {
public:
	
	UPGMA(double indel_exp, size_t n_threads, bool is_modified) 
		: AbstractTreeGenerator(indel_exp, n_threads), is_modified(is_modified) {}

	void run(std::vector<CSequence>& sequences, tree_structure& tree) override;

	void runPartial(std::vector<CSequence*>& sequences, tree_structure& tree) override;

	void computeDistances(std::vector<CSequence>& sequences, UPGMA_dist_t *dist_matrix);

	template <bool is_modified>
	void computeTree(UPGMA_dist_t* distances, size_t n_seq, tree_structure& tree);

protected:
	const UPGMA_dist_t BIG_DIST = (UPGMA_dist_t) 1e29;
	bool is_modified;
};