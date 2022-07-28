/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

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

// *******************************************************************
// Queue for UPGMA
// *******************************************************************
class CUPGMAQueue
{
	std::vector<CSequence*> *sequences;
	uint32_t n_rows;
	UPGMA_dist_t *dist_matrix;
	uint32_t lowest_uncomputed_row;
	bool eoq_flag;

	std::mutex mtx;
	

public:
	CUPGMAQueue(std::vector<CSequence*> *_sequences, uint32_t _n_rows, UPGMA_dist_t *_dist_matrix) : 
		sequences(_sequences), n_rows(_n_rows), dist_matrix(_dist_matrix), 
		lowest_uncomputed_row(0), eoq_flag(false) 
	{}

	~CUPGMAQueue() {}

	bool GetTask(int& row_id, std::vector<CSequence*>*& _sequences, UPGMA_dist_t*& dist_row) {
		unique_lock<mutex> lck(mtx);

		if (eoq_flag)
			return false;	// End of data in the profiles queue

		row_id = lowest_uncomputed_row++;

		if (lowest_uncomputed_row >= n_rows)
			eoq_flag = true;

		_sequences = sequences;
		dist_row = dist_matrix + TriangleMatrix::access(row_id, 0);

		return true;
	}
};

// *******************************************************************
// UPGMA algo
// *******************************************************************
template <Distance _distance>
class UPGMA : public AbstractTreeGenerator, public IPartialGenerator {
public:
	
	UPGMA(int n_threads, instruction_set_t instruction_set, bool is_modified)
		: AbstractTreeGenerator(n_threads, instruction_set), is_modified(is_modified) {}

	void run(std::vector<CSequence*>& sequences, tree_structure& tree) override;

	void runPartial(std::vector<CSequence*>& sequences, tree_structure& tree) override;

	void computeDistances(std::vector<CSequence*>& sequences, UPGMA_dist_t *dist_matrix);

	template <bool is_modified>
	void computeTree(UPGMA_dist_t* distances, int n_seq, tree_structure& tree);

protected:
	const UPGMA_dist_t BIG_DIST = (UPGMA_dist_t) 1e29;
	bool is_modified;
};