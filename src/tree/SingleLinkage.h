/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#pragma once

#include "AbstractTreeGenerator.h"
#include "IPartialGenerator.h"

#include <vector>
#include <stack>
#include <utility>
#include <mutex>
#include <condition_variable>


class CSingleLinkageQueue
{
	std::vector<CSequence> *sequences;

	std::vector<std::vector<double>> sim_vector_2d;
	//	vector<vector<uint32_t>> lcs_len_2d;

	std::vector<std::pair<int, bool>> ready_rows;
	std::stack<int, std::vector<int>> available_buffers;

	uint32_t lowest_uncomputed_row;
	uint32_t n_rows;
	uint32_t max_buffered_rows;

	bool eoq_flag;

	std::mutex mtx;
	std::condition_variable cv_tasks, cv_rows;

public:
	CSingleLinkageQueue(vector<CSequence> *_sequences, uint32_t _n_rows, uint32_t _max_buffered_rows);
	~CSingleLinkageQueue();

	bool GetTask(int &row_id, vector<CSequence> *&_sequences, vector<double> *&sim_vector);
	void RegisterSolution(int row_id);
	bool GetSolution(int row_id, vector<double> *&sim_vector);
	void ReleaseSolution(int row_id);
};



class SingleLinkage : public AbstractTreeGenerator, public IPartialGenerator {
public:

	SingleLinkage(double indel_exp, size_t n_threads) : AbstractTreeGenerator(indel_exp, n_threads) {}

	void run(std::vector<CSequence>& sequences, tree_structure& tree) override;

	void runPartial(std::vector<CSequence*>& sequences, tree_structure& tree) override;
};

