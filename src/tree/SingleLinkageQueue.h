#pragma once

#include "../core/sequence.h"

#include <vector>
#include <stack>
#include <utility>
#include <mutex>
#include <condition_variable>

template <class distance_type>
class CSingleLinkageQueue
{
	std::vector<CSequence*>* sequences;

	std::vector<std::vector<distance_type>> sim_vector_2d;
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

	// *******************************************************************
	// CSingleLinkageQueue
	// *******************************************************************
	CSingleLinkageQueue(vector<CSequence*>* _sequences, uint32_t _n_rows, uint32_t _max_buffered_rows)
	{
		sequences = _sequences;
		n_rows = _n_rows;
		max_buffered_rows = min(n_rows, _max_buffered_rows);

		sim_vector_2d.resize(max_buffered_rows);
		for (auto& x : sim_vector_2d)
			x.resize(n_rows);

		ready_rows.resize(n_rows, make_pair(-1, false));

		lowest_uncomputed_row = 0;

		for (int i = 0; i < (int) max_buffered_rows; ++i)
			available_buffers.push(i);

		eoq_flag = false;
	}

	// *******************************************************************
	~CSingleLinkageQueue()
	{
	}

	// *******************************************************************
	bool GetTask(int& row_id, vector<CSequence*>*& _sequences, vector<distance_type>*& sim_vector)
	{
		unique_lock<mutex> lck(mtx);
		cv_tasks.wait(lck, [this] {return !this->available_buffers.empty() || this->eoq_flag; });

		if (eoq_flag)
			return false;	// End of data in the profiles queue

		row_id = lowest_uncomputed_row++;

		if (lowest_uncomputed_row >= n_rows)
			eoq_flag = true;

		_sequences = sequences;

		int buffer_row_id = available_buffers.top();
		available_buffers.pop();

		sim_vector = &sim_vector_2d[buffer_row_id];

		ready_rows[row_id].first = buffer_row_id;

#ifdef PRODUCE_LOG
		cerr << "GetTask : " << row_id << "\n";
		cerr << available_buffers.size() << "  " << lowest_uncomputed_row << "\n";
#endif

		return true;
	}

	// *******************************************************************
	void RegisterSolution(int row_id)
	{
		unique_lock<mutex> lck(mtx);

		ready_rows[row_id].second = true;

#ifdef PRODUCE_LOG
		cerr << "Registered : " << row_id << "\n";
		cerr << available_buffers.size() << "  " << lowest_uncomputed_row << "\n";
#endif

		cv_rows.notify_one();
	}

	// *******************************************************************
	bool GetSolution(int row_id, vector<distance_type>*& sim_vector)
	{
		unique_lock<mutex> lck(mtx);
		cv_rows.wait(lck, [this, row_id] {return this->ready_rows[row_id].second; });

		int buffer_row_id = ready_rows[row_id].first;

		sim_vector = &sim_vector_2d[buffer_row_id];

#ifdef PRODUCE_LOG
		cerr << "GetSol : " << row_id << "\n";
		cerr << available_buffers.size() << "  " << lowest_uncomputed_row << "\n";
#endif

		return true;
	}

	// *******************************************************************
	void ReleaseSolution(int row_id)
	{
		unique_lock<mutex> lck(mtx);

		int buffer_row_id = ready_rows[row_id].first;

		available_buffers.push(buffer_row_id);

#ifdef PRODUCE_LOG
		cerr << "Release : " << row_id << "\n";
		cerr << available_buffers.size() << "  " << lowest_uncomputed_row << "\n";
#endif

		cv_tasks.notify_all();
	}
};