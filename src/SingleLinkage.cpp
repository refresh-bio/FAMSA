/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#include "SingleLinkage.h"
#include "AbstractTreeGenerator.hpp"
#include "lcsbp.h"

#include "log.h"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <random>
#include <thread>
#include <xmmintrin.h>

using namespace std;

// *******************************************************************
void SingleLinkage::run(std::vector<CSequence>& sequences, tree_structure& tree) {
	int next;
	int n_seq = sequences.size();

	int prefetch_offset = 64 * 2;
	int prefetch_offset_2nd = 128 * 2;

	vector<int> pi(n_seq + max(prefetch_offset, prefetch_offset_2nd), 0);
	vector<double> lambda(n_seq);
	vector<double> *sim_vector;

	CSingleLinkageQueue slq(&sequences, n_seq, n_threads * 8);
	vector<thread *> workers(n_threads, nullptr);

	mutex mtx;

	// Calculation of similarities is made in working threads
	for (uint32_t i = 0; i < n_threads; ++i)
		workers[i] = new thread([&] {
		CLCSBP lcsbp(instruction_set);
		int row_id;
		vector<CSequence> *sequences;
		vector<double> *sim_vector;

		while (slq.GetTask(row_id, sequences, sim_vector))
		{
			calculateSimilarityVector<CSequence, double, Measure::SimilarityDefault>(
				(*sequences)[row_id],
				sequences->data(),
				row_id,
				sim_vector->data(),
				lcsbp);

			slq.RegisterSolution(row_id);
		}
	});

	// Single linkage algorithm is here
	for (int i = 0; i < n_seq; ++i)
	{
		pi[i] = i;
		lambda[i] = -infty_double;

		if (i % (100) == 0) {
			LOG_DEBUG << "Computing guide tree - " << fixed << setprecision(1)
				<< 100.0 * ((double)i * (i + 1) / 2) / ((double)n_seq * (n_seq + 1) / 2) << "\%    (" << i << " of " << n_seq << ")  \r";
		}

		slq.GetSolution(i, sim_vector);

		auto p_lambda = lambda.begin();
		auto p_sim_vector = (*sim_vector).begin();
		auto p_pi = pi.begin();

		for (int j = 0; j < i; ++j)
		{
			next = pi[j];

#ifdef _MSC_VER					// Visual C++
			_mm_prefetch((const char*)&(*sim_vector)[*(p_pi + prefetch_offset)], 2);
#endif
#ifdef __GNUC__
			//			__builtin_prefetch((&(*sim_vector)[pi[j + prefetch_offset]]), 1, 2);
			__builtin_prefetch((&(*sim_vector)[*(p_pi + prefetch_offset)]), 1, 2);
#endif

			auto &x = (*sim_vector)[next];

			if (isgreater(*p_lambda, *p_sim_vector))
			{
				x = max(x, *p_sim_vector);
			}
			else
			{
				x = max(x, *p_lambda);
				*p_pi = i;
				*p_lambda = *p_sim_vector;
			}

			++p_pi;
			++p_lambda;
			++p_sim_vector;
		}

		slq.ReleaseSolution(i);

		p_pi = pi.begin();
		p_lambda = lambda.begin();
		for (int j = 0; j < i; ++j)
		{
#ifdef _MSC_VER					// Visual C++
			_mm_prefetch((const char*)&lambda[*(p_pi + prefetch_offset_2nd)], 0);
#endif
#ifdef __GNUC__
			__builtin_prefetch((&lambda[*(p_pi + prefetch_offset_2nd)]), 1, 0);
#endif

			next = *p_pi;
			if (isgreaterequal(lambda[next], *p_lambda))
				*p_pi = i;

			++p_pi;
			++p_lambda;
		}
	}

	for (auto p : workers)
	{
		p->join();
		delete p;
	}
	workers.clear();

	LOG_DEBUG << "Computing guide tree - 100.0\%                                        \r";

	vector<int> elements(n_seq - 1);
	for (int i = 0; i < n_seq - 1; ++i)
		elements[i] = i;

#ifdef DEBUG_MODE
	identity /= n_seq * (n_seq - 1) / 2.0;
#endif

	stable_sort(elements.begin(), elements.end(), [&](int x, int y) {
		return lambda[x] > lambda[y];
	});

	vector<int> index(n_seq);
	for (int i = 0; i < n_seq; ++i)
		index[i] = i;

	for (int i = 0; i < n_seq - 1; ++i)
	{
		int j = elements[i];
		next = pi[j];
		tree.push_back(make_pair(index[j], index[next]));
		index[next] = n_seq + i;
	}
}

// *******************************************************************
void SingleLinkage::runPartial(std::vector<CSequence*>& sequences, tree_structure& tree)
{
	int next;
	int n_seq = sequences.size();

	int prefetch_offset = 64 * 2;
	int prefetch_offset_2nd = 128 * 2;

	vector<int> pi(n_seq + max(prefetch_offset, prefetch_offset_2nd), 0);
	vector<double> lambda(n_seq);
	vector<double> sim_vector(n_seq);

	CLCSBP lcsbp(instruction_set);


	// Single linkage algorithm is here
	for (int i = 0; i < n_seq; ++i)
	{
		pi[i] = i;
		lambda[i] = -infty_double;

/*		if (i % (100) == 0) {
			LOG_DEBUG << "Computing guide tree - " << fixed << setprecision(1)
				<< 100.0 * ((double)i * (i + 1) / 2) / ((double)n_seq * (n_seq + 1) / 2) << "\%    (" << i << " of " << n_seq << ")  \r";
		}
*/
		calculateSimilarityVector<CSequence*, double, Measure::SimilarityDefault>(
			sequences[i],
			sequences.data(),
			i,
			sim_vector.data(),
			lcsbp);

		auto p_lambda = lambda.begin();
		auto p_sim_vector = sim_vector.begin();
		auto p_pi = pi.begin();

		for (int j = 0; j < i; ++j)
		{
			next = pi[j];

#ifdef _MSC_VER					// Visual C++
			_mm_prefetch((const char*)&(sim_vector)[*(p_pi + prefetch_offset)], 2);
#endif
#ifdef __GNUC__
			//			__builtin_prefetch((&(*sim_vector)[pi[j + prefetch_offset]]), 1, 2);
			__builtin_prefetch((&(sim_vector)[*(p_pi + prefetch_offset)]), 1, 2);
#endif

			auto &x = (sim_vector)[next];

			if (isgreater(*p_lambda, *p_sim_vector))
			{
				x = max(x, *p_sim_vector);
			}
			else
			{
				x = max(x, *p_lambda);
				*p_pi = i;
				*p_lambda = *p_sim_vector;
			}

			++p_pi;
			++p_lambda;
			++p_sim_vector;
		}

		p_pi = pi.begin();
		p_lambda = lambda.begin();
		for (int j = 0; j < i; ++j)
		{
#ifdef _MSC_VER					// Visual C++
			_mm_prefetch((const char*)&lambda[*(p_pi + prefetch_offset_2nd)], 0);
#endif
#ifdef __GNUC__
			__builtin_prefetch((&lambda[*(p_pi + prefetch_offset_2nd)]), 1, 0);
#endif

			next = *p_pi;
			if (isgreaterequal(lambda[next], *p_lambda))
				*p_pi = i;

			++p_pi;
			++p_lambda;
		}
	}

//	LOG_DEBUG << "Computing guide tree - 100.0\%                                        \r";

	vector<int> elements(n_seq - 1);
	for (int i = 0; i < n_seq - 1; ++i)
		elements[i] = i;

#ifdef DEBUG_MODE
	identity /= n_seq * (n_seq - 1) / 2.0;
#endif

	stable_sort(elements.begin(), elements.end(), [&](int x, int y) {
		return lambda[x] > lambda[y];
	});

	vector<int> index(n_seq);
	for (int i = 0; i < n_seq; ++i)
		index[i] = i;

	for (int i = 0; i < n_seq - 1; ++i)
	{
		int j = elements[i];
		next = pi[j];
		tree.push_back(make_pair(index[j], index[next]));
		index[next] = n_seq + i;
	}
}


// *******************************************************************
// CSingleLinkageQueue
// *******************************************************************
CSingleLinkageQueue::CSingleLinkageQueue(vector<CSequence> *_sequences, uint32_t _n_rows, uint32_t _max_buffered_rows)
{
	sequences = _sequences;
	n_rows = _n_rows;
	max_buffered_rows = min(n_rows, _max_buffered_rows);

	sim_vector_2d.resize(max_buffered_rows);
	for (auto &x : sim_vector_2d)
		x.resize(n_rows);

	ready_rows.resize(n_rows, make_pair(-1, false));

	lowest_uncomputed_row = 0;

	for (int i = 0; i < max_buffered_rows; ++i)
		available_buffers.push(i);

	eoq_flag = false;
}

// *******************************************************************
CSingleLinkageQueue::~CSingleLinkageQueue()
{
}

// *******************************************************************
bool CSingleLinkageQueue::GetTask(int &row_id, vector<CSequence> *&_sequences, vector<double> *&sim_vector)
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
void CSingleLinkageQueue::RegisterSolution(int row_id)
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
bool CSingleLinkageQueue::GetSolution(int row_id, vector<double> *&sim_vector)
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
void CSingleLinkageQueue::ReleaseSolution(int row_id)
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
