/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#include "SingleLinkage.h"
#include "AbstractTreeGenerator.hpp"
#include "SingleLinkageQueue.h"
#include "../lcs/lcsbp.h"

#include "../utils/log.h"

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
	vector<slink_similarity_t> lambda(n_seq);
	vector<slink_similarity_t> *sim_vector;

	CSingleLinkageQueue<slink_similarity_t> slq(&sequences, n_seq, n_threads * 8);
	vector<thread *> workers(n_threads, nullptr);

	mutex mtx;

	// Calculation of similarities is made in working threads
	for (uint32_t i = 0; i < n_threads; ++i)
		workers[i] = new thread([&] {
		CLCSBP lcsbp(instruction_set);
		int row_id;
		vector<CSequence> *sequences;
		vector<slink_similarity_t> *sim_vector;
		Transform<double, Measure::SimilarityDefault> transform;
//		Transform<double, Measure::LCS_sqrt_indel> transform;

		vector<double> loc_sim_vector;

		while (slq.GetTask(row_id, sequences, sim_vector))
		{
			loc_sim_vector.resize(sim_vector->size());

			calculateSimilarityVector<CSequence, double, decltype(transform)>(
				transform,
				(*sequences)[row_id],
				sequences->data(),
				row_id,
				loc_sim_vector.data(),
				lcsbp);

#ifdef SLINK_HANDLE_TIES
			for (size_t i = 0; i < loc_sim_vector.size(); ++i)
			{
				(*sim_vector)[i].first = loc_sim_vector[i];
				(*sim_vector)[i].second = ids_to_uint64(i, row_id);
			}
#else
			swap(*sim_vector, loc_sim_vector);
#endif

			slq.RegisterSolution(row_id);
		}
	});

	// Single linkage algorithm is here
	for (int i = 0; i < n_seq; ++i)
	{
		pi[i] = i;
#ifdef SLINK_HANDLE_TIES
		lambda[i] = make_pair(-infty_double, 0);
#else
		lambda[i] = -infty_double;
#endif

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

//			if (isgreater(*p_lambda, *p_sim_vector))
			if (*p_lambda > *p_sim_vector)
//			if (*p_lambda >= *p_sim_vector)
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
//			if (isgreaterequal(lambda[next], *p_lambda))
			if (lambda[next] >= *p_lambda)
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
		tree.emplace_back(index[j], index[next]);
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
	vector<slink_similarity_t> lambda(n_seq);
	vector<slink_similarity_t> sim_vector(n_seq);
	vector<double> loc_sim_vector(sim_vector.size());
	
	Transform<double, Measure::SimilarityDefault> transform;
//	Transform<double, Measure::LCS_sqrt_indel> transform;

	CLCSBP lcsbp(instruction_set);

	// Single linkage algorithm is here
	for (int i = 0; i < n_seq; ++i)
	{
		pi[i] = i;

#ifdef SLINK_HANDLE_TIES
		lambda[i] = make_pair(-infty_double, 0);
#else
		lambda[i] = -infty_double;
#endif

/*		if (i % (100) == 0) {
			LOG_DEBUG << "Computing guide tree - " << fixed << setprecision(1)
				<< 100.0 * ((double)i * (i + 1) / 2) / ((double)n_seq * (n_seq + 1) / 2) << "\%    (" << i << " of " << n_seq << ")  \r";
		}
*/
		calculateSimilarityVector<CSequence*, double, decltype(transform)>(
			transform,
			sequences[i],
			sequences.data(),
			i,
			loc_sim_vector.data(),
			lcsbp);

#ifdef SLINK_HANDLE_TIES
		for (size_t j = 0; j < loc_sim_vector.size(); ++j)
		{
			sim_vector[j].first = loc_sim_vector[j];
			sim_vector[j].second = ids_to_uint64(j, i);
		}
#else
		swap(*sim_vector, loc_sim_vector);
#endif

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

//			if (isgreater(*p_lambda, *p_sim_vector))
			if (*p_lambda > *p_sim_vector)
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
//			if (isgreaterequal(lambda[next], *p_lambda))
			if (lambda[next] >= *p_lambda)
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
		tree.emplace_back(index[j], index[next]);
		index[next] = n_seq + i;
	}
}
