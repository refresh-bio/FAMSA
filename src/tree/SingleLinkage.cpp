/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

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
#include <array>

#ifdef _MSC_VER	
	#include <immintrin.h>
#endif

using namespace std;

// *******************************************************************
template <Distance _distance>
void SingleLinkage<_distance>::run(std::vector<CSequence*>& sequences, tree_structure& tree) {
	int next;
	int n_seq = (int)sequences.size();

	int prefetch_offset = 64 * 2;
	int prefetch_offset_2nd = 128 * 2;

	vector<int> pi(n_seq + max(prefetch_offset, prefetch_offset_2nd), 0);
	vector<slink_dist_t> lambda(n_seq);
	vector<slink_dist_t> *dist_vector;

	CSingleLinkageQueue<slink_dist_t> slq(&sequences, n_seq, n_threads * 8);
	vector<thread *> workers(n_threads, nullptr);

	mutex mtx;

	// Calculation of similarities is made in working threads
	for (int i = 0; i < n_threads; ++i)
		workers[i] = new thread([&] {
		CLCSBP lcsbp(instruction_set);
		int row_id;
		vector<CSequence*> *sequences;
		vector<slink_dist_t> *dist_vector;
		Transform<double, _distance> transform;

		vector<double> loc_dist_vector;

		while (slq.GetTask(row_id, sequences, dist_vector))
		{
			loc_dist_vector.resize(dist_vector->size());

			calculateDistanceVector<CSequence*, double, decltype(transform)>(
				transform,
				(*sequences)[row_id],
				sequences->data(),
				row_id,
				loc_dist_vector.data(),
				lcsbp);

#ifdef SLINK_HANDLE_TIES
			for (size_t i = 0; i < loc_dist_vector.size(); ++i)
			{
				(*dist_vector)[i].first = loc_dist_vector[i];
				(*dist_vector)[i].second = ids_to_uint64(i, row_id);
			}
#else
			swap(*dist_vector, loc_dist_vector);
#endif

			slq.RegisterSolution(row_id);
		}
	});

	// Single linkage algorithm is here
	for (int i = 0; i < n_seq; ++i)
	{
		pi[i] = i;
#ifdef SLINK_HANDLE_TIES
		lambda[i] = slink_dist_t{ std::numeric_limits<double>::max(), 0 };
#else
		lambda[i] = std::numeric_limits<double>::max();
#endif

		if (i % (100) == 0) {
			LOG_DEBUG << "Computing guide tree - " << fixed << setprecision(1)
				<< 100.0 * ((double)i * (i + 1) / 2) / ((double)n_seq * (n_seq + 1) / 2) << "%    (" << i << " of " << n_seq << ")  \r";
		}

		slq.GetSolution(i, dist_vector);

		auto p_lambda = lambda.begin();
		auto p_dist_vector = (*dist_vector).begin();
		auto p_pi = pi.begin();

		for (int j = 0; j < i; ++j)
		{
			next = pi[j];

#ifdef _MSC_VER					// Visual C++
			_mm_prefetch((const char*)&(*dist_vector)[*(p_pi + prefetch_offset)], 2);
#endif
#ifdef __GNUC__
			//			__builtin_prefetch((&(*dist_vector)[pi[j + prefetch_offset]]), 1, 2);
			__builtin_prefetch((&(*dist_vector)[*(p_pi + prefetch_offset)]), 1, 2);
#endif

			auto &x = (*dist_vector)[next];

			if (*p_lambda < *p_dist_vector)
			{
				x = min(x, *p_dist_vector);
			}
			else
			{
				x = min(x, *p_lambda);
				*p_pi = i;
				*p_lambda = *p_dist_vector;
			}

			++p_pi;
			++p_lambda;
			++p_dist_vector;
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
			if (lambda[next] <= *p_lambda)
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

	LOG_DEBUG << "Computing guide tree - 100.0%                                        \r";

	vector<int> elements(n_seq - 1);
	for (int i = 0; i < n_seq - 1; ++i)
		elements[i] = i;

#ifdef DEBUG_MODE
	identity /= n_seq * (n_seq - 1) / 2.0;
#endif

	stable_sort(elements.begin(), elements.end(), [&](int x, int y) {
		return lambda[x] < lambda[y];
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
template <Distance _distance>
void SingleLinkage<_distance>::runPartial(std::vector<CSequence*>& sequences, tree_structure& tree)
{
	int next;
	int n_seq = sequences.size();

	int prefetch_offset = 64 * 2;
	int prefetch_offset_2nd = 128 * 2;

	vector<int> pi(n_seq + max(prefetch_offset, prefetch_offset_2nd), 0);
	vector<slink_dist_t> lambda(n_seq);
	vector<slink_dist_t> dist_vector(n_seq);
	vector<double> loc_dist_vector(dist_vector.size());
	
	Transform<double, _distance> transform;

	CLCSBP lcsbp(instruction_set);

	// Single linkage algorithm is here
	for (int i = 0; i < n_seq; ++i)
	{
		pi[i] = i;

#ifdef SLINK_HANDLE_TIES
		lambda[i] = slink_dist_t{ std::numeric_limits<double>::max(), 0 };
#else
		lambda[i] = std::numeric_limits<double>::max();
#endif

/*		if (i % (100) == 0) {
			LOG_DEBUG << "Computing guide tree - " << fixed << setprecision(1)
				<< 100.0 * ((double)i * (i + 1) / 2) / ((double)n_seq * (n_seq + 1) / 2) << "\%    (" << i << " of " << n_seq << ")  \r";
		}
*/
		calculateDistanceVector<CSequence*, double, decltype(transform)>(
			transform,
			sequences[i],
			sequences.data(),
			i,
			loc_dist_vector.data(),
			lcsbp);

#ifdef SLINK_HANDLE_TIES
		for (int j = 0; j < (int) loc_dist_vector.size(); ++j)
		{
			dist_vector[j].first = loc_dist_vector[j];
			dist_vector[j].second = ids_to_uint64(j, i);
		}
#else
		swap(dist_vector, loc_dist_vector);
#endif

		auto p_lambda = lambda.begin();
		auto p_dist_vector = dist_vector.begin();
		auto p_pi = pi.begin();

		for (int j = 0; j < i; ++j)
		{
			next = pi[j];

#ifdef _MSC_VER					// Visual C++
			_mm_prefetch((const char*)&(dist_vector)[*(p_pi + prefetch_offset)], 2);
#endif
#ifdef __GNUC__
			//			__builtin_prefetch((&(*dist_vector)[pi[j + prefetch_offset]]), 1, 2);
			__builtin_prefetch((&(dist_vector)[*(p_pi + prefetch_offset)]), 1, 2);
#endif

			auto &x = (dist_vector)[next];

			if (*p_lambda < *p_dist_vector)
			{
				x = min(x, *p_dist_vector);
			}
			else
			{
				x = min(x, *p_lambda);
				*p_pi = i;
				*p_lambda = *p_dist_vector;
			}

			++p_pi;
			++p_lambda;
			++p_dist_vector;
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
			if (lambda[next] <= *p_lambda)
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
		return lambda[x] < lambda[y];
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
// Explicit template specializations for specified distance measures

template class SingleLinkage<Distance::indel_div_lcs>;
template class SingleLinkage<Distance::sqrt_indel_div_lcs>;
