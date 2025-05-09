/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#include "MSTPrim.h"
#include "AbstractTreeGenerator.hpp"
#include "SingleLinkage.h"

#include "../utils/log.h"
#include "../core/queues.h"

#include "../libs/refresh/active_thread_pool/lib/utils.h"

#include <algorithm>
#include <atomic>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <random>
#include <thread>
#include <utility>
#include <limits>
#include <functional>

#include <chrono>

using namespace std;

#define USE_THEORETICAL_BEST_POSSIBLE_DIST

#ifndef USE_THEORETICAL_BEST_POSSIBLE_DIST	
//#pragma GCC optimize("align-loops=16")
// *******************************************************************
template <Distance _distance>
void MSTPrim<_distance>::run_view(std::vector<CSequence*>& sequences, tree_structure& tree)
{
	run_view_new(sequences, tree);
	return;

	int n_seq = (int)sequences.size();
	int c_seq = 0;

	CSequence ref_seq ("", "");		// Object just to simplify implementation of variants of CLCSBP

	prepare_sequences_view(sequences);

	v_distances.resize(n_seq, dist_t{ std::numeric_limits<dist_value_t>::max(), 0 });

	vector<mst_edge_t> mst_edges;
	vector<int> v_prim_orders(n_seq, n_seq);

	vector<int> candidates(n_threads, -1);

	vector<thread> workers;

//	MSTPartitioner mst_partitioner(n_threads, n_threads * 16, 8, n_threads * 2);
	MSTPartitionerNew mst_partitioner(n_threads, n_threads * 16, 8, n_threads * 2);

	atomic<uint32_t> part_id{ 0 };
	int cur_seq_id = 0;
	//	int cur_seq_id = n_seq-1;
	int cur_prim_order = 0;
	int cur_no_parts = 0;

	v_prim_orders[cur_seq_id] = cur_prim_order++;

	mst_partitioner.InitPartition(n_seq);

	mst_partitioner.Remove(cur_seq_id);
	prepare_bit_masks_for_sequence(*sequences[cur_seq_id], ref_seq.p_bit_masks, ref_seq.p_bv_len);		// the only allocation of ref_bm is here
	ref_seq.length = sequences[cur_seq_id]->length;
	ref_seq.sequence_no = sequences[cur_seq_id]->sequence_no;

	workers.reserve(n_threads);

	atomic<int> a_cnt0{ (int)n_threads - 1 };
	atomic_flag b_flag0;

	cur_no_parts = mst_partitioner.GetNoParts();

//	atomic<uint64_t> global_no_comp_necessary = 0;
//	atomic<uint64_t> global_no_comp_useless = 0;

	b_flag0.test_and_set();

	for (int i = 0; i < n_threads; ++i)
	{
		workers.push_back(thread([&, i] {
			auto id_worker = i;
			CLCSBP lcsbp(instruction_set);
			Transform<double, _distance> transform;

//			uint64_t no_comp_necessary = 0;
//			uint64_t no_comp_useless = 0;

			auto sv = sequence_views;
			auto v_distances_data = v_distances.data();

			vector<dist_value_t> part_dist;

			int& candidate_worker = candidates[id_worker];

			bool flag_value = true;

			vector<int> ids_sub_range;

			int seq_id = cur_seq_id;
			int no_parts = cur_no_parts;

			while (true)
			{
				int loc_part_id = part_id++;

				if (loc_part_id >= no_parts)
				{
					if (seq_id < 0)
						break;

					if (!a_cnt0.fetch_sub(1, memory_order_relaxed))
					{
						a_cnt0 = n_threads - 1;

						// All threads are over, so we can pick a new sequence

						int best_candidate_id = -1;
						for (int j = 0; j < n_threads; ++j)
							if (best_candidate_id == -1)
								best_candidate_id = candidates[j];
							else if (candidates[j] != -1)
							{
								if (v_distances_data[candidates[j]] < v_distances_data[best_candidate_id])
									best_candidate_id = candidates[j];
							}

						fill_n(candidates.begin(), n_threads, -1);

 						auto pids = uint64_to_id(~v_distances[best_candidate_id].second);

						mst_edges.emplace_back(pids.first, pids.second, cur_prim_order, -v_distances[best_candidate_id].first);	// We use MaxRangeQueries so change here to negative value

						if (v_prim_orders[pids.first] == n_seq)
							v_prim_orders[pids.first] = cur_prim_order++;
						else
							v_prim_orders[pids.second] = cur_prim_order++;

						if (!mst_partitioner.IsAlmostEmpty())
						{
							cur_seq_id = best_candidate_id;
							part_id = 0;

							mst_partitioner.Remove(cur_seq_id);

							cur_no_parts = mst_partitioner.GetNoParts();

							prepare_bit_masks_for_sequence(*sequences[cur_seq_id], ref_seq.p_bit_masks, ref_seq.p_bv_len);
							ref_seq.length = sequences[cur_seq_id]->length;
							ref_seq.sequence_no = sequences[cur_seq_id]->sequence_no;

							if (++c_seq % (100) == 0)
							{
								LOG_DEBUG << "Computing guide tree - " << fixed << setprecision(1)
									<< 100.0 * ((double)c_seq * (2 * n_seq - c_seq + 1) / 2) / ((double)n_seq * (n_seq + 1) / 2) << "%    (" << c_seq << " of " << n_seq << ")  \r";
							}
						}
						else
						{
							++c_seq;
							LOG_DEBUG << "Computing guide tree - " << fixed << setprecision(1)
								<< 100.0 * ((double)c_seq * (2 * n_seq - c_seq + 1) / 2) / ((double)n_seq * (n_seq + 1) / 2) << "\%    (" << c_seq << " of " << n_seq << ")  \r";

							cur_seq_id = -1;
							part_id = no_parts;
						}

						if (flag_value)
							b_flag0.clear(memory_order_relaxed);
						else
							b_flag0.test_and_set(memory_order_relaxed);
					}

					//last_size = -1;

					while (b_flag0.test(memory_order_relaxed) == flag_value)
						;
					flag_value = !flag_value;

					seq_id = cur_seq_id;
					no_parts = cur_no_parts;

					continue;
				}

				auto ids_range = mst_partitioner.GetPart(loc_part_id);

				if (ids_range.first == ids_range.second)
					continue;

				//ids_sub_range.assign(ids_range.first, ids_range.second);
				ids_sub_range.clear();
				for (auto it = ids_range.first; it != ids_range.second; ++it)
					if (*it >= 0)
						ids_sub_range.push_back(*it);

				if (!ids_sub_range.empty())
				{
					part_dist.resize(distance(ids_sub_range.begin(), ids_sub_range.end()));

					calculateDistanceRangeSV<CSequence, CSequenceView, dist_value_t, vector<int>::iterator, decltype(transform)>(
						transform,
						ref_seq,
						sv,
						make_pair(ids_sub_range.begin(), ids_sub_range.end()),
						part_dist.data(),
						lcsbp
						);

					auto p_part_dist = part_dist.begin();
					for (auto p_id = ids_sub_range.begin(); p_id != ids_sub_range.end(); ++p_id, ++p_part_dist)
					{
/*	// Consider future optimization
						double theoretical_best_possible_dist = transform(0.75 * std::min(ref_seq.length, sequences[*p_id]->length), ref_seq.length, sequences[*p_id]->length);

						if (theoretical_best_possible_dist > v_distances_data[*p_id].first)
						{
							no_comp_useless++;
							continue;
						}
						else
							no_comp_necessary++;*/


						if (*p_part_dist <= v_distances_data[*p_id].first)
						{
							//							dist_t s{ *p_part_dist, ~ids_to_uint64(seq_id, *p_id) };
							dist_t s{ *p_part_dist, ~ids_to_uint64(seq_id, *p_id) };

							if (s < v_distances_data[*p_id])
								v_distances_data[*p_id] = s;
						}
					}
				}

				int best_id = *(ids_range.first);
				dist_t best_dist = v_distances_data[best_id];

				for (auto p_id = ids_range.first + 1; p_id != ids_range.second; ++p_id)
					if (*p_id >= 0 && v_distances_data[*p_id] < best_dist)
					{
						best_id = *p_id;
						best_dist = v_distances_data[best_id];
					}

				if (candidate_worker == -1 || best_dist < v_distances_data[candidate_worker])
					candidate_worker = best_id;
			}

//			global_no_comp_necessary += no_comp_necessary;
//			global_no_comp_useless += no_comp_useless;

			}));
	}

	for (auto& t : workers)
		t.join();
	workers.clear();

//	std::cout << "\nNo. comp. necessary: " << global_no_comp_necessary << std::endl;	
//	std::cout << "No. comp. useless  : " << global_no_comp_useless << std::endl;
//	std::cout << "Potential gain     : " << 100.0 * ((double)global_no_comp_useless / (global_no_comp_necessary + global_no_comp_useless)) << "%" << std::endl;

	mst_to_dendogram(mst_edges, v_prim_orders, tree);
}
#else
// *******************************************************************
template <Distance _distance>
void MSTPrim<_distance>::run_view(std::vector<CSequence*>& sequences, tree_structure& tree)
{
	int n_seq = (int)sequences.size();
	int c_seq = 0;

	CSequence ref_seq ("", "");		// Object just to simplify implementation of variants of CLCSBP

	prepare_sequences_view(sequences);

	v_distances.resize(n_seq, dist_t{ std::numeric_limits<dist_value_t>::max(), 0 });

	vector<mst_edge_t> mst_edges;
	mst_edges.reserve(n_seq);	// Just to avoid reallocations

	vector<int> v_prim_orders(n_seq, n_seq);

	vector<int> candidates(n_threads, -1);

	vector<thread> workers;

//	MSTPartitioner mst_partitioner(n_threads, n_threads * 16, 8, n_threads * 2);
	MSTPartitionerNew2 mst_partitioner(n_threads, n_threads * 16, 8, n_threads * 2);

	atomic<uint32_t> part_id{ 0 };
	int cur_seq_id = 0;
	//	int cur_seq_id = n_seq-1;
	int cur_prim_order = 0;
	int cur_no_parts = 0;

	v_prim_orders[cur_seq_id] = cur_prim_order++;

	mst_partitioner.InitPartition(n_seq);

	mst_partitioner.Remove(cur_seq_id);
	prepare_bit_masks_for_sequence(*sequences[cur_seq_id], ref_seq.p_bit_masks, ref_seq.p_bv_len);		// the only allocation of ref_bm is here
	ref_seq.length = sequences[cur_seq_id]->length;
	ref_seq.sequence_no = sequences[cur_seq_id]->sequence_no;

	workers.reserve(n_threads);

	atomic<int> a_cnt0{ (int)n_threads - 1 };
	atomic_flag b_flag0;

	cur_no_parts = mst_partitioner.GetNoParts();

	atomic<uint64_t> global_no_comp_necessary = 0;
	atomic<uint64_t> global_no_comp_useless = 0;

	b_flag0.test_and_set();

	for (int i = 0; i < n_threads; ++i)
	{
		workers.push_back(thread([&, i] {
			auto id_worker = i;
			CLCSBP lcsbp(instruction_set);
			Transform<double, _distance> transform;

			uint64_t no_comp_necessary = 0;
			uint64_t no_comp_useless = 0;

			auto sv = sequence_views;
			auto v_distances_data = v_distances.data();

			vector<dist_value_t> part_dist;

			int& candidate_worker = candidates[id_worker];

			bool flag_value = true;

			vector<int> ids_sub_range;
			vector<int> ids_sub_range_omitted;
			vector<int> ids_sub_range_remaining;

			int seq_id = cur_seq_id;
			int no_parts = cur_no_parts;

			while (true)
			{
				int loc_part_id = part_id++;

//				if (loc_part_id >= no_parts)
				if (loc_part_id >= no_parts && ids_sub_range_remaining.empty())
				{
					if (seq_id < 0)
						break;

					if (!a_cnt0.fetch_sub(1, memory_order_relaxed))
					{
						a_cnt0 = n_threads - 1;

						// All threads are over, so we can pick a new sequence

						int best_candidate_id = -1;
						for (int j = 0; j < n_threads; ++j)
							if (best_candidate_id == -1)
								best_candidate_id = candidates[j];
							else if (candidates[j] != -1)
							{
								if (v_distances_data[candidates[j]] < v_distances_data[best_candidate_id])
									best_candidate_id = candidates[j];
							}

						fill_n(candidates.begin(), n_threads, -1);

 						auto pids = uint64_to_id(~v_distances[best_candidate_id].second);

						mst_edges.emplace_back(pids.first, pids.second, cur_prim_order, -v_distances[best_candidate_id].first);	// We use MaxRangeQueries so change here to negative value

						if (v_prim_orders[pids.first] == n_seq)
							v_prim_orders[pids.first] = cur_prim_order++;
						else
							v_prim_orders[pids.second] = cur_prim_order++;

						bool show_log = false;

						if (!mst_partitioner.IsAlmostEmpty()) [[likely]]
						{
							cur_seq_id = best_candidate_id;
							part_id = 0;

							mst_partitioner.Remove(cur_seq_id);

							cur_no_parts = mst_partitioner.GetNoParts();

							prepare_bit_masks_for_sequence(*sequences[cur_seq_id], ref_seq.p_bit_masks, ref_seq.p_bv_len);
							ref_seq.length = sequences[cur_seq_id]->length;
							ref_seq.sequence_no = sequences[cur_seq_id]->sequence_no;

							if (++c_seq % 100 == 0)
								show_log = true;
						}
						else
						{
							++c_seq;
							show_log = true;

							cur_seq_id = -1;
							part_id = no_parts;
						}

						if (flag_value)
							b_flag0.clear(memory_order_relaxed);
						else
							b_flag0.test_and_set(memory_order_relaxed);

						if (show_log)
						{
							LOG_DEBUG << "Computing guide tree - " << fixed << setprecision(1)
								<< 100.0 * ((double)c_seq * (2 * n_seq - c_seq + 1) / 2) / ((double)n_seq * (n_seq + 1) / 2) << "%    (" << c_seq << " of " << n_seq << ")  \r";
						}
					}

					while (b_flag0.test(memory_order_relaxed) == flag_value)
						;
//						refresh::utils::noop();
					
					flag_value = !flag_value;

					seq_id = cur_seq_id;
					no_parts = cur_no_parts;

					continue;
				}

				auto ids_range = mst_partitioner.GetPart(loc_part_id);

				ids_sub_range.assign(ids_sub_range_remaining.begin(), ids_sub_range_remaining.end());			// avoid std::move() to stop frequent reallocations

//				auto& v_dist_estim_ref = v_dist_estimates[ref_seq.length];

				double theoretical_best_possible_dist = 0;
				uint32_t prev_seq_len = 0;

				for (auto it = ids_range.first; it != ids_range.second; ++it)
					if (*it >= 0)
					{
						if (sv[*it].length != prev_seq_len)		[[unlikely]]
						{
							theoretical_best_possible_dist = transform(std::min(ref_seq.length, sv[*it].length), ref_seq.length, sv[*it].length);
							prev_seq_len = sv[*it].length;
						}
//						double theoretical_best_possible_dist = v_dist_estim_ref[sequences[*it]->length];

						if (theoretical_best_possible_dist <= v_distances_data[*it].first) [[likely]]
							ids_sub_range.push_back(*it);
						else
							ids_sub_range_omitted.push_back(*it);
					}

				size_t batch_size = (loc_part_id < no_parts) ? (ids_sub_range.size() / 8 * 8) : ids_sub_range.size();

				ids_sub_range_remaining.assign(ids_sub_range.begin() + batch_size, ids_sub_range.end());
				ids_sub_range.resize(batch_size);

				if (!ids_sub_range.empty())		[[unlikely]]
				{
					part_dist.resize(distance(ids_sub_range.begin(), ids_sub_range.end()));

					calculateDistanceRangeSV<CSequence, CSequenceView, dist_value_t, vector<int>::iterator, decltype(transform)>(
						transform,
						ref_seq,
						sv,
						make_pair(ids_sub_range.begin(), ids_sub_range.end()),
						part_dist.data(),
						lcsbp
					);

					no_comp_necessary += (uint64_t) ids_sub_range.size();

					int best_id = ids_sub_range.front();
					dist_t best_dist = v_distances_data[best_id];

					auto p_part_dist = part_dist.begin();
					for (auto p_id = ids_sub_range.begin(); p_id != ids_sub_range.end(); ++p_id, ++p_part_dist)
					{
						if (*p_part_dist <= v_distances_data[*p_id].first)		[[unlikely]]
						{
							//							dist_t s{ *p_part_dist, ~ids_to_uint64(seq_id, *p_id) };
							dist_t s{ *p_part_dist, ~ids_to_uint64(seq_id, *p_id) };

							if (s < v_distances_data[*p_id])
								v_distances_data[*p_id] = s;
						}

						if (v_distances_data[*p_id] < best_dist)	[[unlikely]]
						{
							best_id = *p_id;
							best_dist = v_distances_data[best_id];
						}
					}
					
					if (candidate_worker == -1 || best_dist < v_distances_data[candidate_worker])	[[unlikely]]
						candidate_worker = best_id;
				}

				if (!ids_sub_range_omitted.empty())
				{
					no_comp_useless += (uint64_t)ids_sub_range_omitted.size();

					int best_id = ids_sub_range_omitted.front();
					dist_t best_dist = v_distances_data[best_id];

					for (auto p_id = ids_sub_range_omitted.begin() + 1; p_id != ids_sub_range_omitted.end(); ++p_id)
						if (v_distances_data[*p_id] < best_dist)	[[unlikely]]
						{
							best_id = *p_id;
							best_dist = v_distances_data[best_id];
						}
					ids_sub_range_omitted.clear();

					if (candidate_worker == -1 || best_dist < v_distances_data[candidate_worker])	[[unlikely]]
						candidate_worker = best_id;
				}
			}

			global_no_comp_necessary += no_comp_necessary;
			global_no_comp_useless += no_comp_useless;
			}));
	}

	for (auto& t : workers)
		t.join();
	workers.clear();

	LOG_DEBUG << "\nNo. comp. necessary: " << global_no_comp_necessary.load() << std::endl
		<< "No. comp. useless  : " << global_no_comp_useless.load() << std::endl
		<< "Potential gain     : " << 100.0 * ((double) global_no_comp_useless.load() / (global_no_comp_necessary.load() + global_no_comp_useless.load())) << "%" << std::endl;

	mst_to_dendogram(mst_edges, v_prim_orders, tree);
}
#endif

// *******************************************************************
template <Distance _distance>
void MSTPrim<_distance>::prepare_bit_masks_for_sequence(CSequence& seq, bit_vec_t*& bm, uint32_t &p_bv_len)
{
	p_bv_len = (uint32_t)((seq.data_size + bv_size - 1) / bv_size);

	if (bm == nullptr)
		bm = new bit_vec_t[p_bv_len * NO_SYMBOLS];

	fill_n(bm, p_bv_len * NO_SYMBOLS, (bit_vec_t)0);

	for (size_t i = 0; i < seq.length; ++i)
		if (seq.data[i] >= 0 && seq.data[i] < NO_VALID_AMINOACIDS)
			bm[seq.data[i] * p_bv_len + i / bv_size] |= ((bit_vec_t)1) << (i % bv_size);
}

// *******************************************************************
template <Distance _distance>
void MSTPrim<_distance>::run(std::vector<CSequence*>& sequences, tree_structure& tree)
{
	run_view(sequences, tree);
	return;

	int n_seq = (int)sequences.size();
	int c_seq = 0;

	v_distances.resize(n_seq, dist_t{ std::numeric_limits<dist_value_t>::max(), 0});

	vector<mst_edge_t> mst_edges;
	vector<int> v_prim_orders(n_seq, n_seq);

	vector<int> candidates(n_threads, -1);

	vector<thread> workers;

	//	MSTPartitioner mst_partitioner(n_threads, n_threads * 32, 4, n_threads * 4);
	//	MSTPartitioner mst_partitioner(n_threads, n_threads * 16, 4, n_threads * 8);
	MSTPartitioner mst_partitioner(n_threads, n_threads * 16, 4, n_threads * 2);

	atomic<uint32_t> part_id{ 0 };
	int cur_seq_id = 0;
	//	int cur_seq_id = n_seq-1;
	int cur_prim_order = 0;
	int cur_no_parts = 0;

	v_prim_orders[cur_seq_id] = cur_prim_order++;

	mst_partitioner.InitPartition(n_seq);

	mst_partitioner.Remove(cur_seq_id);
	sequences[cur_seq_id]->ComputeBitMasks();

	workers.reserve(n_threads);

	atomic<int> a_cnt0{ (int)n_threads - 1 };
	atomic_flag b_flag0;

	cur_no_parts = mst_partitioner.GetNoParts();

	b_flag0.test_and_set();

	for (int i = 0; i < n_threads; ++i)
	{
		workers.push_back(thread([&, i] {
			auto id_worker = i;
			CLCSBP lcsbp(instruction_set);
//			DistanceToSimilarity<double, _distance> transform;
			Transform<double, _distance> transform;
		
			auto sequences_data = sequences.data();
			auto v_distances_data = v_distances.data();
			CSequence* cur_sequence_ptr = nullptr;

			vector<dist_value_t> part_dist;

			int& candidate_worker = candidates[id_worker];

			//int last_size = -1;

			bool flag_value = true;

			vector<int> ids_sub_range;

			int seq_id = cur_seq_id;
			int no_parts = cur_no_parts;

			cur_sequence_ptr = sequences_data[seq_id];

			while (true)
			{
				int loc_part_id = part_id++;

				if (loc_part_id >= no_parts)
				{
					if (seq_id < 0)
						break;

					if (!a_cnt0.fetch_sub(1, memory_order_relaxed))
					{
						a_cnt0 = n_threads - 1;

						// All threads are over, so we can pick a new sequence
						sequences[cur_seq_id]->ReleaseBitMasks();

						int best_candidate_id = -1;
						for (int j = 0; j < n_threads; ++j)
							if (best_candidate_id == -1)
								best_candidate_id = candidates[j];
							else if (candidates[j] != -1)
							{
								if (v_distances_data[candidates[j]] < v_distances_data[best_candidate_id])
									best_candidate_id = candidates[j];
							}

						fill_n(candidates.begin(), n_threads, -1);

						auto pids = uint64_to_id(~v_distances[best_candidate_id].second);

						mst_edges.emplace_back(pids.first, pids.second, cur_prim_order, -v_distances[best_candidate_id].first);	// We use MaxRangeQueries so change here to negative value

						if (v_prim_orders[pids.first] == n_seq)
							v_prim_orders[pids.first] = cur_prim_order++;
						else
							v_prim_orders[pids.second] = cur_prim_order++;

						if (!mst_partitioner.IsAlmostEmpty())
						{
							cur_seq_id = best_candidate_id;
							part_id = 0;

							mst_partitioner.Remove(cur_seq_id);

							cur_no_parts = mst_partitioner.GetNoParts();

							sequences[cur_seq_id]->ComputeBitMasks();

							if (++c_seq % (100) == 0)
							{
								LOG_DEBUG << "Computing guide tree - " << fixed << setprecision(1)
									<< 100.0 * ((double)c_seq * (2 * n_seq - c_seq + 1) / 2) / ((double)n_seq * (n_seq + 1) / 2) << "%    (" << c_seq << " of " << n_seq << ")  \r";
							}
						}
						else
						{
							++c_seq;
							LOG_DEBUG << "Computing guide tree - " << fixed << setprecision(1)
								<< 100.0 * ((double)c_seq * (2 * n_seq - c_seq + 1) / 2) / ((double)n_seq * (n_seq + 1) / 2) << "\%    (" << c_seq << " of " << n_seq << ")  \r";

							cur_seq_id = -1;
							part_id = no_parts;
						}

						if (flag_value)
							b_flag0.clear(memory_order_relaxed);
						else
							b_flag0.test_and_set(memory_order_relaxed);
					}

					//last_size = -1;

					while (b_flag0.test(memory_order_relaxed) == flag_value)
						;

					flag_value = !flag_value;

					seq_id = cur_seq_id;
					no_parts = cur_no_parts;

					cur_sequence_ptr = sequences_data[seq_id];

					continue;
				}

				auto ids_range = mst_partitioner.GetPart(loc_part_id);

				if (ids_range.first == ids_range.second)
					continue;

				ids_sub_range.assign(ids_range.first, ids_range.second);

				if (!ids_sub_range.empty())
				{
					part_dist.resize(distance(ids_sub_range.begin(), ids_sub_range.end()));

					calculateDistanceRange<CSequence*, dist_value_t, vector<int>::iterator, decltype(transform)>(
						transform,
						cur_sequence_ptr,
						sequences_data,
						make_pair(ids_sub_range.begin(), ids_sub_range.end()),
						part_dist.data(),
						lcsbp
						);

					auto p_part_dist = part_dist.begin();
					for (auto p_id = ids_sub_range.begin(); p_id != ids_sub_range.end(); ++p_id, ++p_part_dist)
					{
						if (*p_part_dist <= v_distances_data[*p_id].first)
						{
							dist_t s{ *p_part_dist, ~ids_to_uint64(seq_id, *p_id) };	// ~ here just for backward compatibility (same order of pairs in case of draws)

							if (s < v_distances_data[*p_id])
								v_distances_data[*p_id] = s;
						}
					}
				}

				int best_id = *(ids_range.first);
				dist_t best_dist = v_distances_data[best_id];

				for (auto p_id = ids_range.first + 1; p_id != ids_range.second; ++p_id)
					if (v_distances_data[*p_id] < best_dist)
					{
						best_id = *p_id;
						best_dist = v_distances_data[best_id];
					}

				if (candidate_worker == -1 || best_dist < v_distances_data[candidate_worker])
					candidate_worker = best_id;
			}
			}));
	}

	for (auto& t : workers)
		t.join();
	workers.clear();


	mst_to_dendogram(mst_edges, v_prim_orders, tree);
}

// *******************************************************************
template <Distance _distance>
void MSTPrim<_distance>::mst_to_dendogram(vector<mst_edge_t>& mst_edges, vector<int>& v_prim_orders, tree_structure& tree)
{
	vector<int> v_rev_prim_orders(v_prim_orders.size());

	mst_edges.emplace(mst_edges.begin(), 0, 0, 0, 0);

	for (int i = 0; i < (int) v_prim_orders.size(); ++i)
		v_rev_prim_orders[v_prim_orders[i]] = i;

	int n_seq = (int)mst_edges.size();

	tree.resize(n_seq * 2 - 1);

	queue<dend_range_t> q_ranges;
	int cur_id = 2 * n_seq - 2;

	q_ranges.emplace(cur_id--, 0, n_seq);

	CMaxRangeQueries<MSTPrim::mst_edge_t> mrq;
	mrq.Init(mst_edges);

	while (!q_ranges.empty())
	{
		auto r = q_ranges.front();
		q_ranges.pop();

		auto p = mrq.MaxElement(r.prim_from + 1, r.prim_to);
		int prim_split = p.prim_order;

		int id_left, id_right;
		 
		if (r.prim_from + 1 == prim_split)
			id_left = v_rev_prim_orders[r.prim_from];
		else
		{
			id_left = cur_id--;
			q_ranges.emplace(id_left, min(r.prim_from, prim_split), max(r.prim_from, prim_split));
		}

		if (prim_split + 1 == r.prim_to)
			id_right = v_rev_prim_orders[prim_split];
		else
		{
			id_right = cur_id--;
			q_ranges.emplace(id_right, min(prim_split, r.prim_to), max(prim_split, r.prim_to));
		}

		tree[r.id] = make_pair(id_left, id_right);
	}
}

// *******************************************************************
template <Distance _distance>
void MSTPrim<_distance>::prepare_sequences_view(std::vector<CSequence*>& sequences)
{
	if (raw_sequence_views)
		free(raw_sequence_views);

	size_t size_sequence_views = sequences.size() * sizeof(CSequenceView);
	size_t raw_size_sequence_views = size_sequence_views + 64;
	auto p = raw_sequence_views = malloc(raw_size_sequence_views);
	sequence_views = (CSequenceView*) my_align(64, size_sequence_views, p, raw_size_sequence_views);

	for (size_t i = 0; i < sequences.size(); ++i)
	{
		sequence_views[i].length = sequences[i]->length;
		sequence_views[i].data = sequences[i]->data;
	}
}

// *******************************************************************
// MSTPartitioner
// *******************************************************************
void MSTPartitioner::InitPartition(int n_elements)
{
	//int cur_part_size;

	min_part_size &= ~0x7;			// round down to mult. of 8
	if (min_part_size == 0)
		min_part_size = 8;

	vd_parts.emplace_back(vector<int>(), 0u, 0u);

	if (min_part_size * n_parts >= n_elements)
	{
		int cur_part_size = (int) min_part_size;

		for (int i = 0; i < n_elements; ++i)
		{
			if ((int) vd_parts.back().data.size() == cur_part_size)
				vd_parts.emplace_back(vector<int>(), 0u, 0u);
			vd_parts.back().data.emplace_back(i);
			++vd_parts.back().i_end;
		}
	}
	else
	{
		int cur_part_size = (int) min_part_size;
	
		double max_part_size = 2.0 * ((double)n_elements - (double)min_part_size * n_tail_parts) / ((double)n_parts - n_tail_parts) - min_part_size;
		double incr = (max_part_size - min_part_size) / ((double)n_parts - n_tail_parts - 1.0);

		double d_cur_part_size = cur_part_size;

		for (int i = 0; i < n_elements; ++i)
		{
			if ((int) vd_parts.back().data.size() == cur_part_size)
			{
				vd_parts.emplace_back(vector<int>(), 0u, 0u);

				if((int64_t) vd_parts.size() > n_tail_parts)
					d_cur_part_size += incr;
				cur_part_size = (int)d_cur_part_size;
				cur_part_size &= ~0x7;
			}

			vd_parts.back().data.emplace_back(i);
			++vd_parts.back().i_end;
		}
	}
}

// *******************************************************************
void MSTPartitioner::Remove(int id)
{
	auto p = lower_bound(vd_parts.begin(), vd_parts.end(), id, [](const part_elem_t& x, int id) {return x.data[x.i_end - 1] < id; });
	int elem_part_id = (int) (p - vd_parts.begin());

	auto q = lower_bound(p->data.begin() + p->i_begin, p->data.begin() + p->i_end, id, [](const int& x, int id) {return x < id; });

	vd_parts[elem_part_id].data.erase(q);
	--vd_parts[elem_part_id].i_end;

	for (size_t i = elem_part_id; i + 1u < vd_parts.size(); ++i)
	{
		vd_parts[i].data.push_back(vd_parts[i + 1u].data[vd_parts[i + 1u].i_begin]);
		++vd_parts[i].i_end;
		++vd_parts[i + 1u].i_begin;

		if (2 * vd_parts[i + 1u].i_begin >= vd_parts[i + 1u].i_end)
		{
			copy(vd_parts[i + 1u].data.begin() + vd_parts[i + 1u].i_begin, vd_parts[i + 1u].data.begin() + vd_parts[i + 1u].i_end, vd_parts[i + 1u].data.begin());
			vd_parts[i + 1u].i_end = vd_parts[i + 1u].i_end - vd_parts[i + 1u].i_begin;
			vd_parts[i + 1u].i_begin = 0u;
			vd_parts[i + 1u].data.resize(vd_parts[i + 1u].i_end);
		}
	}

	if (vd_parts.back().i_begin == vd_parts.back().i_end)
	{
		vd_parts.pop_back();

		if (vd_parts.empty())
			return;

		auto p = max_element(vd_parts.begin(), vd_parts.end(), [](const part_elem_t& x, const part_elem_t& y) {
			return (x.i_end - x.i_begin) < (y.i_end - y.i_begin);
			});

		uint32_t p_size = p->i_end - p->i_begin;
		if (p_size < 16)
			return;

		uint32_t s1 = (p_size / 2u) & ~0x7;

		part_elem_t new_part;
		new_part.data.assign(p->data.begin() + p->i_begin + s1, p->data.end());
		new_part.i_begin = 0;
		new_part.i_end = new_part.data.size();

		p->i_end = (uint32_t)(p->i_begin + s1);
		p->data.resize(p->i_end);
		vd_parts.emplace(next(p), new_part);
	}
}

// *******************************************************************
pair<MSTPartitioner::iterator, MSTPartitioner::iterator> MSTPartitioner::GetPart(int part_id)
{
	if (part_id >= (int) vd_parts.size())
		return make_pair(vd_parts.front().data.begin(), vd_parts.front().data.begin());

	int r_part_id = (int) vd_parts.size() - 1 - part_id;		// Parts are counted from the last

	auto& part = vd_parts[r_part_id];

	return make_pair(part.data.begin() + part.i_begin, part.data.begin() + part.i_end);
}

// *******************************************************************
int MSTPartitioner::GetNoParts()
{
	return (int)vd_parts.size();
}

// *******************************************************************
// True if only one element
bool MSTPartitioner::IsAlmostEmpty()
{
	if (vd_parts.size() > 1)
		return false;

	if (vd_parts.empty())
		return true;

	return vd_parts.front().i_begin + 1 >= vd_parts.front().i_end;
}

// *******************************************************************
// Explicit template specializations for specified distance measures

template class MSTPrim<Distance::indel_div_lcs>;
template class MSTPrim<Distance::indel075_div_lcs>;

// EOF
