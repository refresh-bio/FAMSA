/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#include "MSTPrim.h"
#include "AbstractTreeGenerator.hpp"
#include "SingleLinkage.h"

#include "../utils/log.h"
#include "../core/queues.h"

#include <algorithm>
#include <atomic>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <random>
#include <thread>
#include <utility>


//#define ATOMIC_WAIT
#define USE_HISTOGRAMS

//#define OLD_ATOMIC_FLAG

//#define USE_VERTEX_HEIGHT

#include <chrono>

using namespace std;

#if 1
//#pragma GCC optimize("align-loops=16")
// *******************************************************************
template <Distance _distance>
void MSTPrim<_distance>::run(std::vector<CSequence>& sequences, tree_structure& tree)
{
	int n_seq = (int)sequences.size();
	int c_seq = 0;

	//auto t1 = chrono::high_resolution_clock::now();
	//auto sc1 = chrono::system_clock::now();

#ifdef MANY_CAND
	sim_value_t s0;
	for (auto& x : s0)
		x = -1.0;

	v_similarities.resize(n_seq, sim_t{ s0, 0 });
#else
	v_similarities.resize(n_seq, sim_t{ -1.0, 0 });
#endif

#ifdef USE_VERTEX_HEIGHT
	vector<uint32_t> v_vertex_height(sequences.size(), 0u);
#endif

	for (auto& seq : sequences)
		seq.PrepareHistogram();

	vector<mst_edge_t> mst_edges;
	vector<int> v_prim_orders(n_seq, n_seq);

	vector<int> candidates(n_threads, -1);
	atomic<int> n_thr_ready{ 0 };

	vector<thread> workers;

	//	MSTPartitioner mst_partitioner(n_threads, n_threads * 32, 4, n_threads * 4);
	//	MSTPartitioner mst_partitioner(n_threads, n_threads * 16, 4, n_threads * 8);
	MSTPartitioner mst_partitioner(n_threads, n_threads * 16, 4, n_threads * 2);

	atomic<uint32_t> part_id{ 0 };
	int cur_seq_id = 0;
	//	int cur_seq_id = n_seq-1;
	int cur_prim_order = 0;
	int cur_no_parts = 0;

	map<pair<int, int>, double> m_sim;

	v_prim_orders[cur_seq_id] = cur_prim_order++;

	mst_partitioner.InitPartition(n_seq);

	mst_partitioner.Remove(cur_seq_id);
	sequences[cur_seq_id].ComputeBitMasks();

	workers.reserve(n_threads);

	atomic<int> a_cnt0{ (int)n_threads - 1 };
#ifndef OLD_ATOMIC_FLAG
	atomic_flag b_flag0;
#else
	atomic<bool> b_flag0;
#endif
	// barrier bar_flag0(n_threads);

	atomic<uint64_t> n_filtered_pos{ 0 };
	atomic<uint64_t> n_filtered_neg{ 0 };

	cur_no_parts = mst_partitioner.GetNoParts();

#ifndef OLD_ATOMIC_FLAG
	b_flag0.test_and_set();
#else
	b_flag0.store(true);
#endif

	for (uint32_t i = 0; i < n_threads; ++i)
	{
		workers.push_back(thread([&, i] {
			auto id_worker = i;
			CLCSBP lcsbp(instruction_set);
			DistanceToSimilarity<double, _distance> transform;
		
			auto sequences_data = sequences.data();
			auto v_similarities_data = v_similarities.data();
			CSequence* cur_sequence_ptr = nullptr;

			vector<double> part_sim;

			int& candidate_worker = candidates[id_worker];

			//int last_size = -1;

			uint64_t loc_n_filtered_pos = 0;
			uint64_t loc_n_filtered_neg = 0;

			bool flag_value = true;

			vector<int> ids_sub_range;

			int seq_id = cur_seq_id;
			int no_parts = cur_no_parts;

			cur_sequence_ptr = &sequences_data[seq_id];

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
						sequences[cur_seq_id].ReleaseBitMasks();

						int best_candidate_id = -1;
						for (int j = 0; j < n_threads; ++j)
							if (best_candidate_id == -1)
								best_candidate_id = candidates[j];
							else if (candidates[j] != -1)
							{
								if (v_similarities_data[candidates[j]] > v_similarities_data[best_candidate_id])
									best_candidate_id = candidates[j];
							}

						fill_n(candidates.begin(), n_threads, -1);

						auto pids = uint64_to_id(v_similarities[best_candidate_id].second);

						mst_edges.emplace_back(pids.first, pids.second, cur_prim_order, v_similarities[best_candidate_id].first);

						if (v_prim_orders[pids.first] == n_seq)
							v_prim_orders[pids.first] = cur_prim_order++;
						else
						{
							v_prim_orders[pids.second] = cur_prim_order++;
#ifdef USE_VERTEX_HEIGHT
							v_vertex_height[pids.second] = v_vertex_height[pids.first] + 1;
#endif
						}

						if (!mst_partitioner.IsAlmostEmpty())
						{
							cur_seq_id = best_candidate_id;
							part_id = 0;

							mst_partitioner.Remove(cur_seq_id);

							cur_no_parts = mst_partitioner.GetNoParts();

							sequences[cur_seq_id].ComputeBitMasks();

							if (++c_seq % (100) == 0)
							{
								LOG_DEBUG << "Computing guide tree - " << fixed << setprecision(1)
									<< 100.0 * ((double)c_seq * (2 * n_seq - c_seq + 1) / 2) / ((double)n_seq * (n_seq + 1) / 2) << "\%    (" << c_seq << " of " << n_seq << ")  \r";
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

#ifndef OLD_ATOMIC_FLAG
						if (flag_value)
							b_flag0.clear(memory_order_relaxed);
						else
							b_flag0.test_and_set(memory_order_relaxed);
#else
						if (flag_value)
							b_flag0.store(false, memory_order_relaxed);
						else
							b_flag0.store(true, memory_order_relaxed);
//						b_flag0.store(!flag_value, memory_order_relaxed);
#endif
#ifdef ATOMIC_WAIT
						b_flag0.notify_all();
#endif
					}

					//last_size = -1;

#ifdef ATOMIC_WAIT
					b_flag0.wait(flag_value, memory_order_relaxed);
#else
#ifndef OLD_ATOMIC_FLAG
					while (b_flag0.test(memory_order_relaxed) == flag_value)
#else
					while (b_flag0.load(memory_order_relaxed) == flag_value)
#endif
						;
#endif
					flag_value = !flag_value;

					seq_id = cur_seq_id;
					no_parts = cur_no_parts;

					cur_sequence_ptr = &sequences_data[seq_id];

					continue;
				}

				auto ids_range = mst_partitioner.GetPart(loc_part_id);

				if (ids_range.first == ids_range.second)
					continue;

				ids_sub_range.clear();
				for (auto p = ids_range.first; p != ids_range.second; ++p)
				{
#ifdef USE_HISTOGRAMS
#ifdef MANY_CAND
					if (est_from_hist(transform, sequences_data[seq_id], sequences_data[*p]) >= v_similarities_data[*p].first.back())
#else
					if (est_from_hist(transform, sequences_data[seq_id], sequences_data[*p]) >= v_similarities_data[*p].first)
#endif
#endif
					{
						ids_sub_range.emplace_back(*p);
						loc_n_filtered_pos++;
					}
#ifdef USE_HISTOGRAMS
					else
						loc_n_filtered_neg++;
#endif
				}

				if (!ids_sub_range.empty())
				{
					part_sim.resize(distance(ids_sub_range.begin(), ids_sub_range.end()));

					calculateDistanceRange<CSequence, double, vector<int>::iterator, decltype(transform)>(
						transform,
						*cur_sequence_ptr,
						sequences_data,
						make_pair(ids_sub_range.begin(), ids_sub_range.end()),
						part_sim.data(),
						lcsbp
						);

					auto p_part_sim = part_sim.begin();
					for (auto p_id = ids_sub_range.begin(); p_id != ids_sub_range.end(); ++p_id, ++p_part_sim)
					{
#ifdef MANY_CAND
						if (*p_part_sim >= v_similarities_data[*p_id].first.back())
						{
							if (*p_part_sim > v_similarities_data[*p_id].first.front() ||
								(*p_part_sim == v_similarities_data[*p_id].first.front() && ids_to_uint64(seq_id, *p_id) > v_similarities_data[*p_id].second))
							{
//								copy_backward(v_similarities_data[*p_id].first.begin(), v_similarities_data[*p_id].first.end() - 1, v_similarities_data[*p_id].first.end() - 1);
								copy_backward(v_similarities_data[*p_id].first.begin(), v_similarities_data[*p_id].first.end() - 1, v_similarities_data[*p_id].first.end());
								v_similarities_data[*p_id].first.front() = *p_part_sim;
								v_similarities_data[*p_id].second = ids_to_uint64(seq_id, *p_id);
							}
							else
							{
								v_similarities_data[*p_id].first.back() = *p_part_sim;
								for (int i = N_CAND - 1; i > 0; --i)
									if (v_similarities_data[*p_id].first[i] > v_similarities_data[*p_id].first[i - 1])
										swap(v_similarities_data[*p_id].first[i], v_similarities_data[*p_id].first[i - 1]);
									else
										break;
							}
						}

#else
						if (*p_part_sim >= v_similarities_data[*p_id].first)
						{
							sim_t s{ *p_part_sim, ids_to_uint64(seq_id, *p_id) };

							if (s > v_similarities_data[*p_id])
								v_similarities_data[*p_id] = s;
						}
#endif
					}
				}

				int best_id = *(ids_range.first);
				sim_t best_sim = v_similarities_data[best_id];

				for (auto p_id = ids_range.first + 1; p_id != ids_range.second; ++p_id)
					if (v_similarities_data[*p_id] > best_sim)
					{
						best_id = *p_id;
						best_sim = v_similarities_data[best_id];
					}

				if (candidate_worker == -1 || best_sim > v_similarities_data[candidate_worker])
					candidate_worker = best_id;
			}

			n_filtered_pos += loc_n_filtered_pos;
			n_filtered_neg += loc_n_filtered_neg;
			}));
	}

	//auto c1 = chrono::high_resolution_clock::now();

	for (auto& t : workers)
		t.join();
	workers.clear();

	//auto t2 = chrono::high_resolution_clock::now();
	//auto dur = chrono::duration_cast<chrono::duration<double>>(t2 - t1);

	LOG_VERBOSE 
		<< "No. filtered pos.   : " << (uint64_t)n_filtered_pos << endl
		<< "No. filtered neg.   : " << (uint64_t)n_filtered_neg << endl
		<< "Fract. filtered neg.: " << (double)n_filtered_neg / (n_filtered_pos + n_filtered_neg) << endl;

	mst_to_dendogram(mst_edges, v_prim_orders, tree);
}
#endif

// *******************************************************************
template <Distance _distance>
void MSTPrim<_distance>::mst_to_dendogram(vector<mst_edge_t>& mst_edges, vector<int>& v_prim_orders, tree_structure& tree)
{
	vector<int> v_rev_prim_orders(v_prim_orders.size());

#ifdef MANY_CAND
	mst_edges.emplace(mst_edges.begin(), 0, 0, 0, sim_value_empty);
#else
	mst_edges.emplace(mst_edges.begin(), 0, 0, 0, 0);
#endif

	for (int i = 0; i < v_prim_orders.size(); ++i)
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

#if 0
// *******************************************************************
// CSingleLinkageQueue
// *******************************************************************
//CMSTPrimQueue::CMSTPrimQueue(vector<CSequence>* _sequences, uint32_t _n_rows, uint32_t _max_buffered_rows)
CMSTPrimQueue::CMSTPrimQueue(int _n_threads)
{
	n_threads = _n_threads;
	eoq_flag = false;
}

// *******************************************************************
CMSTPrimQueue::~CMSTPrimQueue()
{
}

// *******************************************************************
bool CMSTPrimQueue::PopTask(int& seq_id, int &from, int &to)
{
	unique_lock<mutex> lck(mtx_tasks);
	cv_tasks.wait(lck, [this] {return !q_tasks.empty() || this->eoq_flag; });

	if (eoq_flag)
		return false;	// End of data in the profiles queue

	seq_id = get<0>(q_tasks.top());
	from   = get<1>(q_tasks.top());
	to     = get<2>(q_tasks.top());

	q_tasks.pop();

#ifdef PRODUCE_LOG
	cerr << "GetTask : " << row_id << "\n";
	cerr << available_buffers.size() << "  " << lowest_uncomputed_row << "\n";
#endif

	return true;
}

// *******************************************************************
void CMSTPrimQueue::CMSTPrimQueue::PushTasks(int seq_id, int n_seq, int part_size)
{
	unique_lock<mutex> lck(mtx_tasks);

//	for(int i = 0; i < n_seq; i += part_size)
//		q_tasks.emplace(seq_id, i, (i + part_size < n_seq) ? i + part_size : n_seq);
	
	for(int i = n_seq; i >= 0; i -= part_size)
		q_tasks.emplace(seq_id, (i - part_size < 0) ? 0 : i - part_size, i);

	n_parts_to_process = (int) q_tasks.size();

	cv_tasks.notify_all();
}

// *******************************************************************
void CMSTPrimQueue::MarkTaskReady()
{
	unique_lock<mutex> lck(mtx_tasks);

	if (--n_parts_to_process == 0)
		cv_ready.notify_one();
}

// *******************************************************************
void CMSTPrimQueue::AllTasksReady()
{
	unique_lock<mutex> lck(mtx_tasks);
	cv_ready.wait(lck, [this] {return n_parts_to_process == 0; });
}

// *******************************************************************
void CMSTPrimQueue::SetEoq()
{
	unique_lock<mutex> lck(mtx_tasks);

	eoq_flag = true;
	cv_tasks.notify_all();
}
#endif

// *******************************************************************
// MSTPartitioner
// *******************************************************************
void MSTPartitioner::InitPartition(int n_elements)
{
	//int cur_part_size;

	min_part_size &= ~0x3;			// round down to mult. of 4
	if (min_part_size == 0)
		min_part_size = 4;

	vd_parts.emplace_back(vector<int>(), 0u, 0u);

	if (min_part_size * n_parts >= n_elements)
	{
		int cur_part_size = min_part_size;

		for (int i = 0; i < n_elements; ++i)
		{
			if (vd_parts.back().data.size() == cur_part_size)
				vd_parts.emplace_back(vector<int>(), 0u, 0u);
			vd_parts.back().data.emplace_back(i);
			++vd_parts.back().i_end;
		}
	}
	else
	{
		int cur_part_size = min_part_size;
	
		double max_part_size = 2.0 * ((double)n_elements - (double)min_part_size * n_tail_parts) / ((double)n_parts - n_tail_parts) - min_part_size;
		double incr = (max_part_size - min_part_size) / ((double)n_parts - n_tail_parts - 1.0);

		double d_cur_part_size = cur_part_size;

		for (int i = 0; i < n_elements; ++i)
		{
			if (vd_parts.back().data.size() == cur_part_size)
			{
				vd_parts.emplace_back(vector<int>(), 0u, 0u);

				if(vd_parts.size() > n_tail_parts)
					d_cur_part_size += incr;
				cur_part_size = (int)d_cur_part_size;
				cur_part_size &= ~0x3;
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
	int elem_part_id = p - vd_parts.begin();

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

		uint32_t s1 = (p_size / 2u) & ~0x3;

		auto q = next(p);
		vd_parts.emplace(q, part_elem_t());
		q->data.assign(p->data.begin() + p->i_begin + s1, p->data.end());
		q->i_end = q->data.size();
		p->i_end = p->i_begin + s1;
		p->data.resize(p->i_end);
	}
}

// *******************************************************************
pair<MSTPartitioner::iterator, MSTPartitioner::iterator> MSTPartitioner::GetPart(int part_id)
{
	if (part_id >= vd_parts.size())
		return make_pair(vd_parts.front().data.begin(), vd_parts.front().data.begin());

	int r_part_id = vd_parts.size() - 1 - part_id;		// Parts are counted from the last

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
template class MSTPrim<Distance::sqrt_indel_div_lcs>;



// EOF
