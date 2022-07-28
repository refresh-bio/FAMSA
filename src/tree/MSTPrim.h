/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#pragma once

#include "AbstractTreeGenerator.h"
#include "IPartialGenerator.h"
#include "../lcs/lcsbp.h"

#include <array>
#include <vector>
#include <stack>
#include <utility>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <tuple>
#include <list>
#include <limits>
#include <functional>
#include "math.h"

//#define MANY_CAND	10

class MSTPartitioner
{
public:
	typedef vector<int>::iterator iterator;

private:
	int64_t n_threads;
	int64_t n_parts;
	int64_t min_part_size;
	int64_t n_tail_parts;
	
	struct part_elem_t { 
		vector<int> data; 
		uint32_t i_begin;
		uint32_t i_end;

		part_elem_t() : i_begin(0), i_end(0) {};

		part_elem_t(const vector<int> &_data, uint32_t _i_begin, uint32_t _i_end) : data(_data), i_begin(_i_begin), i_end(_i_end) {};

		part_elem_t(const part_elem_t&) = default;
		part_elem_t(part_elem_t&&) = default;

		part_elem_t& operator=(const part_elem_t& x) noexcept = delete;
		part_elem_t& operator=(part_elem_t&& x) noexcept = default;
	};
	vector<part_elem_t> vd_parts;

public:
	MSTPartitioner(int _n_threads, int _n_parts, int _min_part_size, int _n_tail_parts) : n_threads(_n_threads), n_parts(_n_parts), min_part_size(_min_part_size), n_tail_parts(_n_tail_parts) {};
	void InitPartition(int n_elements);
	void Remove(int id);
	pair<iterator, iterator> GetPart(int part_id);
	int GetNoParts();
	bool IsAlmostEmpty();
};

template<typename T>
class CMaxRangeQueries
{
	int n_levels;
	vector<vector<T>> vv_data;

	void prepare(const vector<T>& v_data) {
		size_t n = v_data.size();

		n_levels = 1 + (int)log2(n);

		vv_data.resize(n_levels);

		vv_data[0] = v_data;

		for (int i = 1; i < n_levels; ++i)
		{
			auto& vc = vv_data[i];
			auto& vp = vv_data[i - 1u];

			vc.resize(n - (1ull << i) + 1ull);

			for (uint32_t j = 0; j + (1u << i) - 1 < n; ++j)
				if (vp[j] > vp[j + (1u << (i - 1))])
					vc[j] = vp[j];
				else
					vc[j] = vp[j + (1u << (i - 1))];
		}
	}

public:
	CMaxRangeQueries() {
		n_levels = 0;
	}

	void Init(const vector<T>& v_data) {
		prepare(v_data);
	}

	T MaxElement(int begin, int end) {
		uint32_t lev = (uint32_t)log2(end - begin);

		if (vv_data[lev][begin] > vv_data[lev][end - (1u << lev)])
			return vv_data[lev][begin];
		else
			return vv_data[lev][end - (1u << lev)];
	}
};

template <Distance _distance>
class MSTPrim : public AbstractTreeGenerator {
#ifdef MANY_CAND
	static const int N_CAND = MANY_CAND;

	using dist_value_t = array<double, N_CAND>;
#else
	using dist_value_t = double;
//	using dist_value_t = uint64_t;
#endif

	using dist_t = pair<dist_value_t, uint64_t>;
	vector<dist_t> v_distances;
	vector<bool> v_processed;

	void* raw_sequence_views;
	CSequenceView* sequence_views;

	static uint64_t ids_to_uint64(int id1, int id2)
	{
		if (id1 < 0 || id2 < 0)
			return 0u;
		if (id1 > id2)
			return (((uint64_t)id2) << 32) + (uint64_t)id1;
		return (((uint64_t)id1) << 32) + (uint64_t)id2;
	}

	constexpr pair<int, int> uint64_to_id(uint64_t packed_ids) 
	{
		int id1 = (int)(packed_ids >> 32);
		int id2 = (int)(packed_ids & 0xffffffffull);

		if (id1 < id2)
			return make_pair(id1, id2);
		else
			return make_pair(id2, id1);
	}

	struct mst_edge_t {
		int seq_from;
		int seq_to;
		int prim_order;
		dist_value_t dist;

#ifdef MANY_CAND
		mst_edge_t(int _seq_from, int _seq_to, int _prim_order, dist_value_t _dist) : seq_from(_seq_from), seq_to(_seq_to), prim_order(_prim_order), dist(_dist) {}
		mst_edge_t() {
			seq_from = -1;
			seq_to = -1;
			prim_order = -1;
			fill(sim.begin(), sim.end(), numeric_limits<double>::max());
		}
#else
		mst_edge_t(int _seq_from = -1, int _seq_to = -1, int _prim_order = -1, dist_value_t _dist = 0) : seq_from(_seq_from), seq_to(_seq_to), prim_order(_prim_order), dist(_dist) {}
#endif

		bool is_less(const mst_edge_t& x, const mst_edge_t& y)
		{
			if (x.dist != y.dist)
				return x.dist > y.dist;

			return ids_to_uint64(x.seq_from, x.seq_to) > ids_to_uint64(y.seq_from, y.seq_to);
		}

		bool operator<(const mst_edge_t& x) {
			return is_less(*this, x);
		}

		bool operator>(const mst_edge_t& x) {
			return is_less(x, *this);
		}

		bool operator==(const mst_edge_t& x) {
			return !is_less(*this, x) && !is_less(x, *this);
		}
		
		bool operator!=(const mst_edge_t& x) {
			return is_less(*this, x) || is_less(x, *this);
		}
	};

	struct dend_range_t {
		int id;
		int prim_from;
		int prim_to;

		dend_range_t(int _id, int _prim_from, int _prim_to) : id(_id), prim_from(_prim_from), prim_to(_prim_to) {}
	};

	void mst_to_dendogram(vector<mst_edge_t>& mst_edges, vector<int>& v_prim_orders, tree_structure& tree);
	void prepare_sequences_view(std::vector<CSequence*>& sequences);
	void prepare_bit_masks_for_sequence(CSequence& seq, bit_vec_t*& bm, uint32_t& p_bv_len);

public:
	MSTPrim(int n_threads, instruction_set_t instruction_set) : AbstractTreeGenerator(n_threads, instruction_set) {
#ifdef MANY_CAND
//		fill(sim_value_empty.begin(), sim_value_empty.end(), 0.0);
#else
//		sim_value_empty = 0.0;
#endif

		sequence_views = nullptr;
		raw_sequence_views = nullptr;
	}

	~MSTPrim()
	{
		if (raw_sequence_views)
			free(raw_sequence_views);
	}

	void run(std::vector<CSequence*>& sequences, tree_structure& tree) override;

	void run_view(std::vector<CSequence*>& sequences, tree_structure& tree);
};

