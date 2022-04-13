/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

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

	using sim_value_t = array<double, N_CAND>;
#else
	using sim_value_t = double;
#endif

	using sim_t = pair<sim_value_t, uint64_t>;
	vector<sim_t> v_similarities;
	vector<bool> v_processed;

	sim_value_t sim_value_empty;

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
		sim_value_t sim;

#ifdef MANY_CAND
		mst_edge_t(int _seq_from, int _seq_to, int _prim_order, sim_value_t _sim) : seq_from(_seq_from), seq_to(_seq_to), prim_order(_prim_order), sim(_sim) {}
		mst_edge_t() {
			seq_from = -1;
			seq_to = -1;
			prim_order = -1;
			fill(sim.begin(), sim.end(), 0.0);
		}
#else
		mst_edge_t(int _seq_from = -1, int _seq_to = -1, int _prim_order = -1, sim_value_t _sim = 0.0) : seq_from(_seq_from), seq_to(_seq_to), prim_order(_prim_order), sim(_sim) {}
#endif

		bool is_less(const mst_edge_t& x, const mst_edge_t& y)
		{
			if (x.sim != y.sim)
				return x.sim > y.sim;

/*			auto dif_x = abs(x.seq_from - x.seq_to);
			auto dif_y = abs(y.seq_from - y.seq_to);

			if (dif_x != dif_y)
				return dif_x < dif_y;*/

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

//#define EXTENDED_EST_ESTIMATION

	template<typename Transform>
	inline double est_from_hist(Transform &transform, const CSequence& s1, const CSequence& s2)
	{
		int est_lcs = 0;

		auto h1 = s1.hist;
		auto h2 = s2.hist;

		for (int i = 0; i < NO_AMINOACIDS; ++i, ++h1, ++h2)
			est_lcs += min(*h1, *h2);

#ifdef EXTENDED_EST_ESTIMATION
		int no_min1 = 0;
		int no_min2 = 0;

		uint32_t mask = 0;

		h1 = s1.hist;
		h2 = s2.hist;
		int i1, i2;
		uint32_t x;

		for (i1 = i2 = 0; i1 < 8 && i2 < 8;)
		{
			for (; i1 < 8 && no_min1 == no_min2; ++i1)
			{
				symbol_t c1 = s1.data[i1];

				x = 1u << c1;
				if ((mask | x) == mask)
					goto stop_label1;
				mask |= x;

				if (h1[c1] <= h2[c1])
					no_min1++;
			}

			for (; i2 < 8 && no_min1 != no_min2; ++i2)
			{
				symbol_t c2 = s2.data[i2];

				x = 1u << c2;
				if ((mask | x) == mask)
					goto stop_label1;
				mask |= x;

				if (h1[c2] >= h2[c2])
					no_min2++;
			}
		}

	stop_label1:
		est_lcs -= min(no_min1, no_min2);

		no_min1 = no_min2 = 0;
		mask = 0u;

		for (i1 = i2 = 0; i1 < 8 && i2 < 8;)
		{
			for (; i1 < 8 && no_min1 == no_min2; ++i1)
			{
				symbol_t c1 = s1.data[s1.length - i1 - 1];

				x = 1u << c1;
				if ((mask | x) == mask)
					goto stop_label2;
				mask |= x;

				if (h1[c1] <= h2[c1])
					no_min1++;
			}

			for (; i2 < 8 && no_min1 != no_min2; ++i2)
			{
				symbol_t c2 = s2.data[s2.length - i2 - 1];

				x = 1u << c2;
				if ((mask | x) == mask)
					goto stop_label2;
				mask |= x;

				if (h1[c2] >= h2[c2])
					no_min2++;
			}
		}

	stop_label2:
		est_lcs -= min(no_min1, no_min2);
#endif

/*		if (min(no_min1, no_min2))
			cerr << "no_min1: " + to_string(no_min1) + "  no_min2: " + to_string(no_min2) + "\n";*/

		return transform(est_lcs, s1.length, s2.length);

//		est_lcs = min(s1.length, s2.length);

//		est_lcs = lcsbp.EstimateLCS(s1, s2);

/*		int est_indel = s1.length + s2.length - 2 * est_lcs;

		if (est_indel == 0)
			return 1000.0 * est_lcs;

		return (double)est_lcs / est_indel;*/
	}

public:
	MSTPrim(int n_threads, instruction_set_t instruction_set) : AbstractTreeGenerator(n_threads, instruction_set) {
#ifdef MANY_CAND
		fill(sim_value_empty.begin(), sim_value_empty.end(), 0.0);
#else
		sim_value_empty = 0.0;
#endif
	}

	void run(std::vector<CSequence>& sequences, tree_structure& tree) override;
};

