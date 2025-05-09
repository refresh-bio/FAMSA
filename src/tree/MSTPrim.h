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
#include <iterator>
#include <algorithm>
#include "math.h"


class FenwTree 
{
public:
	using iterator = std::vector<int>::iterator;

	FenwTree() : n(0), domain_size(0), removed_count(0), active_count(0) {}

	void init(int n) 
	{
		domain_size = n;
		this->n = n;    
		active_count = n;
		removed_count = 0;
		tree.assign(n + 1, 0);
		elements.resize(n);
		pos.assign(domain_size, -1);

		for (int i = 0; i < n; ++i) {
			elements[i] = i;
			pos[i] = i;
		}

		for (int i = 1; i <= n; ++i)
			update(i, 1);
	}

	void remove(int x) 
	{
		if (x < 0 || x >= domain_size)
			return; 
		int pos_idx = pos[x];
		if (pos_idx == -1)
			return; 

		update(pos_idx + 1, -1);

		elements[pos_idx] = -1;
		pos[x] = -1;
		++removed_count;
		--active_count;

		if (removed_count > elements.size() / 2)
			densify();
	}

	iterator find_k_th(int k) 
	{
		if (k < 0 || k >= active_count)
			return elements.end();

		int idx = findKthIndex(k + 1);
		return elements.begin() + (idx - 1);
	}

	int size() const 
	{
		return active_count;
	}

private:
	int n;                 
	int domain_size;       
	std::vector<int> tree; 
	std::vector<int> elements;
	std::vector<int> pos;     
	int removed_count;     
	int active_count;      

	void update(int i, int delta) 
	{
		int sizeBIT = tree.size() - 1;
		for (; i <= sizeBIT; i += i & -i)
			tree[i] += delta;
	}

	int query(int i) const 
	{
		int sum = 0;
		for (; i > 0; i -= i & -i)
			sum += tree[i];
		return sum;
	}

	int findKthIndex(int k) const 
	{
		int idx = 0;
		int sizeBIT = tree.size() - 1;
		int bit = 1;
		while (bit <= sizeBIT)
			bit <<= 1;
		bit >>= 1;
		for (; bit; bit >>= 1) {
			int next = idx + bit;
			if (next <= sizeBIT && tree[next] < k) {
				k -= tree[next];
				idx = next;
			}
		}
		return idx + 1;
	}

	void densify() 
	{
		std::vector<int> new_elements;

		for (int val : elements)
			if (val != -1)
				new_elements.push_back(val);

		elements = std::move(new_elements);
		n = elements.size();
		removed_count = 0;
		active_count = n;
		tree.assign(n + 1, 0);

		for (int i = 0; i < elements.size(); ++i) 
		{
			int x = elements[i];
			update(i + 1, 1);
			pos[x] = i;
		}
	}
};

class MSTPartitionerNew
{
public:
	using iterator = std::vector<int>::iterator;

private:
	FenwTree fenw_tree;
	int64_t n_threads;
	int64_t n_parts;
	int64_t min_part_size;
	int64_t n_tail_parts;

	int cur_part_size = 0;

	void update_cur_part_size()
	{
		cur_part_size = (fenw_tree.size() / n_parts) / min_part_size * min_part_size;
		if (cur_part_size < min_part_size)
			cur_part_size = min_part_size;
	}

public:
	MSTPartitionerNew(int _n_threads, int _n_parts, int _min_part_size, int _n_tail_parts) : n_threads(_n_threads), n_parts(_n_parts), min_part_size(_min_part_size), n_tail_parts(_n_tail_parts)
	{
	}

	void InitPartition(int n_elements)
	{
		fenw_tree.init(n_elements);
		update_cur_part_size();
	}

	void Remove(int id)
	{
		fenw_tree.remove(id);
		update_cur_part_size();
	}

	std::pair<iterator, iterator> GetPart(int part_id)
	{
		int begin = part_id * cur_part_size;
		int end = (part_id + 1) * cur_part_size;
		if (end > (int)fenw_tree.size())
			end = (int)fenw_tree.size();

		return std::make_pair(fenw_tree.find_k_th(begin), fenw_tree.find_k_th(end));
	}

	int GetNoParts()
	{
		return (fenw_tree.size() + cur_part_size - 1) / cur_part_size;
	}

	bool IsAlmostEmpty()
	{
		return fenw_tree.size() <= 1;
	}
};

class MSTPartitionerNew2
{
public:
	using iterator = std::vector<int>::iterator;

private:
	vector<int> items;
	vector<int> pos;
	vector<int> part_boundaries;

	int64_t n_present{};
	int64_t n_threads;
	int64_t n_parts;
	int64_t min_part_size;

//	int cur_part_size = 0;

	void update_part_boundaries()
	{
		part_boundaries.clear();
		part_boundaries.reserve(n_parts + 2);
		part_boundaries.emplace_back(0);

		if (n_parts * min_part_size >= items.size())
		{
			for(int i = 1; i < n_parts && part_boundaries.back() < items.size(); ++i)
				part_boundaries.emplace_back(i * min_part_size);
		}
		else
		{
			int64_t N = int64_t(items.size());
			int64_t a = min_part_size;

			for(int64_t n = n_parts; n; --n)
			{
				int64_t b = 2 * N / n - a;
				if (b >= a)
				{
					part_boundaries.emplace_back(part_boundaries.back() + int(b));
					N -= b;
				}
			}
		}

		while (part_boundaries.back() > items.size())
			part_boundaries.pop_back();
		if (part_boundaries.back() != items.size())
			part_boundaries.emplace_back(int(items.size()));
	}

	void densify()
	{
		vector<int> new_items;
		new_items.reserve(n_present);

		for (int val : items)
			if (val != -1)
				new_items.push_back(val);

		items = std::move(new_items);

		fill(pos.begin(), pos.end(), -1);
		for (int i = 0; i < n_present; ++i)
			pos[items[i]] = i;

		update_part_boundaries();
	}

public:
	MSTPartitionerNew2(int _n_threads, int _n_parts, int _min_part_size, int _n_tail_parts) : n_threads(_n_threads), n_parts(_n_parts), min_part_size(_min_part_size)
	{
	}

	void InitPartition(int n_elements)
	{
		n_present = n_elements;

		items.resize(n_elements);
		pos.resize(n_elements);

		for (int i = 0; i < n_elements; ++i) 
		{
			items[i] = i;
			pos[i] = i;
		}

		update_part_boundaries();
	}

	void Remove(int id)
	{
		items[pos[id]] = -1;
		pos[id] = -1;

		--n_present;

		if(n_present * 2 <= items.size())
			densify();
	}

	std::pair<iterator, iterator> GetPart(int part_id)
	{
		if (part_id + 1 >= (int)part_boundaries.size())
			return std::make_pair(items.end(), items.end());

		return std::make_pair(items.begin() + part_boundaries[part_id], items.begin() + part_boundaries[part_id + 1]);
	}

	int GetNoParts()
	{
		return int(part_boundaries.size()) - 1;
	}

	bool IsAlmostEmpty()
	{
		return items.size() <= 1;
	}
};

class MSTPartitioner
{
public:
	using iterator = vector<int>::iterator;

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
	using dist_value_t = double;
//	using dist_value_t = uint64_t;

	using dist_t = pair<dist_value_t, uint64_t>;
	vector<dist_t> v_distances;
	vector<bool> v_processed;
	bool very_verbose_mode;

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

		mst_edge_t(int _seq_from = -1, int _seq_to = -1, int _prim_order = -1, dist_value_t _dist = 0) : seq_from(_seq_from), seq_to(_seq_to), prim_order(_prim_order), dist(_dist) {}

		bool is_less(const mst_edge_t& x, const mst_edge_t& y) const
		{
			if (x.dist != y.dist)
				return x.dist > y.dist;

			return ids_to_uint64(x.seq_from, x.seq_to) > ids_to_uint64(y.seq_from, y.seq_to);
		}

		bool operator<(const mst_edge_t& x) const {
			return is_less(*this, x);
		}

		bool operator>(const mst_edge_t& x) const {
			return is_less(x, *this);
		}

		bool operator==(const mst_edge_t& x) const {
			return !is_less(*this, x) && !is_less(x, *this);
		}
		
		bool operator!=(const mst_edge_t& x) const {
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
	MSTPrim(int n_threads, instruction_set_t instruction_set) : 
		AbstractTreeGenerator(n_threads, instruction_set),
		very_verbose_mode(very_verbose_mode)
	{
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

