/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#pragma once
#include "AbstractTreeGenerator.h"

#include "../lcs/lcsbp.h"
#include "../utils/utils.h"

#undef min
#undef max
#include <cmath>
#include <type_traits>
#include <algorithm>

// overloads for converting sequence type to pointer
inline CSequence* seq_to_ptr(CSequence* x) { return x; }
inline CSequence* seq_to_ptr(CSequence& x) { return &x; }

inline CSequenceView* sv_to_ptr(CSequenceView* x) { return x; }
inline CSequenceView* sv_to_ptr(CSequenceView& x) { return &x; }

// dummy implementation
template <class T, Distance measure>
struct Transform {
	//static_assert(0, "Cannot use dummy implementation");
	
	T operator()(uint32_t lcs, uint32_t len1, uint32_t len2) {
		return 0;
	}
};

template <class T>
struct Transform<T, Distance::indel075_div_lcs>{
private:
	std::vector<T> pp_pow075_rec;
	uint32_t cur_pp_size = 0;

	void pp_extend(uint32_t val)
	{
		pp_pow075_rec.resize(val + 1);
		for (; cur_pp_size <= val; ++cur_pp_size)
			pp_pow075_rec[cur_pp_size] = (T) pow(cur_pp_size, 0.75);
	}

public:
	T operator()(uint32_t lcs, uint32_t len1, uint32_t len2) { 
		T indel = (T) (len1 + len2 - 2 * lcs);
		T l = (T) lcs;

		if (indel >= cur_pp_size)
			pp_extend((uint32_t) indel);

		if(l)
			return pp_pow075_rec[indel] / l;
		else
			return nextafter(numeric_limits<T>::max(), 0);
	}
};

template <class T>
struct Transform<T, Distance::indel_div_lcs> {
	T operator()(uint32_t lcs, uint32_t len1, uint32_t len2) {
		T indel = (T) (len1 + len2 - 2 * lcs);

		if (lcs)
			return (T)indel / lcs;
		else
			return nextafter(numeric_limits<T>::max(), 0);
	}
};

template <class T>
struct Transform<T, Distance::pairwise_identity> {
	T operator()(uint32_t lcs, uint32_t len1, uint32_t len2) {
		return (T)lcs / std::min(len1, len2);
	}
};

/*template <class T>
struct Transform<T, Distance::neg_lcs_div_indel> {
	T operator()(uint32_t lcs, uint32_t len1, uint32_t len2) {
		T indel = len1 + len2 - 2 * lcs;
		return indel == 0 ? (-(T)lcs * 1000) : (-(T)lcs / indel); 
	}
};

template <class T>
struct Transform<T, Distance::neg_lcs_div_minlen> {
	T operator()(uint32_t lcs, uint32_t len1, uint32_t len2) {
		return 1.0 - (T)lcs / std::min(len1, len2);
	}
};

template <class T>
struct Transform<T, Distance::neg_lcs_div_len_corrected> {
	T operator()(uint32_t lcs, uint32_t len1, uint32_t len2) {
		// make len1 the longer 
		if (len1 < len2) { std::swap(len1, len2);  }
		
		T d = 1.0 - (T)lcs / len2;

		// MAFFT distance correction
		auto correction = [](T x, T y) ->T { return y / x * 0.1 + 10000 / (x + 10000) + 0.01;  };
		d = d / correction(len1, len2);
		return d;
	}
};*/

template<class T, Distance measure>
struct DistanceToSimilarity {
	Transform<T, measure> transform;
		
	T operator()(uint32_t lcs, uint32_t len1, uint32_t len2) {
		auto val = transform(lcs, len1, len2);
		return (val == 0) ? (lcs * 1000) : (1.0 / val);
	}
};

// *******************************************************************
/*
	seq_type can be:
	- CSequence,
	- CSequence*,
*/
template <class seq_type, class distance_type, typename Transform>
void AbstractTreeGenerator::calculateDistanceVector(
	Transform& transform,
	seq_type& ref,
	seq_type* sequences, 
	int n_seqs, 
	distance_type* out_vector, 
	CLCSBP& lcsbp)
{
	uint32_t lcs_lens[8];
	
	seq_to_ptr(ref)->ComputeBitMasks();

	// process portions of 8 sequences
	for (int j = 0; j < n_seqs / 8; ++j) {
		lcsbp.GetLCSBP(
			seq_to_ptr(ref),
			seq_to_ptr(sequences[j * 8 + 0]),
			seq_to_ptr(sequences[j * 8 + 1]),
			seq_to_ptr(sequences[j * 8 + 2]),
			seq_to_ptr(sequences[j * 8 + 3]),
			seq_to_ptr(sequences[j * 8 + 4]),
			seq_to_ptr(sequences[j * 8 + 5]),
			seq_to_ptr(sequences[j * 8 + 6]),
			seq_to_ptr(sequences[j * 8 + 7]),
			lcs_lens);

		for (int k = 0; k < 8; ++k) {
			out_vector[j * 8 + k] = transform(lcs_lens[k], seq_to_ptr(ref)->length, seq_to_ptr(sequences[j * 8 + k])->length);
		}
	}

	// if there is something left
	int n_processed = n_seqs / 8 * 8;
	if (n_processed < n_seqs) {
		lcsbp.GetLCSBP(
			seq_to_ptr(ref),
			(n_processed + 0 < n_seqs) ? seq_to_ptr(sequences[n_processed + 0]) : nullptr,
			(n_processed + 1 < n_seqs) ? seq_to_ptr(sequences[n_processed + 1]) : nullptr,
			(n_processed + 2 < n_seqs) ? seq_to_ptr(sequences[n_processed + 2]) : nullptr,
			(n_processed + 3 < n_seqs) ? seq_to_ptr(sequences[n_processed + 3]) : nullptr,
			(n_processed + 4 < n_seqs) ? seq_to_ptr(sequences[n_processed + 4]) : nullptr,
			(n_processed + 5 < n_seqs) ? seq_to_ptr(sequences[n_processed + 5]) : nullptr,
			(n_processed + 6 < n_seqs) ? seq_to_ptr(sequences[n_processed + 6]) : nullptr,
			(n_processed + 7 < n_seqs) ? seq_to_ptr(sequences[n_processed + 7]) : nullptr,
			lcs_lens);

		for (int k = 0; k < 8 && n_processed + k < n_seqs; ++k)
			out_vector[n_processed + k] = transform(lcs_lens[k], seq_to_ptr(ref)->length, seq_to_ptr(sequences[n_processed + k])->length);
	}

	seq_to_ptr(ref)->ReleaseBitMasks();
}

// *******************************************************************
/*
	seq_type can be:
	- CSequence,
	- CSequence*,
*/
template <class seq_type, class distance_type, typename Iter, typename Transform>
void AbstractTreeGenerator::calculateDistanceRange(
	Transform &transform,
	seq_type& ref, 
	seq_type* sequences, 
	pair<Iter, Iter> ids_range,
	distance_type* out_vector,
	CLCSBP& lcsbp)
{
	uint32_t lcs_lens[8];
	int n_seqs = distance(ids_range.first, ids_range.second);

	auto p_ids = ids_range.first;

	if (n_seqs >= 16)
	{
		tpl_prefetch(seq_to_ptr(sequences[*(p_ids + 8)]));
		tpl_prefetch(seq_to_ptr(sequences[*(p_ids + 9)]));
		tpl_prefetch(seq_to_ptr(sequences[*(p_ids + 10)]));
		tpl_prefetch(seq_to_ptr(sequences[*(p_ids + 11)]));
		tpl_prefetch(seq_to_ptr(sequences[*(p_ids + 12)]));
		tpl_prefetch(seq_to_ptr(sequences[*(p_ids + 13)]));
		tpl_prefetch(seq_to_ptr(sequences[*(p_ids + 14)]));
		tpl_prefetch(seq_to_ptr(sequences[*(p_ids + 15)]));
	}

	// process portions of 8 sequences
	for (int j = 0; j < n_seqs / 8; ++j, p_ids += 8) {
		if (j + 3 < n_seqs / 8)
		{
			tpl_prefetch(seq_to_ptr(sequences[*(p_ids + 16)]));
			tpl_prefetch(seq_to_ptr(sequences[*(p_ids + 17)]));
			tpl_prefetch(seq_to_ptr(sequences[*(p_ids + 18)]));
			tpl_prefetch(seq_to_ptr(sequences[*(p_ids + 19)]));
			tpl_prefetch(seq_to_ptr(sequences[*(p_ids + 20)]));
			tpl_prefetch(seq_to_ptr(sequences[*(p_ids + 21)]));
			tpl_prefetch(seq_to_ptr(sequences[*(p_ids + 22)]));
			tpl_prefetch(seq_to_ptr(sequences[*(p_ids + 23)]));
		}

		if (j + 2 < n_seqs / 8)
		{
			tpl_prefetch<symbol_t>(seq_to_ptr(sequences[*(p_ids + 8)])->data);
			tpl_prefetch<symbol_t>(seq_to_ptr(sequences[*(p_ids + 9)])->data);
			tpl_prefetch<symbol_t>(seq_to_ptr(sequences[*(p_ids + 10)])->data);
			tpl_prefetch<symbol_t>(seq_to_ptr(sequences[*(p_ids + 11)])->data);
			tpl_prefetch<symbol_t>(seq_to_ptr(sequences[*(p_ids + 12)])->data);
			tpl_prefetch<symbol_t>(seq_to_ptr(sequences[*(p_ids + 13)])->data);
			tpl_prefetch<symbol_t>(seq_to_ptr(sequences[*(p_ids + 14)])->data);
			tpl_prefetch<symbol_t>(seq_to_ptr(sequences[*(p_ids + 15)])->data);
		}

		lcsbp.GetLCSBP(
			seq_to_ptr(ref),
			seq_to_ptr(sequences[*(p_ids + 0)]),
			seq_to_ptr(sequences[*(p_ids + 1)]),
			seq_to_ptr(sequences[*(p_ids + 2)]),
			seq_to_ptr(sequences[*(p_ids + 3)]),
			seq_to_ptr(sequences[*(p_ids + 4)]),
			seq_to_ptr(sequences[*(p_ids + 5)]),
			seq_to_ptr(sequences[*(p_ids + 6)]),
			seq_to_ptr(sequences[*(p_ids + 7)]),
			lcs_lens);

		auto ref_length = seq_to_ptr(ref)->length;

		for (int k = 0; k < 8; ++k) {
			out_vector[j * 8 + k] = transform(lcs_lens[k], ref_length, seq_to_ptr(sequences[*(p_ids + k)])->length);
		}
	}

	// if there is something left
	int n_processed = n_seqs / 8 * 8;
	if (n_processed < n_seqs) {
		lcsbp.GetLCSBP(
			seq_to_ptr(ref),
			seq_to_ptr(sequences[*(p_ids + 0)]),
			(n_processed + 1 < n_seqs) ? seq_to_ptr(sequences[*(p_ids + 1)]) : seq_to_ptr(sequences[*(p_ids + 0)]),
			(n_processed + 2 < n_seqs) ? seq_to_ptr(sequences[*(p_ids + 2)]) : seq_to_ptr(sequences[*(p_ids + 0)]),
			(n_processed + 3 < n_seqs) ? seq_to_ptr(sequences[*(p_ids + 3)]) : seq_to_ptr(sequences[*(p_ids + 0)]),
			(n_processed + 4 < n_seqs) ? seq_to_ptr(sequences[*(p_ids + 4)]) : seq_to_ptr(sequences[*(p_ids + 0)]),
			(n_processed + 5 < n_seqs) ? seq_to_ptr(sequences[*(p_ids + 5)]) : seq_to_ptr(sequences[*(p_ids + 0)]),
			(n_processed + 6 < n_seqs) ? seq_to_ptr(sequences[*(p_ids + 6)]) : seq_to_ptr(sequences[*(p_ids + 0)]),
			(n_processed + 7 < n_seqs) ? seq_to_ptr(sequences[*(p_ids + 7)]) : seq_to_ptr(sequences[*(p_ids + 0)]),
			lcs_lens);

		for (int k = 0; k < 8 && n_processed + k < n_seqs; ++k)
			out_vector[n_processed + k] = transform(lcs_lens[k], seq_to_ptr(ref)->length, seq_to_ptr(sequences[*(p_ids + k)])->length);
	}
}

// *******************************************************************
/*
	seq_type can be:
	- CSequenceView,
	- CSequenceView*,
*/
template <class seq_type, class sv_type, class distance_type, typename Iter, typename Transform>
void AbstractTreeGenerator::calculateDistanceRangeSV(
	Transform& transform,
	seq_type& ref,
	sv_type* sv,
	pair<Iter, Iter> ids_range,
	distance_type* out_vector,
	CLCSBP& lcsbp)
{
	uint32_t lcs_lens[8];
	int n_seqs = distance(ids_range.first, ids_range.second);

	auto p_ids = ids_range.first;

	if (n_seqs >= 16)
	{
		tpl_prefetch(sv_to_ptr(sv[*(p_ids + 8)]));
		tpl_prefetch(sv_to_ptr(sv[*(p_ids + 9)]));
		tpl_prefetch(sv_to_ptr(sv[*(p_ids + 10)]));
		tpl_prefetch(sv_to_ptr(sv[*(p_ids + 11)]));
		tpl_prefetch(sv_to_ptr(sv[*(p_ids + 12)]));
		tpl_prefetch(sv_to_ptr(sv[*(p_ids + 13)]));
		tpl_prefetch(sv_to_ptr(sv[*(p_ids + 14)]));
		tpl_prefetch(sv_to_ptr(sv[*(p_ids + 15)]));
	}

	// process portions of 8 sequences
	for (int j = 0; j < n_seqs / 8; ++j, p_ids += 8) {
		if (j + 3 < n_seqs / 8)
		{
			tpl_prefetch(sv_to_ptr(sv[*(p_ids + 16)]));
			tpl_prefetch(sv_to_ptr(sv[*(p_ids + 17)]));
			tpl_prefetch(sv_to_ptr(sv[*(p_ids + 18)]));
			tpl_prefetch(sv_to_ptr(sv[*(p_ids + 19)]));
			tpl_prefetch(sv_to_ptr(sv[*(p_ids + 20)]));
			tpl_prefetch(sv_to_ptr(sv[*(p_ids + 21)]));
			tpl_prefetch(sv_to_ptr(sv[*(p_ids + 22)]));
			tpl_prefetch(sv_to_ptr(sv[*(p_ids + 23)]));
		}

		if (j + 2 < n_seqs / 8)
		{
			tpl_prefetch<symbol_t>(sv_to_ptr(sv[*(p_ids + 8)])->data);
			tpl_prefetch<symbol_t>(sv_to_ptr(sv[*(p_ids + 9)])->data);
			tpl_prefetch<symbol_t>(sv_to_ptr(sv[*(p_ids + 10)])->data);
			tpl_prefetch<symbol_t>(sv_to_ptr(sv[*(p_ids + 11)])->data);
			tpl_prefetch<symbol_t>(sv_to_ptr(sv[*(p_ids + 12)])->data);
			tpl_prefetch<symbol_t>(sv_to_ptr(sv[*(p_ids + 13)])->data);
			tpl_prefetch<symbol_t>(sv_to_ptr(sv[*(p_ids + 14)])->data);
			tpl_prefetch<symbol_t>(sv_to_ptr(sv[*(p_ids + 15)])->data);
		}

		lcsbp.GetLCSBP(
			seq_to_ptr(ref),
			sv_to_ptr(sv[*(p_ids + 0)]),
			sv_to_ptr(sv[*(p_ids + 1)]),
			sv_to_ptr(sv[*(p_ids + 2)]),
			sv_to_ptr(sv[*(p_ids + 3)]),
			sv_to_ptr(sv[*(p_ids + 4)]),
			sv_to_ptr(sv[*(p_ids + 5)]),
			sv_to_ptr(sv[*(p_ids + 6)]),
			sv_to_ptr(sv[*(p_ids + 7)]),
			lcs_lens);

		auto ref_length = seq_to_ptr(ref)->length;

		for (int k = 0; k < 8; ++k)
			out_vector[j * 8 + k] = transform(lcs_lens[k], ref_length, sv_to_ptr(sv[*(p_ids + k)])->length);
	}

	// if there is something left
	int n_processed = n_seqs / 8 * 8;
	if (n_processed < n_seqs) {
		lcsbp.GetLCSBP(
			seq_to_ptr(ref),
			sv_to_ptr(sv[*(p_ids + 0)]),
			(n_processed + 1 < n_seqs) ? sv_to_ptr(sv[*(p_ids + 1)]) : sv_to_ptr(sv[*(p_ids + 0)]),
			(n_processed + 2 < n_seqs) ? sv_to_ptr(sv[*(p_ids + 2)]) : sv_to_ptr(sv[*(p_ids + 0)]),
			(n_processed + 3 < n_seqs) ? sv_to_ptr(sv[*(p_ids + 3)]) : sv_to_ptr(sv[*(p_ids + 0)]),
			(n_processed + 4 < n_seqs) ? sv_to_ptr(sv[*(p_ids + 4)]) : sv_to_ptr(sv[*(p_ids + 0)]),
			(n_processed + 5 < n_seqs) ? sv_to_ptr(sv[*(p_ids + 5)]) : sv_to_ptr(sv[*(p_ids + 0)]),
			(n_processed + 6 < n_seqs) ? sv_to_ptr(sv[*(p_ids + 6)]) : sv_to_ptr(sv[*(p_ids + 0)]),
			(n_processed + 7 < n_seqs) ? sv_to_ptr(sv[*(p_ids + 7)]) : sv_to_ptr(sv[*(p_ids + 0)]),
			lcs_lens);

		for (int k = 0; k < 8 && n_processed + k < n_seqs; ++k)
			out_vector[n_processed + k] = transform(lcs_lens[k], seq_to_ptr(ref)->length, sv_to_ptr(sv[*(p_ids + k)])->length);
	}
}

// *******************************************************************
template <class seq_type, class distance_type, typename Transform>
void AbstractTreeGenerator::calculateDistanceMatrix(
	Transform& transform,
	seq_type* sequences,
	int n_seq, 
	distance_type *out_matrix, 
	CLCSBP& lcsbp) {

	for (int row_id = 0; row_id < n_seq; ++row_id) {
		
		size_t row_offset = TriangleMatrix::access(row_id, 0);
		
		calculateDistanceVector<seq_type, distance_type, decltype(transform)>(
			transform,
			sequences[row_id],
			sequences,
			row_id,
			out_matrix + row_offset,
			lcsbp);
	}	
}