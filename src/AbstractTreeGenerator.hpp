/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#pragma once
#include "AbstractTreeGenerator.h"

#include "lcsbp.h"

#include <cmath>
#include <type_traits>

// overloads for converting sequence type to pointer
inline CSequence* seq_to_ptr(CSequence* x) { return x; }
inline CSequence* seq_to_ptr(CSequence& x) { return &x; }

// dummy implementation
template <class T, Measure measure>
struct Transform {
	//static_assert(0, "Cannot use dummy implementation");

	T operator()(uint32_t lcs, uint32_t len1, uint32_t len2) {
		return 0;
	}
};

template <class T>
struct Transform<T, Measure::SimilarityDefault>{
	T operator()(uint32_t lcs, uint32_t len1, uint32_t len2) {
		T indel = len1 + len2 - 2 * lcs;
		return indel == 0 ? ((T)lcs * 1000) : ((T)lcs / indel);
	}
};

template <class T>
struct Transform<T, Measure::DistanceReciprocal> {
	T operator()(uint32_t lcs, uint32_t len1, uint32_t len2) {
		T indel = len1 + len2 - 2 * lcs;
		return (T)indel / lcs;
	}
};

template <class T>
struct Transform<T, Measure::DistanceInverse> {
	T operator()(uint32_t lcs, uint32_t len1, uint32_t len2) {
		T indel = len1 + len2 - 2 * lcs;
		return indel == 0 ? (-(T)lcs * 1000) : (-(T)lcs / indel);
	}
};

template <class T>
struct Transform<T, Measure::DistanceLCSByLength> {
	T operator()(uint32_t lcs, uint32_t len1, uint32_t len2) {
		return 1.0 - (T)lcs / std::min(len1, len2);
	}
};

template <class T>
struct Transform<T, Measure::DistanceLCSByLengthCorrected> {
	T operator()(uint32_t lcs, uint32_t len1, uint32_t len2) {
		// make len1 the longer
		if (len1 < len2) { std::swap(len1, len2);  }

		T d = 1.0 - (T)lcs / len2;

		// MAFFT distance correction
		auto correction = [](T x, T y) ->T { return y / x * 0.1 + 10000 / (x + 10000) + 0.01;  };
		d = d / correction(len1, len2);
		return d;
	}
};

// *******************************************************************
/*
similarity_type can be:
	- CSequence,
	- CSequence*,
*/
template <class seq_type, class similarity_type, Measure measure>
void AbstractTreeGenerator::calculateSimilarityVector(
	seq_type& ref,
	seq_type* sequences,
	size_t n_seqs,
	similarity_type* out_vector,
	CLCSBP& lcsbp)
{
	uint32_t lcs_lens[4];
	Transform<similarity_type, measure> transform;

	seq_to_ptr(ref)->ComputeBitMasks();

	// process portions of 4 sequences
	for (int j = 0; j < n_seqs / 4; ++j) {
		lcsbp.GetLCSBP(
			seq_to_ptr(ref),
			seq_to_ptr(sequences[j * 4 + 0]),
			seq_to_ptr(sequences[j * 4 + 1]),
			seq_to_ptr(sequences[j * 4 + 2]),
			seq_to_ptr(sequences[j * 4 + 3]),
			lcs_lens[0], lcs_lens[1], lcs_lens[2], lcs_lens[3]);

		for (int k = 0; k < 4; ++k) {
			out_vector[j * 4 + k] = transform(lcs_lens[k], seq_to_ptr(ref)->length, seq_to_ptr(sequences[j * 4 + k])->length);
		}
	}

	// if there is something left
	size_t n_processed = n_seqs / 4 * 4;
	if (n_processed < n_seqs) {
		lcsbp.GetLCSBP(
			seq_to_ptr(ref),
			(n_processed + 0 < n_seqs) ? seq_to_ptr(sequences[n_processed + 0]) : nullptr,
			(n_processed + 1 < n_seqs) ? seq_to_ptr(sequences[n_processed + 1]) : nullptr,
			(n_processed + 2 < n_seqs) ? seq_to_ptr(sequences[n_processed + 2]) : nullptr,
			(n_processed + 3 < n_seqs) ? seq_to_ptr(sequences[n_processed + 3]) : nullptr,
			lcs_lens[0], lcs_lens[1], lcs_lens[2], lcs_lens[3]);

		for (int k = 0; k < 4 && n_processed + k < n_seqs; ++k)
		{

			out_vector[n_processed + k] = transform(lcs_lens[k], seq_to_ptr(ref)->length, seq_to_ptr(sequences[n_processed + k])->length);
		}
	}

	seq_to_ptr(ref)->ReleaseBitMasks();

}


// *******************************************************************
template <class seq_type, class similarity_type, Measure measure>
void AbstractTreeGenerator::calculateSimilarityMatrix(
	seq_type* sequences,
	size_t n_seq,
	similarity_type *out_matrix,
	CLCSBP& lcsbp) {

	for (size_t row_id = 0; row_id < n_seq; ++row_id) {

		size_t row_offset = TriangleMatrix::access(row_id, 0);

		calculateSimilarityVector<seq_type, similarity_type, measure>(
			sequences[row_id],
			sequences,
			row_id,
			out_matrix + row_offset,
			lcsbp);
	}

}
