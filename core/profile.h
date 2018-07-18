/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#ifndef _PROFILE_H
#define _PROFILE_H

#include "../core/defs.h"
#include "../core/sequence.h"
#include "../core/params.h"

#include <vector>
#include <tuple>
#include <array>
#include <algorithm>

#include "../libs/asmlib.h"

using namespace std;

enum class direction_t { D, H, V };

// *********************************************************************************
class CDPMatrix {
	size_t n_rows;
	size_t n_cols;

	unsigned char* raw_data;
	unsigned char** data;

public:
	CDPMatrix(size_t _n_rows, size_t _n_cols) : n_rows(_n_rows), n_cols(_n_cols)
	{
		raw_data = new unsigned char[n_rows*n_cols];
		data = new unsigned char*[n_rows];
		for(size_t i = 0; i < n_rows; ++i)
			data[i] = raw_data + i*n_cols;
	}

	~CDPMatrix()
	{
		delete[] raw_data;
		delete[] data;
	}

	void set_zeros(void)
	{
		A_memset(raw_data, 0, n_rows * n_cols);
	}

	unsigned char *get_row(size_t row_id)
	{
		return data[row_id];
	}

	unsigned char *get_cell(size_t row_id, size_t col_id)
	{
		return &data[row_id][col_id];
	}

	direction_t get_dir_D(size_t row_id, size_t col_id)
	{
		return (direction_t) (data[row_id][col_id] & 0x03);
	}

	direction_t get_dir_H(size_t row_id, size_t col_id)
	{
		return (direction_t) ((data[row_id][col_id] >> 2) & 0x03);
	}

	direction_t get_dir_V(size_t row_id, size_t col_id)
	{
		return (direction_t) ((data[row_id][col_id] >> 4) & 0x03);
	}

	void set_dir_D(size_t row_id, size_t col_id, direction_t dir)
	{
		data[row_id][col_id] = (data[row_id][col_id] & (0xff - 0x03)) + (unsigned char) dir;
	}

	// Assumption: this is the first attempt to store in this field
	void set_dir_D(unsigned char *cell, direction_t dir)
	{
		*cell += (unsigned char) dir;
	}

	void set_dir_H(size_t row_id, size_t col_id, direction_t dir)
	{
		data[row_id][col_id] = (data[row_id][col_id] & (0xff - (0x03 << 2))) + (((unsigned char) dir) << 2);
	}

	// Assumption: this is the first attempt to store in this field
	void set_dir_H(unsigned char *cell, direction_t dir)
	{
		*cell += ((unsigned char) dir) << 2;
	}

	void set_dir_V(size_t row_id, size_t col_id, direction_t dir)
	{
		data[row_id][col_id] = (data[row_id][col_id] & (0xff - (0x03 << 4))) + (((unsigned char) dir) << 4);
	}

	// Assumption: this is the first attempt to store in this field
	void set_dir_V(unsigned char *cell, direction_t dir)
	{
		*cell += ((unsigned char) dir) << 4;
	}

	void set_dir_all(size_t row_id, size_t col_id, direction_t dir)
	{
		unsigned char x = (unsigned char) dir;
		data[row_id][col_id] = x + (x << 2) + (x << 4);
	}
};

// *********************************************************************************
// Elements are organized in columns (not in rows)!
template<typename T, unsigned N_ROWS> class CProfileValues {
	size_t n_cols;

	T* raw_data;
	T** data;

	void allocate(size_t _n_cols)
	{
		n_cols = _n_cols;

		if(data)
		{
			delete[] raw_data;
			delete[] data;
		}

		if(!n_cols)
		{
			raw_data = nullptr;
			data     = nullptr;
			return;
		}

		raw_data = new T[N_ROWS * n_cols];
		data = new T*[n_cols];

		for(size_t i = 0; i < n_cols; ++i)
			data[i] = raw_data + i * N_ROWS;
	}

	void deallocate()
	{
		if(data)
		{
			delete[] raw_data;
			delete[] data;
			
			raw_data = nullptr;
			data     = nullptr;
			
			n_cols = 0;
		}
	}

public:
	CProfileValues(size_t _n_cols = 0)
	{
		data = nullptr;
		raw_data = nullptr;

		allocate(_n_cols);
	}

	CProfileValues(const CProfileValues &x)
	{
		data = nullptr;
		raw_data = nullptr;

		allocate(x.n_cols);

		copy(x.raw_data, x.raw_data + N_ROWS * n_cols, raw_data);
	}

	~CProfileValues()
	{
		deallocate();
	}

	CProfileValues& operator=(const CProfileValues &x)
	{
		allocate(x.n_cols);

		copy(x.raw_data, x.raw_data + N_ROWS * n_cols, raw_data);

		return *this;
	}

	void swap(CProfileValues<T, N_ROWS> &x)
	{
		std::swap(data, x.data);
		std::swap(raw_data, x.raw_data);
		std::swap(n_cols, x.n_cols);
	}

	bool empty()
	{
		return n_cols == 0;
	}

	void resize(size_t _n_cols)
	{
		deallocate();
		allocate(_n_cols);
	}

	void clear(void)
	{
		deallocate();
	}

	void set_zeros(void)
	{
		fill(raw_data, raw_data + N_ROWS * n_cols, T());
	}

	size_t get_num_of_non_zeros(void)
	{
		return n_cols * N_ROWS - count(raw_data, raw_data + n_cols * N_ROWS, T());
	}

	T get_value(size_t col_id, size_t row_id)
	{
		return data[col_id][row_id];		// elements are organized in columns!
	}

	void set_value(size_t col_id, size_t row_id, T value)
	{
		data[col_id][row_id] = value;
	}

	void add_value(size_t col_id, size_t row_id, T value)
	{
 		data[col_id][row_id] += value;
	}

	void sub_value(size_t col_id, size_t row_id, T value)
	{
		data[col_id][row_id] -= value;
	}

	void add_column(size_t col_id, T* source)
	{
		T* dest = data[col_id];
		
		for(size_t i = 0; i < N_ROWS; ++i)
			dest[i] += source[i];
	}

	void add_value_to_column_part(size_t col_id, size_t max_row_id, T value)
	{
		for(size_t i = 0; i < max_row_id; ++i)
 			data[col_id][i] += value;
	}

	void add_column_part(size_t col_id, size_t n, vector<score_t> &source)
	{
		T* dest = data[col_id];
		
		for(size_t i = 0; i < n; ++i)
			dest[i] += source[i];
	}

	void add_column_part_mult(size_t col_id, size_t n, vector<score_t> &source, score_t mult)
	{
		T* dest = data[col_id];
		
		for(size_t i = 0; i < n; ++i)
			dest[i] += source[i] * mult;
	}

	T* get_column(size_t col_id)
	{
		return data[col_id];
	}

	void set_column(size_t col_id, T* column)
	{
		copy(column, column + N_ROWS, data[col_id]);
	}
};

// *********************************************************************************
class CProfile {
	CParams *params;

	struct dp_matrix_elem_t {
		direction_t dir_D, dir_H, dir_V;

		dp_matrix_elem_t()
		{};
	};

	struct dp_row_elem_t {
		score_t D, H, V;

		dp_row_elem_t(score_t _D = 0.0, score_t _H = 0.0, score_t _V = 0.0) :
			D(_D), H(_H), V(_V)
		{};
	};

	typedef vector<dp_row_elem_t> dp_row_t;

	struct dp_row_gap_info_elem {
		size_t n_gap_open_at_left, n_gap_ext_at_left, n_gap_term_open_at_left, n_gap_term_ext_at_left;

		dp_row_gap_info_elem(size_t _n_gap_open_at_left = 0, size_t _n_gap_ext_at_left = 0,
			size_t _n_gap_term_open_at_left = 0, size_t _n_gap_term_ext_at_left = 0) :
			n_gap_open_at_left(_n_gap_open_at_left), n_gap_ext_at_left(_n_gap_ext_at_left),
			n_gap_term_open_at_left(_n_gap_term_open_at_left), n_gap_term_ext_at_left(_n_gap_term_ext_at_left)
		{};
	};

	struct dp_gap_corrections {
		size_t n_gap_start_open, n_gap_start_ext, n_gap_start_term_open, n_gap_start_term_ext;
		size_t n_gap_cont_ext, n_gap_cont_term_ext;

		dp_gap_corrections() : n_gap_start_open(0), n_gap_start_ext(0), n_gap_start_term_open(0), n_gap_start_term_ext(0),
			n_gap_cont_ext(0), n_gap_cont_term_ext()
		{};
	};

	struct dp_gap_costs {
		score_t open, ext, term_open, term_ext;
	};

	typedef vector<dp_row_gap_info_elem> dp_row_gap_info_t;

	bool cumulate_gap_inserts;
	int32_t no_cumulated_gap_inserts;

	void FindRowRanges(vector<int> *column_mapping1, vector<int> *column_mapping2, vector<pair<int, int>> &row_ranges);

	void AlignSeqSeq(CProfile *profile1, CProfile *profile2);
	void AlignSeqProf(CProfile *profile1, CProfile *profile2, vector<int> *column_mapping1, vector<int> *column_mapping2);
	void AlignProfProf(CProfile *profile1, CProfile *profile2, vector<int> *column_mapping1, vector<int> *column_mapping2);
	
	
	void ConstructProfile(CProfile *profile1, CProfile *profile2, CDPMatrix &matrix, dp_row_elem_t &last_row);
	
	void InsertGaps(size_t prof_col_id, CProfile *profile, size_t col_id, size_t n_gap_open, size_t n_gap_ext, size_t n_gap_term_open, size_t n_gap_term_ext);
	void InsertColumn(size_t prof_col_id, CProfile *profile, size_t col_id);

	void CalculateCounters(CGappedSequence *gs);
    void CalculateScores();

	void SolveGapsProblemWhenContinuing(size_t source_col_id, size_t prof_width, size_t prof_size, size_t &n_gap_to_transfer, size_t &n_gap_term_to_transfer, size_t &n_gap_open, size_t &n_gap_ext, size_t &n_gap_term_open, size_t &n_gap_term_ext, size_t n_gap_open_at_left, size_t n_gap_ext_at_left, size_t n_gap_term_open_at_left, size_t n_gap_term_ext_at_left);
	void SolveGapsProblemWhenStarting(size_t source_col_id, size_t prof_width, size_t prof_size, CProfile *profile, size_t &n_gap_to_transfer, size_t &n_gap_term_to_transfer, size_t &n_gap_open, size_t &n_gap_ext, size_t &n_gap_term_open, size_t &n_gap_term_ext);
	void DP_SolveGapsProblemWhenStarting(size_t source_col_id, size_t prof_width, size_t prof_size, CProfile *profile, size_t &n_gap_open, size_t &n_gap_ext, size_t &n_gap_term_open, size_t &n_gap_term_ext);
	void DP_SolveGapsProblemWhenContinuing(size_t source_col_id, size_t prof_width, size_t prof_size, CProfile *profile, size_t &n_gap_ext, size_t &n_gap_term_ext);

public:
    vector<CGappedSequence*> data;
	CProfileValues<score_t, NO_SYMBOLS> scores;
	CProfileValues<counter_t, NO_SYMBOLS> counters;

	size_t width;
	score_t total_score;

public:
	CProfile(CParams *_params = nullptr);
	CProfile(const CGappedSequence &gapped_sequence, CParams *_params);
	CProfile(const CProfile &profile);
	CProfile(CProfile *profile1, CProfile *profile2, CParams *_params);
	~CProfile();

	bool operator==(const CProfile &profile);
	bool operator!=(const CProfile &profile);

	void Align(CProfile *profile1, CProfile *profile2, vector<int> *column_mapping1 = nullptr, vector<int> *column_mapping2 = nullptr);

	void AppendRawSequence(const CGappedSequence &gs);
	bool Condense(vector<int> &column_mapping);
	bool OptimizeGaps();

	void Clear();
	bool Empty();
	size_t Size();

	void CalculateCountersScores();

	void GetGapStats(vector<size_t> &stats);

	score_t CalculateTotalScore(void);

	void Swap(CProfile &profile);
};


#endif