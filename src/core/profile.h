/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#ifndef _PROFILE_H
#define _PROFILE_H

#include "../core/defs.h"
#include "../core/params.h"
#include "../utils/utils.h"

#include <vector>
#include <tuple>
#include <array>
#include <algorithm>
#include <cstring>

#ifdef WIN32
#include <xmmintrin.h>
#include <mmintrin.h>
#endif


using namespace std;

class CGappedSequence;

enum class direction_t { D, H, V };

// *********************************************************************************
class CDPMatrix {
	size_t n_rows;
	size_t n_cols;

	unsigned char* raw_data;

public:
	CDPMatrix(size_t _n_rows, size_t _n_cols) : n_rows(_n_rows), n_cols(_n_cols)
	{
		raw_data = new unsigned char[n_rows*n_cols];
	}

	~CDPMatrix()
	{
		delete[] raw_data;
	}

	void set_zeros(instruction_set_t instruction_set = instruction_set_t::none)
	{
#if SIMD==SIMD_NONE
		mem_clear(raw_data, n_rows * n_cols);
#elif SIMD==SIMD_AVX1
		if (instruction_set < instruction_set_t::avx) {
			mem_clear(raw_data, n_rows * n_cols);
		}
		else {
			mem_clear_avx(raw_data, n_rows * n_cols);
		}
#elif SIMD==SIMD_AVX2 || SIMD==SIMD_AVX512
		if (instruction_set < instruction_set_t::avx) {
			mem_clear(raw_data, n_rows * n_cols);
		}
		else if (instruction_set < instruction_set_t::avx2) {
			mem_clear_avx(raw_data, n_rows * n_cols);
		}
		else {
			mem_clear_avx2(raw_data, n_rows * n_cols);
		}		
#elif SIMD==SIMD_NEON
		mem_clear_neon(raw_data, n_rows * n_cols);
#else
		// Impossible to be here
#endif
	}

	unsigned char *get_row(size_t row_id)
	{
		return raw_data + row_id * n_cols;
	}

	unsigned char *get_cell(size_t row_id, size_t col_id)
	{
		return raw_data + row_id * n_cols + col_id;
	}

	direction_t get_dir_D(size_t row_id, size_t col_id)
	{
		return (direction_t) (raw_data[row_id * n_cols + col_id] & 0x03);
	}

	direction_t get_dir_H(size_t row_id, size_t col_id)
	{
		return (direction_t) ((raw_data[row_id * n_cols + col_id] >> 2) & 0x03);
	}

	direction_t get_dir_V(size_t row_id, size_t col_id)
	{
		return (direction_t) ((raw_data[row_id * n_cols + col_id] >> 4) & 0x03);
	}

	void set_dir_D(size_t row_id, size_t col_id, direction_t dir)
	{
		auto p = raw_data + row_id * n_cols + col_id;
		*p = (*p & (0xff - 0x03)) + (unsigned char) dir;
	}

	// Assumption: this is the first attempt to store in this field
	void set_dir_D(unsigned char *cell, direction_t dir)
	{
		*cell += (unsigned char) dir;
	}

	void set_dir_H(size_t row_id, size_t col_id, direction_t dir)
	{
		auto p = raw_data + row_id * n_cols + col_id;
		*p = (*p & (0xff - (0x03 << 2))) + (((unsigned char)dir) << 2);
	}

	// Assumption: this is the first attempt to store in this field
	void set_dir_H(unsigned char *cell, direction_t dir)
	{
		*cell += ((unsigned char) dir) << 2;
	}

	void set_dir_V(size_t row_id, size_t col_id, direction_t dir)
	{
		auto p = raw_data + row_id * n_cols + col_id;
		*p = (*p & (0xff - (0x03 << 4))) + (((unsigned char)dir) << 4);
	}

	// Assumption: this is the first attempt to store in this field
	void set_dir_V(unsigned char *cell, direction_t dir)
	{
		*cell += ((unsigned char) dir) << 4;
	}

	void set_dir_all(size_t row_id, size_t col_id, direction_t dir)
	{
		unsigned char x = (unsigned char) dir;
		raw_data[row_id * n_cols + col_id] = x + (x << 2) + (x << 4);
	}
};

// *********************************************************************************
// Elements are organized in columns (not in rows)!
template<typename T, unsigned N_ROWS> class CProfileValues {
	size_t n_cols;
	size_t allocated_n_cols;

	T* raw_data;

	size_t round_n_cols_imp(size_t n, size_t r)
	{
		return (n + r - 1) / r * r;
	}

	size_t round_n_cols(size_t n)
	{
		return n;

		// For possible future use when extending of profiles will be implemented
/*		if (n < 128)
			return round_n_cols_imp(n, 16);
		if (n < 256)
			return round_n_cols_imp(n, 16);
		if (n < 512)
			return round_n_cols_imp(n, 32);
		if (n < 1024)
			return round_n_cols_imp(n, 64);
		if (n < 2 * 1024)
			return round_n_cols_imp(n, 128);
		if (n < 4 * 1024)
			return round_n_cols_imp(n, 256);
		if (n < 8 * 1024)
			return round_n_cols_imp(n, 512);
		if (n < 16 * 1024)
			return round_n_cols_imp(n, 1024);
		if (n < 32 * 1024)
			return round_n_cols_imp(n, 2 * 1024);
		if (n < 64 * 1024)
			return round_n_cols_imp(n, 4 * 1024);

		return round_n_cols_imp(n, 8 * 1024);*/
	}

	void allocate(size_t _n_cols)
	{
		n_cols = _n_cols;

		allocated_n_cols = round_n_cols(n_cols);

		if(raw_data)
		{
			delete[] raw_data;
		}

		if(!n_cols)
		{
			raw_data = nullptr;
			return;
		}

		raw_data = new T[N_ROWS * allocated_n_cols];
	}

	void deallocate()
	{
		if(raw_data)
		{
			delete[] raw_data;
			
			raw_data = nullptr;
			
			n_cols = 0;
			allocated_n_cols = 0;
		}
	}

public:
	CProfileValues(size_t _n_cols = 0)
	{
		raw_data = nullptr;

		allocate(_n_cols);
	}

	CProfileValues(const CProfileValues &x)
	{
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
		std::swap(raw_data, x.raw_data);
		std::swap(n_cols, x.n_cols);
		std::swap(allocated_n_cols, x.allocated_n_cols);
	}

	bool empty()
	{
		return n_cols == 0;
	}

	void resize(size_t _n_cols)
	{
		if (round_n_cols(_n_cols) != allocated_n_cols)
		{
			deallocate();
			allocate(_n_cols);
		}
		else
			n_cols = _n_cols;
	}

	void clear(void)
	{
		deallocate();
	}

	void set_zeros(instruction_set_t instruction_set = instruction_set_t::none)
	{
#if SIMD==SIMD_NONE
		memset(raw_data, 0, N_ROWS * n_cols * sizeof(T));
#elif SIMD==SIMD_AVX1
		if (instruction_set < instruction_set_t::avx) {
			memset(raw_data, 0, N_ROWS * n_cols * sizeof(T));
		}
		else {
			mem_clear_avx(raw_data, N_ROWS * n_cols * sizeof(T));
		}
#elif SIMD==SIMD_AVX2 || SIMD==SIMD_AVX512
		if (instruction_set < instruction_set_t::avx) {
			memset(raw_data, 0, N_ROWS * n_cols * sizeof(T));
		}
		else if (instruction_set < instruction_set_t::avx2) {
			mem_clear_avx(raw_data, N_ROWS * n_cols * sizeof(T));
		}
		else {
			mem_clear_avx2(raw_data, N_ROWS * n_cols * sizeof(T));
		}
#elif SIMD==SIMD_NEON
		mem_clear_neon(raw_data, N_ROWS * n_cols * sizeof(T));
#endif			
	}

	size_t get_num_of_non_zeros(void)
	{
		return n_cols * N_ROWS - count(raw_data, raw_data + n_cols * N_ROWS, T());
	}

	T get_value(size_t col_id, size_t row_id)
	{
		return raw_data[col_id * N_ROWS + row_id];		// elements are organized in columns!
	}

	void set_value(size_t col_id, size_t row_id, T value)
	{
		raw_data[col_id * N_ROWS + row_id] = value;
	}

	void add_value(size_t col_id, size_t row_id, T value)
	{
 		raw_data[col_id * N_ROWS + row_id] += value;
	}

	void sub_value(size_t col_id, size_t row_id, T value)
	{
		raw_data[col_id * N_ROWS + row_id] -= value;
	}

	void add_column(size_t col_id, T* source)
	{
		T* dest = &raw_data[col_id * N_ROWS];
		
		for (size_t i = 0; i < N_ROWS; ++i)
			*dest++ += *source++;
	}

	void add_value_to_column_part(size_t col_id, size_t max_row_id, T value)
	{
		T* dest = &raw_data[col_id * N_ROWS];

		for(size_t i = 0; i < max_row_id; ++i)
 			*dest++ += value;
	}

	void add_column_part(size_t col_id, size_t n, vector<score_t> &source)
	{
		T* dest = &raw_data[col_id * N_ROWS];
		
		for(size_t i = 0; i < n; ++i)
			dest[i] += source[i];
	}

	void add_column_part_mult(size_t col_id, size_t n, vector<score_t> &source, score_t mult)
	{
		T* dest = &raw_data[col_id * N_ROWS];
		
		for(size_t i = 0; i < n; ++i)
			dest[i] += source[i] * mult;
	}

	T* get_column(size_t col_id)
	{
		return &raw_data[col_id * N_ROWS];
	}

	void prefetch_column(size_t col_id)
	{
#ifdef WIN32
		_mm_prefetch((const char*)(raw_data + N_ROWS * col_id), _MM_HINT_T0);
#else
		__builtin_prefetch(raw_data * N_ROWS + col_id);
#endif
	}

	void set_column(size_t col_id, T* column)
	{
		copy_n(column, N_ROWS, &raw_data[col_id * N_ROWS]);
	}
};

// *********************************************************************************
class CProfile {
	CParams *params;

	struct dp_row_elem_t {
		score_t D, H, V;

		dp_row_elem_t(score_t _D = 0.0, score_t _H = 0.0, score_t _V = 0.0) :
			D(_D), H(_H), V(_V)
		{};
	};

	typedef vector<dp_row_elem_t> dp_row_t;

	struct dp_row_gap_info_elem {
		counter_t n_gap_open_at_left, n_gap_ext_at_left, n_gap_term_open_at_left, n_gap_term_ext_at_left;

		dp_row_gap_info_elem(counter_t _n_gap_open_at_left = 0, counter_t _n_gap_ext_at_left = 0,
			counter_t _n_gap_term_open_at_left = 0, counter_t _n_gap_term_ext_at_left = 0) :
			n_gap_open_at_left(_n_gap_open_at_left), n_gap_ext_at_left(_n_gap_ext_at_left),
			n_gap_term_open_at_left(_n_gap_term_open_at_left), n_gap_term_ext_at_left(_n_gap_term_ext_at_left)
		{};
	};

	struct dp_gap_corrections {
		counter_t n_gap_start_open, n_gap_start_ext, n_gap_start_term_open, n_gap_start_term_ext;
		counter_t n_gap_cont_ext, n_gap_cont_term_ext;

		dp_gap_corrections() : n_gap_start_open(0), n_gap_start_ext(0), n_gap_start_term_open(0), n_gap_start_term_ext(0),
			n_gap_cont_ext(0), n_gap_cont_term_ext(0)
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
	
#ifndef NO_PROFILE_PAR
	// parallel variant of profile alignment
	void ParAlignSeqProf(CProfile* profile1, CProfile* profile2, uint32_t no_threads, uint32_t rows_per_box);
	void ParAlignProfProf(CProfile* profile1, CProfile* profile2, uint32_t no_threads, uint32_t rows_per_box);
#endif

	void ConstructProfile(CProfile *profile1, CProfile *profile2, CDPMatrix &matrix, dp_row_elem_t &last_row, uint32_t no_threads = 1);
	
	void InsertGaps(size_t prof_col_id, CProfile *profile, size_t col_id, 
		counter_t n_gap_open, counter_t n_gap_ext, counter_t n_gap_term_open, counter_t n_gap_term_ext, vector<pair<uint32_t, uint32_t>> &v_gaps_prof);
	void InsertColumn(size_t prof_col_id, CProfile *profile, size_t col_id);
	void FinalizeGaps(CProfile* profile, vector<pair<uint32_t, uint32_t>>& v_gaps_prof, uint32_t no_threads = 1);

	void CalculateCounters(CGappedSequence *gs);
    void CalculateScores();

	void SolveGapsProblemWhenContinuing(size_t source_col_id, size_t prof_width, size_t prof_size, 
		counter_t&n_gap_to_transfer, counter_t&n_gap_term_to_transfer, counter_t&n_gap_open, counter_t&n_gap_ext, counter_t&n_gap_term_open, counter_t&n_gap_term_ext, 
		counter_t n_gap_open_at_left, counter_t n_gap_ext_at_left, counter_t n_gap_term_open_at_left, counter_t n_gap_term_ext_at_left);
	void SolveGapsProblemWhenStarting(size_t source_col_id, size_t prof_width, size_t prof_size, CProfile *profile, 
		counter_t&n_gap_to_transfer, counter_t&n_gap_term_to_transfer, counter_t&n_gap_open, counter_t&n_gap_ext, counter_t&n_gap_term_open, counter_t&n_gap_term_ext);
	void DP_SolveGapsProblemWhenStarting(size_t source_col_id, size_t prof_width, size_t prof_size, CProfile *profile, 
		counter_t&n_gap_open, counter_t&n_gap_ext, counter_t&n_gap_term_open, counter_t&n_gap_term_ext);
	void DP_SolveGapsProblemWhenContinuing(size_t source_col_id, size_t prof_width, size_t prof_size, CProfile *profile, 
		counter_t&n_gap_ext, counter_t&n_gap_term_ext);

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
	CProfile(CProfile *profile1, CProfile *profile2, CParams *_params, uint32_t no_threads, uint32_t no_rows_per_box);
	~CProfile();

	bool operator==(const CProfile &profile) const;

	void Align(CProfile *profile1, CProfile *profile2, uint32_t no_threads, uint32_t no_rows_per_box, vector<int> *column_mapping1 = nullptr, vector<int> *column_mapping2 = nullptr);

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