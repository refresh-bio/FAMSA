/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#ifndef NO_PROFILE_PAR

#include "../core/profile.h"
#include "../core/sequence.h"
#include "../core/queues.h"
#include <assert.h>
#include <array>
#include <set>
#include <math.h>
#include <algorithm>
#include <stdio.h>
#include <thread>
#include <future>

//#include "../utils/pooled_threads.h"

#define max3(x, y, z)	(max((x), max((y), (z))))

#define NINF(x)		((x) < -999 ? -999: (x))

#define USE_ASYNC

#ifdef _MSC_VER
#include <barrier>
#else
#ifdef USE_NATIVE_BARRIERS
#include <barrier>
#else
thread_local size_t __barrier_favorite_hash =
std::hash<std::thread::id>()(std::this_thread::get_id());
#include "../libs/atomic_wait/barrier"
#endif
#endif

// ****************************************************************************
void CProfile::ParAlignSeqProf(CProfile* profile1, CProfile* profile2, uint32_t no_threads, uint32_t rows_per_box)
{
	int no_dp_rows = 3 * rows_per_box;

	size_t prof1_width = profile1->width;
	size_t prof2_width = profile2->width;

	//size_t prof1_card = profile1->data.size();
	//size_t prof1_card = 1;
	size_t prof2_card = profile2->data.size();

	// The profile1 contains a single sequence so there are no gaps here
//	vector<symbol_t>& seq1 = const_cast<vector<symbol_t>&>(profile1->data[0]->symbols);
	symbol_t* seq1 = const_cast<symbol_t*>(profile1->data[0]->symbols);

	score_t gap_open = params->gap_open;
	score_t gap_ext = params->gap_ext;
	score_t gap_term_open = params->gap_term_open;
	score_t gap_term_ext = params->gap_term_ext;

	CDPMatrix matrix(prof1_width + 1, prof2_width + 1);
	vector<dp_row_t> dp_rows(no_dp_rows);

	//CProfileValues<score_t, NO_SYMBOLS>& scores1 = const_cast<CProfileValues<score_t, NO_SYMBOLS>&>(profile1->scores);
	CProfileValues<score_t, NO_SYMBOLS>& scores2 = const_cast<CProfileValues<score_t, NO_SYMBOLS>&>(profile2->scores);

	// Precompute scores for gaps in profile2
	vector<dp_gap_costs> prof2_gaps;
	vector<dp_gap_corrections> gap_corrections;

	// Precompute gap corrections in profile2
	vector<score_t> n_gaps_prof2_to_change;
	vector<score_t> n_gaps_prof2_term_to_change;
	vector<score_t> gaps_prof2_change;

	vector<pair<score_t, score_t>> v_gap_corr;

	if (no_threads > 4)
	{
		auto fut0 = async([&] {
			for (int i = 0; i < no_dp_rows; ++i)
				dp_rows[i].resize(prof2_width + 1);
			});

		auto fut1 = async([&] {
			matrix.set_zeros(params->instruction_set);
			});

		auto fut2 = async([&] {
			gap_corrections.resize(prof2_width + 1);
			});

		auto fut3 = async([&] {
			n_gaps_prof2_to_change.resize(prof2_width + 1);
			n_gaps_prof2_term_to_change.resize(prof2_width + 1);
			gaps_prof2_change.resize(prof2_width + 1);

			v_gap_corr.resize(prof2_width + 1);
			});

		prof2_gaps.resize(prof2_width + 1);

		for (size_t j = 0; j <= prof2_width; ++j)
		{
			prof2_gaps[j].open = scores2.get_value(j, GAP_OPEN);
			prof2_gaps[j].ext = scores2.get_value(j, GAP_EXT);
			prof2_gaps[j].term_open = scores2.get_value(j, GAP_TERM_OPEN);
			prof2_gaps[j].term_ext = scores2.get_value(j, GAP_TERM_EXT);
		}

		fut0.wait();
		fut1.wait();
		fut2.wait();
		fut3.wait();
	}	
	else if (no_threads == 4)
	{
		auto fut0 = async([&] {
			for (int i = 0; i < no_dp_rows; ++i)
				dp_rows[i].resize(prof2_width + 1);
			});

		auto fut1 = async([&] {
			matrix.set_zeros(params->instruction_set);
			});

		auto fut2 = async([&] {
			gap_corrections.resize(prof2_width + 1);

			n_gaps_prof2_to_change.resize(prof2_width + 1);
			n_gaps_prof2_term_to_change.resize(prof2_width + 1);
			gaps_prof2_change.resize(prof2_width + 1);

			v_gap_corr.resize(prof2_width + 1);
			});

		prof2_gaps.resize(prof2_width + 1);

		for (size_t j = 0; j <= prof2_width; ++j)
		{
			prof2_gaps[j].open = scores2.get_value(j, GAP_OPEN);
			prof2_gaps[j].ext = scores2.get_value(j, GAP_EXT);
			prof2_gaps[j].term_open = scores2.get_value(j, GAP_TERM_OPEN);
			prof2_gaps[j].term_ext = scores2.get_value(j, GAP_TERM_EXT);
		}

		fut0.wait();
		fut1.wait();
		fut2.wait();
	}
	else if (no_threads == 3)
	{
		auto fut0 = async([&] {
			for (int i = 0; i < no_dp_rows; ++i)
				dp_rows[i].resize(prof2_width + 1);
			});

		auto fut1 = async([&] {
			matrix.set_zeros(params->instruction_set);
			});
		
		prof2_gaps.resize(prof2_width + 1);
		gap_corrections.resize(prof2_width + 1);

		n_gaps_prof2_to_change.resize(prof2_width + 1);
		n_gaps_prof2_term_to_change.resize(prof2_width + 1);
		gaps_prof2_change.resize(prof2_width + 1);

		v_gap_corr.resize(prof2_width + 1);

		for (size_t j = 0; j <= prof2_width; ++j)
		{
			prof2_gaps[j].open = scores2.get_value(j, GAP_OPEN);
			prof2_gaps[j].ext = scores2.get_value(j, GAP_EXT);
			prof2_gaps[j].term_open = scores2.get_value(j, GAP_TERM_OPEN);
			prof2_gaps[j].term_ext = scores2.get_value(j, GAP_TERM_EXT);
		}

		fut0.wait();
		fut1.wait();
	}
	else if (no_threads == 2)
	{
		auto fut = async([&] {
			for (int i = 0; i < no_dp_rows; ++i)
				dp_rows[i].resize(prof2_width + 1);
			});

		matrix.set_zeros(params->instruction_set);
		prof2_gaps.resize(prof2_width + 1);
		gap_corrections.resize(prof2_width + 1);

		n_gaps_prof2_to_change.resize(prof2_width + 1);
		n_gaps_prof2_term_to_change.resize(prof2_width + 1);
		gaps_prof2_change.resize(prof2_width + 1);

		v_gap_corr.resize(prof2_width + 1);

		for (size_t j = 0; j <= prof2_width; ++j)
		{
			prof2_gaps[j].open = scores2.get_value(j, GAP_OPEN);
			prof2_gaps[j].ext = scores2.get_value(j, GAP_EXT);
			prof2_gaps[j].term_open = scores2.get_value(j, GAP_TERM_OPEN);
			prof2_gaps[j].term_ext = scores2.get_value(j, GAP_TERM_EXT);
		}

		fut.wait();
	}
	else
	{
		matrix.set_zeros(params->instruction_set);

		for (int i = 0; i < no_dp_rows; ++i)
			dp_rows[i].resize(prof2_width + 1);
		
		prof2_gaps.resize(prof2_width + 1);
		gap_corrections.resize(prof2_width + 1);

		n_gaps_prof2_to_change.resize(prof2_width + 1);
		n_gaps_prof2_term_to_change.resize(prof2_width + 1);
		gaps_prof2_change.resize(prof2_width + 1);

		v_gap_corr.resize(prof2_width + 1);

		for (size_t j = 0; j <= prof2_width; ++j)
		{
			prof2_gaps[j].open = scores2.get_value(j, GAP_OPEN);
			prof2_gaps[j].ext = scores2.get_value(j, GAP_EXT);
			prof2_gaps[j].term_open = scores2.get_value(j, GAP_TERM_OPEN);
			prof2_gaps[j].term_ext = scores2.get_value(j, GAP_TERM_EXT);
		}
	}

	// Prepare ranges for threads
	vector<pair<int, int>> v_thr_range;
	
	if (prof2_width < no_threads)
		no_threads = (uint32_t) prof2_width;

	for (int i = 0; i < (int) no_threads; ++i)
		v_thr_range.emplace_back(1 + i * prof2_width / no_threads, 1 + (i + 1) * prof2_width / no_threads);

	// Boundary conditions
	dp_rows[0][0].D = 0;
	dp_rows[0][0].H = -infty;
	dp_rows[0][0].V = -infty;

	if (prof2_width >= 1)
	{
		dp_rows[0][1].D = -infty;

		dp_rows[0][1].H = dp_rows[0][0].D + prof2_gaps[1].term_open;

		dp_rows[0][1].V = -infty;
		matrix.set_dir_all(0, 1, direction_t::H);
	}

	for (size_t j = 2; j <= prof2_width; ++j)
	{
		//prev_row[j].D = -infty;
		dp_rows[0][j].D = -infty;

		dp_rows[0][j].H = dp_rows[0][j - 1].H + prof2_gaps[j].term_ext;

		dp_rows[0][j].V = -infty;
		matrix.set_dir_all(0, j, direction_t::H);
	}

	dp_rows[0][prof2_width].H = -infty;

	score_t prof1_gap_term_open_curr_card = gap_term_open * prof2_card;
	score_t prof1_gap_term_ext_curr_card = gap_term_ext * prof2_card;

	for (size_t j = 1; j <= prof2_width; ++j)
	{
		DP_SolveGapsProblemWhenStarting(j, prof2_width, prof2_card, profile2,
			gap_corrections[j].n_gap_start_open, gap_corrections[j].n_gap_start_ext, gap_corrections[j].n_gap_start_term_open, gap_corrections[j].n_gap_start_term_ext);

		DP_SolveGapsProblemWhenContinuing(j, prof2_width, prof2_card, profile2,
			gap_corrections[j].n_gap_cont_ext, gap_corrections[j].n_gap_cont_term_ext);

#ifndef NO_GAP_CORRECTIONS
		n_gaps_prof2_to_change[j] = profile2->counters.get_value(j, GAP_OPEN);      //the number of gaps to be changed from gap_open into gap_ext 
		n_gaps_prof2_term_to_change[j] = profile2->counters.get_value(j, GAP_TERM_OPEN); //the number of gaps to be changed from term_open into term_ext
#else
		n_gaps_prof2_to_change[j] = 0;
		n_gaps_prof2_term_to_change[j] = 0;
#endif

		gaps_prof2_change[j] = n_gaps_prof2_to_change[j] * (gap_ext - gap_open) +
			n_gaps_prof2_term_to_change[j] * (gap_term_ext - gap_term_open);

		v_gap_corr[j].first = gap_open * gap_corrections[j].n_gap_start_open +
			gap_ext * gap_corrections[j].n_gap_start_ext +
			gap_term_open * gap_corrections[j].n_gap_start_term_open +
			gap_term_ext * gap_corrections[j].n_gap_start_term_ext;

		v_gap_corr[j].second = 
			gap_ext * gap_corrections[j].n_gap_cont_ext +
			gap_term_ext * gap_corrections[j].n_gap_cont_term_ext;
	}

#ifdef USE_ASYNC
	vector<future<void>> v_fut(no_threads);
#else
//	vector<thread> v_thr;
//	v_thr.reserve(no_threads);
	vector<pooled_threads::thread> v_thr;
#endif
	barrier bar(no_threads);

	for (int t = 0; t < (int) no_threads; ++t)
#ifdef USE_ASYNC
		v_fut[t] = async([&, t]{
#else
		v_thr.emplace_back([&, t] {
#endif
			int my_id = t;
			int curr_row_id = 1;
			int my_col_from = v_thr_range[my_id].first;
			int my_col_to = v_thr_range[my_id].second;

			for (int i = 0; i < my_id; ++i)
				bar.arrive_and_wait();

			// Calculate matrix interior
			for (size_t i = 1; i <= prof1_width; ++i)
			{
				auto& curr_row = dp_rows[curr_row_id];
				auto& prev_row = dp_rows[(curr_row_id + no_dp_rows - 1) % no_dp_rows];

				// Boundary conditions
				if (my_id == 0)
				{
					curr_row[0].D = -infty;
					curr_row[0].H = -infty;
					matrix.set_dir_all(i, 0, direction_t::V);

					if (i < prof1_width)
					{
						if (i == 1)
							curr_row[0].V = max(prev_row[0].D, prev_row[0].V) + prof1_gap_term_open_curr_card;
						else
							curr_row[0].V = max(prev_row[0].D, prev_row[0].V) + prof1_gap_term_ext_curr_card;
					}
					else
						curr_row[0].V = -infty;
				}

				unsigned char* matrix_cell = matrix.get_cell(i, my_col_from-1);

				auto ptr_prev_row = &prev_row[my_col_from - 1];
				auto ptr_curr_row = &curr_row[my_col_from - 1];

				for (int j = my_col_from; j < my_col_to; ++j)
				{
					// Get current cell of matrix for faster access
					++matrix_cell;
					//...........................................................................
					// Calcualte score for D
					//...........................................................................

					//analyzing an alignment of the column j of the profile 2 with the column i of the profile 1
					score_t* scores2_column = scores2.get_column(j);

					//an alignment of a column with another column, after an an alignment of two columns  
					score_t t = scores2_column[seq1[i]];
					score_t t_D = ptr_prev_row->D;

					//..........................................................................
					//an alignment of a column with another column, after inserting a column of gaps before, into the profile 1
					score_t t_H = ptr_prev_row->H;

					//..........................................................................
					//an alignment of a column with another column, after inserting a column of gaps before, into the profile 2 
					score_t t_V = ptr_prev_row->V;

					//One has to check if there is a need of exchanging gap_opens/gap_term_opens (from the column of the profile 2 that is being just aligned) 
					//into gap_extensions/gap_term_extensions
					//If so, the value of the variable t has to be corrected

					t_V += gaps_prof2_change[j];

					if (t_D > t_H && t_D > t_V)
					{
						(ptr_curr_row+1)->D = t_D + t;
						matrix.set_dir_D(matrix_cell, direction_t::D);
					}
					else if (t_H > t_V)
					{
						(ptr_curr_row+1)->D = t_H + t;
						matrix.set_dir_D(matrix_cell, direction_t::H);
					}
					else
					{
						(ptr_curr_row+1)->D = t_V + t;
						matrix.set_dir_D(matrix_cell, direction_t::V);
					}

					//...........................................................................
					// Calcualte score for H
					//...........................................................................
					//inserting a columns of gaps into the profile 1

					//the case when a column (i) was aligned with the column (j-1)

					score_t gap_corr = (i < prof1_width ? prof2_gaps[j].open : prof2_gaps[j].term_open);
					t_D = ptr_curr_row->D + gap_corr;

					//the case when a sequence of gaps is continued (another columnn of gaps is inserted)
					t_H = ptr_curr_row->H + (i < prof1_width ? prof2_gaps[j].ext : prof2_gaps[j].term_ext);

#ifdef ALWAYS_3_DIRS
					//the column of gaps is inserted into the profile 1 after a column of gaps in the profile 2
					if ((i > 1) && (j > 1))
					{
						t_V = ptr_curr_row->V + gap_corr;

						if (t_D > t_H && t_D > t_V)
						{
							(ptr_curr_row+1)->H = t_D;
							matrix.set_dir_H(matrix_cell, direction_t::D);
						}
						else if (t_V > t_H)	//Attention!!! Swapping the checking order
						{
							(ptr_curr_row + 1)->H = t_V;
							matrix.set_dir_H(matrix_cell, direction_t::V);
						}
						else
						{
							(ptr_curr_row + 1)->H = t_H;
							matrix.set_dir_H(matrix_cell, direction_t::H);
						}
					}
					else
					{
						if (t_D > t_H)
						{
							(ptr_curr_row + 1)->H = t_D;
							matrix.set_dir_H(matrix_cell, direction_t::D);
						}
						else
						{
							(ptr_curr_row + 1)->H = t_H;
							matrix.set_dir_H(matrix_cell, direction_t::H);
						}
					}
#else // !ALWAYS_3_DIRS
					if (t_D > t_H)
					{
						(ptr_curr_row + 1)->H = t_D;
						matrix.set_dir_H(matrix_cell, direction_t::D);
					}
					else
					{
						(ptr_curr_row + 1)->H = t_H;
						matrix.set_dir_H(matrix_cell, direction_t::H);
					}
#endif // !ALWAYS_3_DIRS
					//...........................................................................
					// Calcualte score for V
					//...........................................................................
					//inserting a columns of gaps into the profile 2

					++ptr_prev_row;
					++ptr_curr_row;

					//the case when a column (i-1) was aligned with the column (j)
					gap_corr = v_gap_corr[j].first;
					t_D = ptr_prev_row->D + gap_corr;

					//the case when a sequence of gaps is continued (another columnn of gaps is inserted)
					t_V = ptr_prev_row->V +	v_gap_corr[j].second;

#ifdef ALWAYS_3_DIRS
					if ((i > 1) && (j > 1))
					{
						t_H = ptr_prev_row->H + gap_corr;

						if (t_D > t_H && t_D > t_V)
						{
							ptr_curr_row->V = t_D;
							matrix.set_dir_V(matrix_cell, direction_t::D);
						}
						else if (t_H > t_V)
						{
							ptr_curr_row->V = t_H;
							matrix.set_dir_V(matrix_cell, direction_t::H);
						}
						else
						{
							ptr_curr_row->V = t_V;
							matrix.set_dir_V(matrix_cell, direction_t::V);
						}
					}
					else
					{
						if (t_D > t_V)
						{
							ptr_curr_row->V = t_D;
							matrix.set_dir_V(matrix_cell, direction_t::D);
						}
						else
						{
							ptr_curr_row->V = t_V;
							matrix.set_dir_V(matrix_cell, direction_t::V);
						}
					}
#else
					if (t_D > t_V)
					{
						ptr_curr_row->V = t_D;
						matrix.set_dir_V(matrix_cell, direction_t::D);
					}
					else
					{
						ptr_curr_row->V = t_V;
						matrix.set_dir_V(matrix_cell, direction_t::V);
					}
#endif //ALWAYS_3_DIRS
				}

				curr_row_id = (curr_row_id + 1) % no_dp_rows;

				if(i % rows_per_box == 0 || i == prof1_width)
					bar.arrive_and_wait();
			}

			for (size_t i = my_id; i < no_threads - 1; ++i)
				bar.arrive_and_wait();
		
			});

#ifdef USE_ASYNC
		for (auto& v : v_fut)
			v.wait();
#else
	for (auto& t : v_thr)
		t.join();
#endif

	// Construct alignment
	ConstructProfile(profile1, profile2, matrix, dp_rows[prof1_width % no_dp_rows].back(), no_threads);
}

// ****************************************************************************
void CProfile::ParAlignProfProf(CProfile* profile1, CProfile* profile2, uint32_t no_threads, uint32_t rows_per_box)
{
	uint32_t no_dp_rows = 3 * rows_per_box;

	size_t prof1_width = profile1->width;
	size_t prof2_width = profile2->width;

	size_t prof1_card = profile1->data.size();
	size_t prof2_card = profile2->data.size();

	score_t gap_open = params->gap_open;
	score_t gap_ext = params->gap_ext;
	score_t gap_term_open = params->gap_term_open;
	score_t gap_term_ext = params->gap_term_ext;

	CDPMatrix matrix(prof1_width + 1, prof2_width + 1);
	vector<dp_row_t> dp_rows(no_dp_rows);

	CProfileValues<score_t, NO_SYMBOLS>& scores1 = const_cast<CProfileValues<score_t, NO_SYMBOLS>&>(profile1->scores);
	CProfileValues<score_t, NO_SYMBOLS>& scores2 = const_cast<CProfileValues<score_t, NO_SYMBOLS>&>(profile2->scores);

	// Precompute scores for gaps for profile2
	vector<dp_gap_costs> prof2_gaps;
	vector<dp_gap_corrections> gap_corrections2;

	if (no_threads > 2)
	{
		auto fut0 = async([&] {
			for (int i = 0; i < (int) no_dp_rows; ++i)
				dp_rows[i].resize(prof2_width + 1);
			});

		auto fut1 = async([&] {
			matrix.set_zeros(params->instruction_set);
			});

		prof2_gaps.resize(prof2_width + 1);
		gap_corrections2.resize(prof2_width + 1);

		fut0.wait();
		fut1.wait();
	}
	else if (no_threads == 2)
	{
		auto fut = async([&] {
			for (int i = 0; i < (int) no_dp_rows; ++i)
				dp_rows[i].resize(prof2_width + 1);
			});

		matrix.set_zeros(params->instruction_set);
		prof2_gaps.resize(prof2_width + 1);
		gap_corrections2.resize(prof2_width + 1);

		fut.wait();
	}
	else
	{
		matrix.set_zeros(params->instruction_set);

		for (int i = 0; i < (int) no_dp_rows; ++i)
			dp_rows[i].resize(prof2_width + 1);

		prof2_gaps.resize(prof2_width + 1);
		gap_corrections2.resize(prof2_width + 1);
	}

	// Prepare ranges for threads
	vector<pair<int, int>> v_thr_range;

	if (prof2_width < no_threads)
		no_threads = (uint32_t) prof2_width;

	for (int i = 0; i < (int) no_threads; ++i)
		v_thr_range.emplace_back(1 + i * prof2_width / no_threads, 1 + (i + 1) * prof2_width / no_threads);

	//scoresX.get_value(j, GAP_OPEN) returns information of how the column j of the profile X aligns with a single gap_open
	//scoresX.get_value(j, GAP_EXT) returns information of how the column j of the profile X aligns with a single gap_ext
	//etc.
	//The task is to calculate how the column j of the profile X aligns with a column of gaps of different categories

	for (size_t j = 0; j <= prof2_width; ++j)
	{
		prof2_gaps[j].open = scores2.get_value(j, GAP_OPEN);
		prof2_gaps[j].ext = scores2.get_value(j, GAP_EXT);
		prof2_gaps[j].term_open = scores2.get_value(j, GAP_TERM_OPEN);
		prof2_gaps[j].term_ext = scores2.get_value(j, GAP_TERM_EXT);
	}

	// Boundary conditions
	dp_rows[0][0].D = 0;
	dp_rows[0][0].H = -infty;
	dp_rows[0][0].V = -infty;

	if (prof2_width >= 1)
	{
		dp_rows[0][1].D = -infty;

		dp_rows[0][1].H = dp_rows[0][0].D + prof2_gaps[1].term_open * prof1_card;

		dp_rows[0][1].V = -infty;
		matrix.set_dir_all(0, 1, direction_t::H);
	}

	for (size_t j = 2; j <= prof2_width; ++j)
	{
		dp_rows[0][j].D = -infty;

		dp_rows[0][j].H = dp_rows[0][j - 1].H + prof2_gaps[j].term_ext * prof1_card;

		dp_rows[0][j].V = -infty;
		matrix.set_dir_all(0, j, direction_t::H);
	}

	dp_rows[0][prof2_width].H = -infty;

	// Precomputing for gap correction in the profile2
	vector<score_t> n_gaps_prof2_to_change(prof2_width + 1);
	vector<score_t> n_gaps_prof2_term_to_change(prof2_width + 1);
	vector<score_t> gaps_prof2_change(prof2_width + 1);

	for (size_t j = 1; j <= prof2_width; ++j)
	{
		DP_SolveGapsProblemWhenStarting(j, prof2_width, prof2_card, profile2,
			gap_corrections2[j].n_gap_start_open, gap_corrections2[j].n_gap_start_ext, gap_corrections2[j].n_gap_start_term_open, gap_corrections2[j].n_gap_start_term_ext);

		DP_SolveGapsProblemWhenContinuing(j, prof2_width, prof2_card, profile2,
			gap_corrections2[j].n_gap_cont_ext, gap_corrections2[j].n_gap_cont_term_ext);

#ifndef NO_GAP_CORRECTION
		n_gaps_prof2_to_change[j] = profile2->counters.get_value(j, GAP_OPEN);      //the number of gaps to be changed from gap_open into gap_ext 
		n_gaps_prof2_term_to_change[j] = profile2->counters.get_value(j, GAP_TERM_OPEN); //the number of gaps to be changed from term_open into term_ext

#else
		n_gaps_prof2_to_change[j] = 0;
		n_gaps_prof2_term_to_change[j] = 0;
#endif

		gaps_prof2_change[j] = n_gaps_prof2_to_change[j] * (gap_ext - gap_open) +
			n_gaps_prof2_term_to_change[j] * (gap_term_ext - gap_term_open);
	}

	vector<dp_gap_corrections> gap_corrections1(prof1_width + 1);
	// Precomputing for gap correction in the profile1
	for (size_t i = 1; i <= prof1_width; ++i)
	{
		DP_SolveGapsProblemWhenStarting(i, prof1_width, prof1_card, profile1,
			gap_corrections1[i].n_gap_start_open, gap_corrections1[i].n_gap_start_ext, gap_corrections1[i].n_gap_start_term_open, gap_corrections1[i].n_gap_start_term_ext);

		DP_SolveGapsProblemWhenContinuing(i, prof1_width, prof1_card, profile1,
			gap_corrections1[i].n_gap_cont_ext, gap_corrections1[i].n_gap_cont_term_ext);
	}

#ifdef USE_ASYNC
	vector<future<void>> v_fut(no_threads);
#else
	//	vector<thread> v_thr;
	//	v_thr.reserve(no_threads);
	vector<pooled_threads::thread> v_thr;
#endif
	barrier bar(no_threads);

	for (int t = 0; t < (int) no_threads; ++t)
#ifdef USE_ASYNC
		v_fut[t] = async([&, t] {
#else
		v_thr.emplace_back([&, t] {
#endif
			int my_id = t;
			int curr_row_id = 1;
			int my_col_from = v_thr_range[my_id].first;
			int my_col_to = v_thr_range[my_id].second;

			// Frequency of symbols at all columns of profile1 and profile2
//			vector<pair<int, int>> col1(32);
			array<pair<int, int>, 32> col1;

			for (int i = 0; i < my_id; ++i)
				bar.arrive_and_wait();

			// Calculate matrix interior
			for (size_t i = 1; i <= prof1_width; ++i)
			{
				auto& curr_row = dp_rows[curr_row_id];
				auto& prev_row = dp_rows[(curr_row_id + no_dp_rows - 1) % no_dp_rows];

				// Precompute scores for gaps for current and previous row of profile1
				//score_t prof1_gap_open_prev = scores1.get_value(i - 1, GAP_OPEN);
				score_t prof1_gap_open_curr = scores1.get_value(i, GAP_OPEN);
				//score_t prof1_gap_term_open_prev = scores1.get_value(i - 1, GAP_TERM_OPEN);
				score_t prof1_gap_term_open_curr = scores1.get_value(i, GAP_TERM_OPEN);
				score_t prof1_gap_ext_curr = scores1.get_value(i, GAP_EXT);
				score_t prof1_gap_term_ext_curr = scores1.get_value(i, GAP_TERM_EXT);

				// Boundary conditions
				if (my_id == 0)
				{
					curr_row[0].D = -infty;
					curr_row[0].H = -infty;
					matrix.set_dir_all(i, 0, direction_t::V);

					if (i < prof1_width)
					{
						if (i == 1)
							curr_row[0].V = max(prev_row[0].D, prev_row[0].V) + prof1_gap_term_open_curr * prof2_card;
						else
							curr_row[0].V = max(prev_row[0].D, prev_row[0].V) + prof1_gap_term_ext_curr * prof2_card;
					}
					else
						curr_row[0].V = -infty;
				}

				// Calculate frequency of symbols in (i-1)-th column of profile1
/*				col1.clear();
				size_t col1_n_non_gaps = 0;
				for (size_t k = 0; k < NO_AMINOACIDS_AND_GAPS; ++k)
					if (profile1->counters.get_value(i, k))
					{
						size_t count = profile1->counters.get_value(i, k);
						col1.emplace_back(k, count);
						if (k < NO_AMINOACIDS)
							col1_n_non_gaps += count;
					}
				size_t col1_size = col1.size();*/

				size_t col1_size = 0;
				size_t col1_n_non_gaps = 0;
				for (size_t k = 0; k < NO_AMINOACIDS_AND_GAPS; ++k)
					if (profile1->counters.get_value(i, k))
					{
						size_t count = profile1->counters.get_value(i, k);
						col1[col1_size++] = make_pair(k, count);
						if (k < NO_AMINOACIDS)
							col1_n_non_gaps += count;
					}

#ifndef NO_GAP_CORRECTION
				size_t n_gaps_prof1_to_change = profile1->counters.get_value(i, GAP_OPEN);   //the number of gaps to be changed from gap_open into gap_ext
				size_t n_gaps_prof1_term_to_change = profile1->counters.get_value(i, GAP_TERM_OPEN); //the number of gaps to be changed from term_open into term_ext
#else
				size_t n_gaps_prof1_to_change = 0;
				size_t n_gaps_prof1_term_to_change = 0;
#endif

				unsigned char* matrix_cell = matrix.get_cell(i, my_col_from - 1);

				auto ptr_prev_row = &prev_row[my_col_from - 1];
				auto ptr_curr_row = &curr_row[my_col_from - 1];

				for (int j = my_col_from; j < my_col_to; ++j)
				{
					// Get current cell of matrix for faster access
					++matrix_cell;

					//...........................................................................
					// Calcualte score for D
					//...........................................................................

					//analyzing an alignment of the column j of the profile 2 with the column i of the profile 1
					score_t* scores2_column = scores2.get_column(j);

					//scores2_column contains information on how the column j aligns with a single character that could appear in the profile 1
					//the variable t is to store a cost of the alignment of the column j with a specific character (col1[k].first) 
					//of a concrete number of occurrences (col1[k].second) of the profile 1

					score_t t = 0;
					switch (col1_size)
					{
					case 3:
						t += col1[2].second * scores2_column[col1[2].first];
						FALL_THROUGH;
					case 2:
						t += col1[1].second * scores2_column[col1[1].first];
						FALL_THROUGH;
					case 1:
						t += col1[0].second * scores2_column[col1[0].first];
					case 0:
						break;
					default:
						for (size_t k = 0; k < col1_size; ++k)
							t += col1[k].second * scores2_column[col1[k].first];
					}

					//an alignment of a column with another column, after an an alignment of two columns  
					score_t t_D = ptr_prev_row->D + t;

					//..........................................................................
					//an alignment of a column with another column, after inserting a column of gaps before, into the profile 1
					score_t t_H = ptr_prev_row->H;

					//One has to check if there is a need of exchanging gap_opens/gap_term_opens (from the column of the profile 1 that is being just aligned) 
					//into gap_extensions/gap_term_extensions
					//If so, the value of the variable t has to be corrected, because it was computed with:
					//t += col1[k].liczba_gap_open * scores2_column[GAP_OPEN];
					//t += col1[k].liczba_gap_ext * scores2_column[GAP_EXT];
					//t += col1[k].liczba_gap_term_open * scores2_column[GAP_TERM_OPEN];
					//t += col1[k].liczba_gap_term_ext * scores2_column[GAP_TERM_EXT];

					if (n_gaps_prof1_to_change || n_gaps_prof1_term_to_change)
					{
						//correction by delta_t
						score_t delta_t = n_gaps_prof1_to_change * (scores2_column[GAP_EXT] - scores2_column[GAP_OPEN]) +
							n_gaps_prof1_term_to_change * (scores2_column[GAP_TERM_EXT] - scores2_column[GAP_TERM_OPEN]);

						t_H = t_H + t + delta_t;
					}
					else
						t_H = t_H + t;

					//..........................................................................
					//an alignment of a column with another column, after inserting a column of gaps before, into the profile 2
					score_t t_V = ptr_prev_row->V;

					t_V += t + gaps_prof2_change[j] * col1_n_non_gaps;

					if (t_D > t_H && t_D > t_V)
					{
						(ptr_curr_row+1)->D = t_D;
						matrix.set_dir_D(matrix_cell, direction_t::D);
					}
					else if (t_H > t_V)
					{
						(ptr_curr_row + 1)->D = t_H;
						matrix.set_dir_D(matrix_cell, direction_t::H);
					}
					else
					{
						(ptr_curr_row + 1)->D = t_V;
						matrix.set_dir_D(matrix_cell, direction_t::V);
					}

					//...........................................................................
					// Calcualte score for H
					//...........................................................................
					//inserting a columns of gaps into the profile 1

					//the case when a column (i) was aligned with the column (j-1)
					score_t gap_corr = prof2_gaps[j].open * gap_corrections1[i].n_gap_start_open + prof2_gaps[j].ext * gap_corrections1[i].n_gap_start_ext +
						prof2_gaps[j].term_open * gap_corrections1[i].n_gap_start_term_open + prof2_gaps[j].term_ext * gap_corrections1[i].n_gap_start_term_ext;

					t_D = ptr_curr_row->D + gap_corr;

					//the case when a sequence of gaps is continued (another columnn of gaps is inserted)
					t_H = ptr_curr_row->H + prof2_gaps[j].ext * gap_corrections1[i].n_gap_cont_ext +
						prof2_gaps[j].term_ext * gap_corrections1[i].n_gap_cont_term_ext;

#ifdef ALWAYS_3_DIRS
					//the column of gaps is inserted into the profile 1 after a column of gaps in the profile 2
					if ((i > 1) && (j > 1))
					{
						t_V = ptr_curr_row->V + gap_corr;

						if (t_D > t_H && t_D > t_V)
						{
							(ptr_curr_row+1)->H = t_D;
							matrix.set_dir_H(matrix_cell, direction_t::D);
						}
						else if (t_V > t_H)	//Attention!!! Swapping the checking order 
						{
							(ptr_curr_row + 1)->H = t_V;
							matrix.set_dir_H(matrix_cell, direction_t::V);
						}
						else
						{
							(ptr_curr_row + 1)->H = t_H;
							matrix.set_dir_H(matrix_cell, direction_t::H);
						}
					}
					else
					{
						if (t_D > t_H)
						{
							(ptr_curr_row + 1)->H = t_D;
							matrix.set_dir_H(matrix_cell, direction_t::D);
						}
						else
						{
							(ptr_curr_row + 1)->H = t_H;
							matrix.set_dir_H(matrix_cell, direction_t::H);
						}
					}
#else // !ALWAYS_3_DIRS
					if (t_D > t_H)
					{
						(ptr_curr_row + 1)->H = t_D;
						matrix.set_dir_H(matrix_cell, direction_t::D);
					}
					else
					{
						(ptr_curr_row + 1)->H = t_H;
						matrix.set_dir_H(matrix_cell, direction_t::H);
					}
#endif // !ALWAYS_3_DIRS
					//...........................................................................
					// Calcualte score for V
					//...........................................................................
					//inserting a columns of gaps into the profile 2

					++ptr_prev_row;
					++ptr_curr_row;

					//the case when a column (i-1) was aligned with the column (j)
					gap_corr = prof1_gap_open_curr * gap_corrections2[j].n_gap_start_open + prof1_gap_ext_curr * gap_corrections2[j].n_gap_start_ext +
						prof1_gap_term_open_curr * gap_corrections2[j].n_gap_start_term_open + prof1_gap_term_ext_curr * gap_corrections2[j].n_gap_start_term_ext;
					t_D = ptr_prev_row->D + gap_corr;

					//the case when a sequence of gaps is continued (another columnn of gaps is inserted)
					t_V = ptr_prev_row->V + prof1_gap_ext_curr * gap_corrections2[j].n_gap_cont_ext +
						prof1_gap_term_ext_curr * gap_corrections2[j].n_gap_cont_term_ext;

#ifdef ALWAYS_3_DIRS
					if ((i > 1) && (j > 1))
					{
						t_H = ptr_prev_row->H + gap_corr;

						if (t_D > t_H && t_D > t_V)
						{
							ptr_curr_row->V = t_D;
							matrix.set_dir_V(matrix_cell, direction_t::D);
						}
						else if (t_H > t_V)
						{
							ptr_curr_row->V = t_H;
							matrix.set_dir_V(matrix_cell, direction_t::H);
						}
						else
						{
							ptr_curr_row->V = t_V;
							matrix.set_dir_V(matrix_cell, direction_t::V);
						}
					}
					else
					{
						if (t_D > t_V)
						{
							ptr_curr_row->V = t_D;
							matrix.set_dir_V(matrix_cell, direction_t::D);
						}
						else
						{
							ptr_curr_row->V = t_V;
							matrix.set_dir_V(matrix_cell, direction_t::V);
						}
					}
#else
					if (t_D > t_V)
					{
						ptr_curr_row->V = t_D;
						matrix.set_dir_V(matrix_cell, direction_t::D);
					}
					else
					{
						ptr_curr_row->V = t_V;
						matrix.set_dir_V(matrix_cell, direction_t::V);
					}
#endif //ALWAYS_3_DIRS
				}

				curr_row_id = (curr_row_id + 1) % no_dp_rows;

				if (i % rows_per_box == 0 || i == prof1_width)
					bar.arrive_and_wait();
			}

			for (size_t i = my_id; i < no_threads - 1; ++i)
				bar.arrive_and_wait();
			
			});

#ifdef USE_ASYNC
	for (auto& v : v_fut)
		v.wait();
#else
	for (auto& t : v_thr)
		t.join();
#endif

	// Construct alignment
	ConstructProfile(profile1, profile2, matrix, dp_rows[prof1_width % no_dp_rows].back(), no_threads);
}

#endif

// EOF
