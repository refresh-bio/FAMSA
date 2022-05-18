/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#include "../core/profile.h"
#include "../core/sequence.h"
#include <assert.h>
#include <array>
#include <set>
#include <math.h>
#include <algorithm>
#include <stdio.h>

#define max3(x, y, z)	(max((x), max((y), (z))))

#define NINF(x)		((x) < -999 ? -999: (x))

// ****************************************************************************
// 
void CProfile::AlignSeqSeq(CProfile* profile1, CProfile* profile2)
{
	size_t prof1_width = profile1->width;
	size_t prof2_width = profile2->width;

	// These profiles contain single sequences so there are no gaps here
//	vector<symbol_t>& seq1 = const_cast<vector<symbol_t>&>(profile1->data[0]->symbols);
//	vector<symbol_t>& seq2 = const_cast<vector<symbol_t>&>(profile2->data[0]->symbols);
	symbol_t* seq1 = const_cast<symbol_t*>(profile1->data[0]->symbols);
	symbol_t* seq2 = const_cast<symbol_t*>(profile2->data[0]->symbols);

	score_t gap_open = params->gap_open;
	score_t gap_ext = params->gap_ext;
	score_t gap_term_open = params->gap_term_open;
	score_t gap_term_ext = params->gap_term_ext;

	CDPMatrix matrix(prof1_width + 1, prof2_width + 1);
	matrix.set_zeros(params->instruction_set);

	dp_row_t curr_row(prof2_width + 1);
	dp_row_t prev_row(prof2_width + 1);

	// Boundary conditions
	prev_row[0].D = 0;
	prev_row[0].H = -infty;
	prev_row[0].V = -infty;

	if (1 <= prof2_width)
	{
		prev_row[1].D = -infty;
		prev_row[1].H = gap_term_open;
		prev_row[1].V = -infty;
		matrix.set_dir_all(0, 1, direction_t::H);
	}

	for (size_t j = 2; j <= prof2_width; ++j)
	{
		prev_row[j].D = -infty;
		prev_row[j].H = max(prev_row[j - 1].H, prev_row[j - 1].D) + gap_term_ext;
		prev_row[j].V = -infty;
		matrix.set_dir_all(0, j, direction_t::H);
	}
	prev_row[prof2_width].H = -infty;

	// Calculate matrix interior
	for (size_t i = 1; i <= prof1_width; ++i)
	{
		// Boundary conditions
		curr_row[0].D = -infty;
		curr_row[0].H = -infty;
		matrix.set_dir_all(i, 0, direction_t::V);

		if (i < prof1_width)
		{

			if (i == 1)
				curr_row[0].V = max(prev_row[0].V, prev_row[0].D) + gap_term_open;
			else
				curr_row[0].V = max(prev_row[0].V, prev_row[0].D) + gap_term_ext;
		}
		else
			curr_row[0].V = -infty;

		auto score_row = params->score_matrix[seq1[i]];

		for (size_t j = 1; j <= prof2_width; ++j)
		{
			//...........................................................................
			// Calcualte score for D
			//...........................................................................
			score_t t_D = prev_row[j - 1].D; // aligning two characters after an alignment of two characters placed before
			score_t t_H = prev_row[j - 1].H; // aligning two characters after an alignment of a character with a gap
			score_t t_V = prev_row[j - 1].V; // aligning two characters after an alignment of a character with a gap

			if (t_D > t_H && t_D > t_V)
			{
				curr_row[j].D = t_D + score_row[seq2[j]];
				matrix.set_dir_D(i, j, direction_t::D);
			}
			else if (t_H >= t_V)
			{
				curr_row[j].D = t_H + score_row[seq2[j]];
				matrix.set_dir_D(i, j, direction_t::H);
				//???
			}
			else
			{
				curr_row[j].D = t_V + score_row[seq2[j]];
				matrix.set_dir_D(i, j, direction_t::V);
			}

			//...........................................................................
			// Calcualte score for H
			//...........................................................................

			//a gap is to be inserted into a profile

			//This gap takes the form of gap_open if a match was before
			t_D = curr_row[j - 1].D + (i < prof1_width ? gap_open : gap_term_open);

			//This gap takes the form of gap_ext if a gap was inserted before
			t_H = curr_row[j - 1].H + (i < prof1_width ? gap_ext : gap_term_ext);

			if (t_D > t_H)
			{
				curr_row[j].H = t_D;
				matrix.set_dir_H(i, j, direction_t::D);
			}
			else
			{
				curr_row[j].H = t_H;
				matrix.set_dir_H(i, j, direction_t::H);
			}

			//...........................................................................
			// Calculate score for V
			//...........................................................................

			t_D = prev_row[j].D + (j < prof2_width ? gap_open : gap_term_open);
			t_V = prev_row[j].V + (j < prof2_width ? gap_ext : gap_term_ext);

			if (t_D > t_V)
			{
				curr_row[j].V = t_D;
				matrix.set_dir_V(i, j, direction_t::D);
			}
			else
			{
				curr_row[j].V = t_V;
				matrix.set_dir_V(i, j, direction_t::V);
			}

		}

		curr_row.swap(prev_row);
	}

	ConstructProfile(profile1, profile2, matrix, prev_row.back());
}

// ****************************************************************************
void CProfile::AlignSeqProf(CProfile* profile1, CProfile* profile2, vector<int>* column_mapping1, vector<int>* column_mapping2)
{
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
	matrix.set_zeros(params->instruction_set);

	dp_row_t curr_row(prof2_width + 1);
	dp_row_t prev_row(prof2_width + 1);

	//CProfileValues<score_t, NO_SYMBOLS>& scores1 = const_cast<CProfileValues<score_t, NO_SYMBOLS>&>(profile1->scores);
	CProfileValues<score_t, NO_SYMBOLS>& scores2 = const_cast<CProfileValues<score_t, NO_SYMBOLS>&>(profile2->scores);

	bool is_guided = column_mapping1 != nullptr && column_mapping2 != nullptr;
	vector<pair<int, int>> row_ranges;

	if (is_guided)
		FindRowRanges(column_mapping1, column_mapping2, row_ranges);
	else
		row_ranges.assign(prof1_width + 1, make_pair(0, prof2_width));

	// Precompute scores for gaps for profile2
	vector<dp_gap_costs> prof2_gaps(prof2_width + 1);
	vector<dp_gap_corrections> gap_corrections(prof2_width + 1);

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
	prev_row[0].D = 0;
	prev_row[0].H = -infty;
	prev_row[0].V = -infty;

	if (prof2_width >= 1)
	{
		prev_row[1].D = -infty;

		//prev_row[1].H = prev_row[0].D + prof2_gaps[1].term_open * prof1_card;
		prev_row[1].H = prev_row[0].D + prof2_gaps[1].term_open;

		prev_row[1].V = -infty;
		matrix.set_dir_all(0, 1, direction_t::H);
	}


	for (size_t j = 2; j <= prof2_width; ++j)
	{
		prev_row[j].D = -infty;

		prev_row[j].H = prev_row[j - 1].H + prof2_gaps[j].term_ext;

		prev_row[j].V = -infty;
		matrix.set_dir_all(0, j, direction_t::H);
	}
	prev_row[prof2_width].H = -infty;

	// Precomputing korekcji gapów w profile2
	vector<score_t> n_gaps_prof2_to_change(prof2_width + 1);
	vector<score_t> n_gaps_prof2_term_to_change(prof2_width + 1);
	vector<score_t> gaps_prof2_change(prof2_width + 1);

	vector<pair<score_t, score_t>> v_gap_corr(prof2_width + 1);

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

	// Calculate matrix interior
	for (size_t i = 1; i <= prof1_width; ++i)
	{
		// Boundary conditions
		curr_row[0].D = -infty;
		curr_row[0].H = -infty;
		matrix.set_dir_all(i, 0, direction_t::V);

		if (row_ranges[i].first)
			curr_row[row_ranges[i].first - 1].D = curr_row[row_ranges[i].first - 1].H = curr_row[row_ranges[i].first - 1].V = -infty;

		if (i < prof1_width)
		{
			if (i == 1)
				curr_row[0].V = max(prev_row[0].D, prev_row[0].V) + prof1_gap_term_open_curr_card;
			else
				curr_row[0].V = max(prev_row[0].D, prev_row[0].V) + prof1_gap_term_ext_curr_card;

			for (int j = row_ranges[i].second + 1; j <= min(row_ranges[i + 1].second, (int)prof2_width); ++j)
				curr_row[j].D = curr_row[j].H = curr_row[j].V = -infty;
		}
		else
			curr_row[0].V = -infty;

		unsigned char* matrix_cell = matrix.get_cell(i, max(1, row_ranges[i].first) - 1);
		size_t max_j = min(row_ranges[i].second, (int)prof2_width);
		size_t min_j = max(1, row_ranges[i].first);

		auto ptr_prev_row = &prev_row[min_j - 1];
		auto ptr_curr_row = &curr_row[min_j - 1];

		auto seq1_i = seq1[i];

		for (size_t j = min_j; j <= max_j; ++j)
		{
			// Get current cell of matrix for faster access
			++matrix_cell;
			//...........................................................................
			// Calcualte score for D
			//...........................................................................

			//analyzing an alignment of the column j of the profile 2 with the column i of the profile 1
			score_t* scores2_column = scores2.get_column(j);

			//an alignment of a column with another column, after an an alignment of two columns  
			score_t t = scores2_column[seq1_i];
			score_t t_D = ptr_prev_row->D;

			//..........................................................................
			//an alignment of a column with another column, after inserting a column of gaps before, into the profile 1
			score_t t_H = ptr_prev_row->H;
//			t_H = t_H + t;

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
				(ptr_curr_row + 1)->D = t_H + t;
				matrix.set_dir_D(matrix_cell, direction_t::H);
			}
			else
			{
				(ptr_curr_row + 1)->D = t_V + t;
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
			t_V = ptr_prev_row->V + v_gap_corr[j].second;

#ifdef ALWAYS_3_DIRS
			if ((i > 1) && (j > 1))
			{
				t_H = prev_row[j].H + gap_corr;

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

		curr_row.swap(prev_row);
	}

	// Construct alignment
	ConstructProfile(profile1, profile2, matrix, prev_row.back());
}

// ****************************************************************************
// 
void CProfile::AlignProfProf(CProfile* profile1, CProfile* profile2, vector<int>* column_mapping1, vector<int>* column_mapping2)
{
	size_t prof1_width = profile1->width;
	size_t prof2_width = profile2->width;

	size_t prof1_card = profile1->data.size();
	size_t prof2_card = profile2->data.size();

	score_t gap_open = params->gap_open;
	score_t gap_ext = params->gap_ext;
	score_t gap_term_open = params->gap_term_open;
	score_t gap_term_ext = params->gap_term_ext;

	CDPMatrix matrix(prof1_width + 1, prof2_width + 1);
	matrix.set_zeros(params->instruction_set);

	dp_row_t curr_row(prof2_width + 1);
	dp_row_t prev_row(prof2_width + 1);

	// Frequency of symbols at all columns of profile1 and profile2
//	vector<pair<int, int>> col1(32);
	array<pair<int, int>, 32> col1;

	CProfileValues<score_t, NO_SYMBOLS>& scores1 = const_cast<CProfileValues<score_t, NO_SYMBOLS>&>(profile1->scores);
	CProfileValues<score_t, NO_SYMBOLS>& scores2 = const_cast<CProfileValues<score_t, NO_SYMBOLS>&>(profile2->scores);

	// Prepare ranges for computed cells in rows
	bool is_guided = column_mapping1 != nullptr && column_mapping2 != nullptr;
	vector<pair<int, int>> row_ranges;

	if (is_guided)
		FindRowRanges(column_mapping1, column_mapping2, row_ranges);
	else
		row_ranges.assign(prof1_width + 1, make_pair(0, prof2_width));

	// Precompute scores for gaps for profile2
	vector<dp_gap_costs> prof2_gaps(prof2_width + 1);

	//size_t n_gap_open = 0;
	//size_t n_gap_ext = 0;
	//size_t n_gap_term_open = 0;
	//size_t n_gap_term_ext = 0;
	
	counter_t n_gap_prof1_start_open, n_gap_prof1_start_ext, n_gap_prof1_start_term_open, n_gap_prof1_start_term_ext,
		n_gap_prof1_cont_ext, n_gap_prof1_cont_term_ext;
	vector<dp_gap_corrections> gap_corrections(prof2_width + 1);

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
	prev_row[0].D = 0;
	prev_row[0].H = -infty;
	prev_row[0].V = -infty;

	if (prof2_width >= 1)
	{
		prev_row[1].D = -infty;

		prev_row[1].H = prev_row[0].D + prof2_gaps[1].term_open * prof1_card;

		prev_row[1].V = -infty;
		matrix.set_dir_all(0, 1, direction_t::H);
	}

	for (size_t j = 2; j <= prof2_width; ++j)
	{
		prev_row[j].D = -infty;

		prev_row[j].H = prev_row[j - 1].H + prof2_gaps[j].term_ext * prof1_card;

		prev_row[j].V = -infty;
		matrix.set_dir_all(0, j, direction_t::H);
	}
	prev_row[prof2_width].H = -infty;

	// Precomputing for gap correction in the profile2
	vector<score_t> n_gaps_prof2_to_change(prof2_width + 1);
	vector<score_t> n_gaps_prof2_term_to_change(prof2_width + 1);
	vector<score_t> gaps_prof2_change(prof2_width + 1);

	for (size_t j = 1; j <= prof2_width; ++j)
	{
		DP_SolveGapsProblemWhenStarting(j, prof2_width, prof2_card, profile2,
			gap_corrections[j].n_gap_start_open, gap_corrections[j].n_gap_start_ext, gap_corrections[j].n_gap_start_term_open, gap_corrections[j].n_gap_start_term_ext);

		DP_SolveGapsProblemWhenContinuing(j, prof2_width, prof2_card, profile2,
			gap_corrections[j].n_gap_cont_ext, gap_corrections[j].n_gap_cont_term_ext);

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

	// Calculate matrix interior
	for (size_t i = 1; i <= prof1_width; ++i)
	{
		// Precompute scores for gaps for current and previous row of profile1
		//score_t prof1_gap_open_prev = scores1.get_value(i - 1, GAP_OPEN);
		score_t prof1_gap_open_curr = scores1.get_value(i, GAP_OPEN);
		//score_t prof1_gap_term_open_prev = scores1.get_value(i - 1, GAP_TERM_OPEN);
		score_t prof1_gap_term_open_curr = scores1.get_value(i, GAP_TERM_OPEN);
		score_t prof1_gap_ext_curr = scores1.get_value(i, GAP_EXT);
		score_t prof1_gap_term_ext_curr = scores1.get_value(i, GAP_TERM_EXT);

		// Boundary conditions
		curr_row[0].D = -infty;
		curr_row[0].H = -infty;
		matrix.set_dir_all(i, 0, direction_t::V);

		if (row_ranges[i].first)
			curr_row[row_ranges[i].first - 1].D = curr_row[row_ranges[i].first - 1].H = curr_row[row_ranges[i].first - 1].V = -infty;

		if (i < prof1_width)
		{
			if (i == 1)
				curr_row[0].V = max(prev_row[0].D, prev_row[0].V) + prof1_gap_term_open_curr * prof2_card;
			else
				curr_row[0].V = max(prev_row[0].D, prev_row[0].V) + prof1_gap_term_ext_curr * prof2_card;

			for (int j = row_ranges[i].second + 1; j <= min(row_ranges[i + 1].second, (int)prof2_width); ++j)
				curr_row[j].D = curr_row[j].H = curr_row[j].V = -infty;
		}
		else
			curr_row[0].V = -infty;

		// Calculate frequency of symbols in (i-1)-th column of profile1
/*		col1.clear();
		size_t col1_n_non_gaps = 0;
		for (size_t k = 0; k < NO_AMINOACIDS_AND_GAPS; ++k)
			if (profile1->counters.get_value(i, k))
			{
				size_t count = profile1->counters.get_value(i, k);
				col1.emplace_back(k, count);
				if (k < NO_AMINOACIDS)
					col1_n_non_gaps += count;
			}
		size_t col1_size = col1.size();
		*/

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

		// Precomputing for gap correction in the profile1
		n_gap_prof1_start_open = n_gap_prof1_start_ext = n_gap_prof1_start_term_open = n_gap_prof1_start_term_ext = 0;
		DP_SolveGapsProblemWhenStarting(i, prof1_width, prof1_card, profile1,
			n_gap_prof1_start_open, n_gap_prof1_start_ext, n_gap_prof1_start_term_open, n_gap_prof1_start_term_ext);

		DP_SolveGapsProblemWhenContinuing(i, prof1_width, prof1_card, profile1,
			n_gap_prof1_cont_ext, n_gap_prof1_cont_term_ext);

#ifndef NO_GAP_CORRECTION
		size_t n_gaps_prof1_to_change = profile1->counters.get_value(i, GAP_OPEN);   //the number of gaps to be changed from gap_open into gap_ext
		size_t n_gaps_prof1_term_to_change = profile1->counters.get_value(i, GAP_TERM_OPEN); //the number of gaps to be changed from term_open into term_ext
#else
		size_t n_gaps_prof1_to_change = 0;
		size_t n_gaps_prof1_term_to_change = 0;
#endif

		unsigned char* matrix_cell = matrix.get_cell(i, max(1, row_ranges[i].first) - 1);
		size_t max_j = min(row_ranges[i].second, (int)prof2_width);
		for (size_t j = max(1, row_ranges[i].first); j <= max_j; ++j)
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
			case 2:
				t += col1[1].second * scores2_column[col1[1].first];
			case 1:
				t += col1[0].second * scores2_column[col1[0].first];
			case 0:
				break;
			default:
				for (size_t k = 0; k < col1_size; ++k)
					t += col1[k].second * scores2_column[col1[k].first];
			}

			//an alignment of a column with another column, after an an alignment of two columns  
			score_t t_D = prev_row[j - 1].D + t;

			//..........................................................................
			//an alignment of a column with another column, after inserting a column of gaps before, into the profile 1
			score_t t_H = prev_row[j - 1].H;

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
			score_t t_V = prev_row[j - 1].V;

			t_V += t + gaps_prof2_change[j] * col1_n_non_gaps;

			if (t_D > t_H && t_D > t_V)
			{
				curr_row[j].D = t_D;
				matrix.set_dir_D(matrix_cell, direction_t::D);
			}
			else if (t_H > t_V)
			{
				curr_row[j].D = t_H;
				matrix.set_dir_D(matrix_cell, direction_t::H);
			}
			else
			{
				curr_row[j].D = t_V;
				matrix.set_dir_D(matrix_cell, direction_t::V);
			}

			//...........................................................................
			// Calcualte score for H
			//...........................................................................
			//inserting a columns of gaps into the profile 1

			//the case when a column (i) was aligned with the column (j-1)
			score_t gap_corr = prof2_gaps[j].open * n_gap_prof1_start_open + prof2_gaps[j].ext * n_gap_prof1_start_ext +
				prof2_gaps[j].term_open * n_gap_prof1_start_term_open + prof2_gaps[j].term_ext * n_gap_prof1_start_term_ext;

			t_D = curr_row[j - 1].D + gap_corr;

			//the case when a sequence of gaps is continued (another columnn of gaps is inserted)
			t_H = curr_row[j - 1].H + prof2_gaps[j].ext * n_gap_prof1_cont_ext +
				prof2_gaps[j].term_ext * n_gap_prof1_cont_term_ext;

#ifdef ALWAYS_3_DIRS
			//the column of gaps is inserted into the profile 1 after a column of gaps in the profile 2
			if ((i > 1) && (j > 1))
			{
				t_V = curr_row[j - 1].V + gap_corr;

				if (t_D > t_H && t_D > t_V)
				{
					curr_row[j].H = t_D;
					matrix.set_dir_H(matrix_cell, direction_t::D);
				}
				else if (t_V > t_H)	//Attention!!! Swapping the checking order 
				{
					curr_row[j].H = t_V;
					matrix.set_dir_H(matrix_cell, direction_t::V);
				}
				else
				{
					curr_row[j].H = t_H;
					matrix.set_dir_H(matrix_cell, direction_t::H);
				}
			}
			else
			{
				if (t_D > t_H)
				{
					curr_row[j].H = t_D;
					matrix.set_dir_H(matrix_cell, direction_t::D);
				}
				else
				{
					curr_row[j].H = t_H;
					matrix.set_dir_H(matrix_cell, direction_t::H);
				}
			}
#else // !ALWAYS_3_DIRS
			if (t_D > t_H)
			{
				curr_row[j].H = t_D;
				matrix.set_dir_H(matrix_cell, direction_t::D);
			}
			else
			{
				curr_row[j].H = t_H;
				matrix.set_dir_H(matrix_cell, direction_t::H);
			}
#endif // !ALWAYS_3_DIRS
			//...........................................................................
			// Calcualte score for V
			//...........................................................................
			//inserting a columns of gaps into the profile 2

			//the case when a column (i-1) was aligned with the column (j)
			gap_corr = prof1_gap_open_curr * gap_corrections[j].n_gap_start_open + prof1_gap_ext_curr * gap_corrections[j].n_gap_start_ext +
				prof1_gap_term_open_curr * gap_corrections[j].n_gap_start_term_open + prof1_gap_term_ext_curr * gap_corrections[j].n_gap_start_term_ext;
			t_D = prev_row[j].D + gap_corr;

			//the case when a sequence of gaps is continued (another columnn of gaps is inserted)
			t_V = prev_row[j].V + prof1_gap_ext_curr * gap_corrections[j].n_gap_cont_ext +
				prof1_gap_term_ext_curr * gap_corrections[j].n_gap_cont_term_ext;

#ifdef ALWAYS_3_DIRS
			if ((i > 1) && (j > 1))
			{
				t_H = prev_row[j].H + gap_corr;

				if (t_D > t_H && t_D > t_V)
				{
					curr_row[j].V = t_D;
					matrix.set_dir_V(matrix_cell, direction_t::D);
				}
				else if (t_H > t_V)
				{
					curr_row[j].V = t_H;
					matrix.set_dir_V(matrix_cell, direction_t::H);
				}
				else
				{
					curr_row[j].V = t_V;
					matrix.set_dir_V(matrix_cell, direction_t::V);
				}
			}
			else
			{
				if (t_D > t_V)
				{
					curr_row[j].V = t_D;
					matrix.set_dir_V(matrix_cell, direction_t::D);
				}
				else
				{
					curr_row[j].V = t_V;
					matrix.set_dir_V(matrix_cell, direction_t::V);
				}
			}
#else
			if (t_D > t_V)
			{
				curr_row[j].V = t_D;
				matrix.set_dir_V(matrix_cell, direction_t::D);
			}
			else
			{
				curr_row[j].V = t_V;
				matrix.set_dir_V(matrix_cell, direction_t::V);
			}
#endif //ALWAYS_3_DIRS

		}

		curr_row.swap(prev_row);
	}

	// Construct alignment
	ConstructProfile(profile1, profile2, matrix, prev_row.back());
}

// EOF
