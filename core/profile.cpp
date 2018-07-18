/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#include "../core/profile.h"
#include <assert.h>
#include <array>
#include <set>
#include <math.h>
#include <algorithm>
#include <stdio.h>

#define max3(x, y, z)	(max((x), max((y), (z))))

#define NINF(x)		((x) < -999 ? -999: (x))

int aacode[32] = {0, 20, 4, 3, 6, 13, 7, 8, 9, 23, 11, 10, 12, 2, 23, 14, 5, 1, 15, 16, 23, 19, 17, 22, 18, 21, 24, 25, 26, 27, 28, 29};
/* {a0,b1,c2,d3,e4,f5,g6,h7,i8,j-1,k9,l10,m11,n12,o23,p13,q14,r15,s16,t17,u17,v18,w19,x20,y21,z22};

 char CSequence::mapping_table[25] = "ARNDCQEGHILKMFPSTWYVBZX*";
                                      0123456789*123456789*123
									  */
// ****************************************************************************
// Construct empty profile
CProfile::CProfile(CParams *_params)
{
	params = _params;
}

// ****************************************************************************
// Convert sequence into a profile
CProfile::CProfile(const CGappedSequence &gapped_sequence, CParams *_params)
{
	params = _params;

	AppendRawSequence(gapped_sequence);
	CalculateCountersScores();
}


// ****************************************************************************
// Copy constructor
CProfile::CProfile(const CProfile &profile)
{
	data     = profile.data;
	counters = profile.counters;
	scores   = profile.scores;
	params   = profile.params;
	width    = profile.width;
	total_score = profile.total_score;
}


// ****************************************************************************
CProfile::CProfile(CProfile *profile1, CProfile *profile2, CParams *_params)
{
	params = _params;

	Align(profile1, profile2);
}


// ****************************************************************************
CProfile::~CProfile()
{
	for (auto &p : data)
		delete p;
}

// ****************************************************************************
bool CProfile::operator==(const CProfile &profile)
{
	if (data.size() != profile.data.size())
		return false;
	if (width != profile.width)
		return false;

	for (size_t i = 0; i < data.size(); ++i)
		if (*(data[i]) != *(profile.data[i]))
			return false;

	return true;
}

// ****************************************************************************
bool CProfile::operator!=(const CProfile &profile)
{
	return !(*this == profile);
}

// ****************************************************************************
void CProfile::CalculateCounters(CGappedSequence *gs)
{
	size_t first_non_gap;
	size_t last_non_gap;

	// Find first non-GAP symbol
	first_non_gap = gs->n_gaps[0] + 1;			// no. of initial gaps + GUARD

	// Find last non-GAP symbol
	last_non_gap = gs->gapped_size - gs->n_gaps[gs->size];			// gapped sequence size - no. of ending gaps

	// Set counters for terminal gaps at front
	if(first_non_gap > 1)
	{
		counters.add_value(1, GAP_TERM_OPEN, 1);

		for(size_t i = 2; i < first_non_gap; ++i)
			counters.add_value(i, GAP_TERM_EXT, 1);
	}
	
	// Set counters for terminal gaps at back
	if (last_non_gap < width)
	{
		counters.add_value(last_non_gap + 1, GAP_TERM_OPEN, 1);

		for (size_t i = width; i > (last_non_gap + 1); --i)
			counters.add_value(i, GAP_TERM_EXT, 1);
	}

	// Set counters for symbols and internal gaps
	bool was_gap = false;
	size_t seq_pos = first_non_gap;

	auto symbols = gs->symbols;
	auto n_gaps = gs->n_gaps;

	for(size_t i = 1; i < gs->size; ++i)
	{
		// Insert symbol
		counters.add_value(seq_pos, symbols[i], 1);
		++seq_pos;

		// Insert gaps after the symbol
		if(n_gaps[i] > 0)
		{
			// Increment gap open and gap close counters
			counters.add_value(seq_pos, GAP_OPEN, 1);
			//counters.add_value(seq_pos+gs->n_gaps[i]-1, GAP_OPEN, 1);

			// Increment gap ext counters
			for(int j = 1; j < n_gaps[i]; ++j)		
				counters.add_value(seq_pos+j, GAP_EXT, 1);
		}
		seq_pos += n_gaps[i];
	}

	// Insert last symbol in a sequence
	counters.add_value(seq_pos, gs->symbols[gs->size], 1);
	++seq_pos;
}

// ****************************************************************************
void CProfile::CalculateScores()
{
	score_t gap_open	   = params->gap_open;
	score_t gap_ext		   = params->gap_ext;
	score_t gap_term_open  = params->gap_term_open;
	score_t gap_term_ext   = params->gap_term_ext;

	size_t prof_size = data.size();

	scores.add_value(0, GAP_OPEN     , prof_size * gap_open);
	scores.add_value(0, GAP_EXT      , prof_size * gap_ext);
	scores.add_value(0, GAP_TERM_EXT , prof_size * gap_term_ext);
	scores.add_value(0, GAP_TERM_OPEN, prof_size * gap_term_open);

	for (size_t i = 1; i <= width; ++i)
	{
		size_t n_gap_open      = counters.get_value(i, GAP_OPEN);
		size_t n_gap_term_open = counters.get_value(i, GAP_TERM_OPEN);
		size_t n_gap_ext       = counters.get_value(i, GAP_EXT);
		size_t n_gap_term_ext  = counters.get_value(i, GAP_TERM_EXT);

		// Increment gap scores
		score_t gap_cost = 0;

		if(n_gap_open)
			gap_cost += n_gap_open * gap_open;
		if(n_gap_term_open)
			gap_cost += n_gap_term_open * gap_term_open;
		if(n_gap_ext)
			gap_cost += n_gap_ext * gap_ext;
		if(n_gap_term_ext)
			gap_cost += n_gap_term_ext * gap_term_ext;

		scores.add_value_to_column_part(i, NO_AMINOACIDS, gap_cost);

		// Increment symbol scores
		size_t tot_n_sym = 0;
		for(symbol_t sym = 0; sym < NO_AMINOACIDS; ++sym)
		{
			size_t n_sym = counters.get_value(i, sym);
			if(n_sym)
			{
				scores.add_column_part_mult(i, NO_AMINOACIDS, params->score_matrix[sym], n_sym);
				tot_n_sym += n_sym;
			}
		}

		scores.add_value(i, GAP_OPEN     , tot_n_sym * gap_open);
		scores.add_value(i, GAP_TERM_OPEN, tot_n_sym * gap_term_open);
		scores.add_value(i, GAP_EXT      , tot_n_sym * gap_ext);
		scores.add_value(i, GAP_TERM_EXT , tot_n_sym * gap_term_ext);
	}
}

// ****************************************************************************
void CProfile::CalculateCountersScores()
{
	if(data.empty())
		return;

	// Calculate counters
	counters.resize(data.front()->gapped_size+1);
	counters.set_zeros();
	for(auto &p : data)
		CalculateCounters(p);

	// Calculate scores
	scores.resize(data.front()->gapped_size+1);
	scores.set_zeros();
	CalculateScores();
}

// ****************************************************************************
// Perform alignment of two profiles
void CProfile::Align(CProfile *profile1, CProfile *profile2, vector<int> *column_mapping1, vector<int> *column_mapping2)
{
	if(profile1->data.size() == 0 || profile2->data.size() == 0)
	{
		assert(0);	// Profiles cannot be empty
	}

	if(profile1->counters.empty())
		profile1->CalculateCountersScores();
	if(profile2->counters.empty())
		profile2->CalculateCountersScores();

	if(profile1->data.size() == 1 && profile2->data.size() == 1)
		AlignSeqSeq(profile1, profile2);
	else if(profile1->data.size() == 1)
	{
		AlignSeqProf(profile1, profile2, column_mapping1, column_mapping2);
	}
	else if (profile2->data.size() == 1)
	{
		AlignSeqProf(profile2, profile1, column_mapping2, column_mapping1);
	}
	else
	{
		// Count the number of non-zero cells in profile and pass for the alignment the profile with fewer no. of non-zeros as the first one
		if (profile1->counters.get_num_of_non_zeros() * profile2->width < profile2->counters.get_num_of_non_zeros() * profile1->width)
			AlignProfProf(profile1, profile2, column_mapping1, column_mapping2);
		else
			AlignProfProf(profile2, profile1, column_mapping2, column_mapping1);
	}
}

// ****************************************************************************
void CProfile::Clear()
{
	data.clear();
	counters.clear();
	scores.clear();
}

// ****************************************************************************
bool CProfile::Empty()
{
	return data.empty();
}

// ****************************************************************************
size_t CProfile::Size()
{
	return data.size();
}

// ****************************************************************************
// Just append sequence in raw format (no alignment!)
void CProfile::AppendRawSequence(const CGappedSequence &gs)
{
	// Check that the sequence is of correct length
	if(data.empty())
	{
		if(gs.symbols.front() != GUARD)
			width = gs.gapped_size;
		else
			width = gs.gapped_size - 1;
		counters.resize(width + 1);
		counters.set_zeros();
		scores.resize(width + 1);
		scores.set_zeros();
	}
	else
	{
		size_t curr_size;
		if(gs.symbols.front() != GUARD)
			curr_size = gs.gapped_size;
		else
			curr_size = gs.gapped_size - 1;
		if (curr_size != width)
			assert(0);
	}

	data.push_back(new CGappedSequence(gs));
	if(gs.symbols.front() != GUARD)
		data.back()->symbols.insert(data.back()->symbols.begin(), GUARD);
}

// ****************************************************************************
void CProfile::GetGapStats(vector<size_t> &stats)
{
	stats.clear();
	stats.resize(width + 1, data.size());

	vector<size_t> stats_term;
	
	for(auto &p : data)
	{
		size_t seq_pos = 0;
		for(size_t i = 0; i <= p->size; ++i)
		{
			--stats[seq_pos];
			seq_pos += p->n_gaps[i] + 1;
		}
	}
}

// ****************************************************************************
// Remove empty columns
bool CProfile::Condense(vector<int> &column_mapping)
{
	bool r = false;

	if(data.empty())
		return r;

	width = data.front()->gapped_size;

	size_t card = data.size();

	vector<size_t> gap_stats;

	// Find non empty columns
	GetGapStats(gap_stats);
	
	vector<pair<uint32_t, uint32_t>> cols_to_remove;

	for(uint32_t i = 1; i <= width; ++i)
		if(gap_stats[i] == card)
		{
			if(gap_stats[i-1] == card)
				cols_to_remove.back().second++;
			else 
				cols_to_remove.push_back(make_pair(i, 1));
		}

	set<int> tmp_set;
	for(int i = 0; i <= (int) width; ++i)
		tmp_set.insert(i);
	for(auto &x: cols_to_remove)
		for(uint32_t i = x.first; i < x.first+x.second; ++i)
			tmp_set.erase(i);
	column_mapping.assign(tmp_set.begin(), tmp_set.end());

	if(!cols_to_remove.empty())
		r = true;

	while(!cols_to_remove.empty())
	{
		auto x = cols_to_remove.back();
		cols_to_remove.pop_back();

		for(auto &p: data)
			p->RemoveGaps(x.first, x.second);
	}
	
	width = data.front()->gapped_size;

	CalculateCountersScores();

	return r;
}

// ****************************************************************************
//
bool CProfile::OptimizeGaps()
{
	const int no_gap = 0;
	const int gap    = 1;

	// Do not optimize gaps if unnecessary!
	if(!params->enable_gap_optimization)
		return false;

	bool r = false;

	// Allocate space for decoded sequences (to speed up the analysis of exchangable columns)
	symbol_t *raw_gap_sequences = new symbol_t[data.size() * (width+1)];
	symbol_t **trans_gap_sequences = new symbol_t*[width+1];

	uint32_t data_size = (uint32_t) data.size();

	// Decode sequences and find boundariers between exchangable columns
	vector<int> exch_cols;
	vector<int> gap_sym_broundaries(width+1, 1);

	fill_n(raw_gap_sequences, data_size * (width+1), gap);

	for(int i = 0; i < (int) width+1; ++i)
		trans_gap_sequences[i] = &raw_gap_sequences[i*data_size];

	for(uint32_t i = 0; i < data_size; ++i)
	{
		auto n_gaps = data[i]->n_gaps;
		auto size = data[i]->size;

		uint32_t seq_pos = 1;

		// Starting gaps
		seq_pos += n_gaps[0];

		for(int j = 1; j <= (int) size; ++j)
		{
			trans_gap_sequences[seq_pos++][i] = no_gap;

			if(n_gaps[j] == 0)
				gap_sym_broundaries[seq_pos-1] = 0;

			seq_pos += n_gaps[j];
		}
	}

	for(int i = 1; i < (int) width; ++i)
		if(gap_sym_broundaries[i])
			exch_cols.push_back(i);

	if(exch_cols.empty())
	{
		delete[] raw_gap_sequences;
		delete[] trans_gap_sequences;

		return r;
	}

	// Find exchangable regions around possible boundaries
	vector<pair<int, int>> exch_ranges;
	int max_possible_col_id = width;					// do not analyze columns beyond this value (they were exchanged before)

	enum class exch_t {left, right, both, both_not_extend, both_to_left, both_to_right};
	vector<exch_t> exch_type;
	exch_type.reserve(data_size);

	while(!exch_cols.empty())
	{
		int curr_col_id = exch_cols.back();
		exch_cols.pop_back();

		if(curr_col_id + 2 > max_possible_col_id)
			continue;

		// Clear statistis of possible exchangable ranges
		exch_ranges.assign(data_size, make_pair(0, 0));
		int left_side, right_side;

		exch_type.clear();

		for(int j = 0; j < (int) data_size; ++j)
			if((trans_gap_sequences[curr_col_id][j]) && trans_gap_sequences[curr_col_id+1][j])
			{
				exch_type.push_back(exch_t::both);
				exch_ranges[j] = make_pair(1, 1);
			}
			else if(trans_gap_sequences[curr_col_id][j])
			{
				exch_type.push_back(exch_t::left);
				exch_ranges[j].first = 1;
			}
			else if(trans_gap_sequences[curr_col_id+1][j])
			{
				exch_type.push_back(exch_t::right);
				exch_ranges[j].second = 1;
			}

		// Find max. exchangable region to left
		bool exch_stop = false;
		for(left_side = 2; (curr_col_id + 1) - left_side > 0; ++left_side)
		{
			int examined_col_id = curr_col_id+1-left_side;
			auto trans_gap_sequences_col = trans_gap_sequences[examined_col_id];

			for(int j = 0; j < (int) data_size; ++j)
			{
				if(exch_type[j] == exch_t::left)
				{
					if(trans_gap_sequences_col[j])
						exch_ranges[j].first = left_side;
					else
						exch_stop = true;
				}
				else if(exch_type[j] == exch_t::both)
				{
					if(trans_gap_sequences_col[j])
						exch_ranges[j].first = left_side;
					else
						exch_type[j] = exch_t::both_not_extend;
				}
			}

			if(exch_stop)
				break;
		}

		// Find max. exchangable region to right
		for(int j = 0; j < (int) data_size; ++j)
			if(exch_type[j] == exch_t::both_not_extend)
				exch_type[j] = exch_t::both;

		exch_stop = false;
		for(right_side = 2; curr_col_id + right_side < max_possible_col_id; ++right_side)
		{
			int examined_col_id = curr_col_id+right_side;
			auto trans_gap_sequences_col = trans_gap_sequences[examined_col_id];

			for(int j = 0; j < (int) data_size; ++j)
			{
				if(exch_type[j] == exch_t::right)
				{
					if(trans_gap_sequences_col[j])
						exch_ranges[j].second = right_side;
					else
						exch_stop = true;
				}
				else if(exch_type[j] == exch_t::both)
				{
					if(trans_gap_sequences_col[j])
						exch_ranges[j].second = right_side;
					else
						exch_type[j] = exch_t::both_not_extend;
				}
			}

			if(exch_stop)
				break;
		}

		// Determine type of sequences of type 'both'
		for(int j = 0; j < (int) data_size; ++j)
			if(exch_type[j] == exch_t::both_not_extend)
				exch_type[j] = exch_t::both;

		bool is_exch_profitable = true;
		for(int j = 0; j < (int) data_size; ++j)
			if(exch_type[j] == exch_t::both)
			{
				if(exch_ranges[j].first < left_side-1 && exch_ranges[j].second < right_side-1)
					is_exch_profitable = false;
				else if(exch_ranges[j].first >= left_side-1 && exch_ranges[j].second >= right_side-1)
					;	// do nothing - gap can be used as left or right
				else if(exch_ranges[j].first >= left_side-1)
					exch_type[j] = exch_t::both_to_left;
				else
					exch_type[j] = exch_t::both_to_right;
			}

		if(!is_exch_profitable)
			continue;

		if(curr_col_id - left_side <= 0 || curr_col_id + right_side >= (int) width)	// do not exchange columns at profile ends
			continue;

		// Determine whether exchanging coolumns would be profitable (would reduce the number of gap opens)
		int n_gap_open_balance = 0;

		int left_col = curr_col_id - (left_side - 1);
		int right_col = curr_col_id + right_side;

		auto trans_gap_sequences_left  = trans_gap_sequences[left_col];
		auto trans_gap_sequences_right = trans_gap_sequences[right_col];

		for(int j = 0; j < (int) data_size; ++j)
		{
			if(exch_type[j] == exch_t::left)
			{
				if(!trans_gap_sequences_left[j])
					n_gap_open_balance--;		// left side gap open would be removed
				if(!trans_gap_sequences_right[j])
					n_gap_open_balance++;		// gap would be added at right side
			}
			else if(exch_type[j] == exch_t::right)
			{
				if(!trans_gap_sequences_right[j])
					n_gap_open_balance--;		// right side gap open would be removed
				if(!trans_gap_sequences_left[j])
					n_gap_open_balance++;		// gap would be added at left side
			}
			else if(exch_type[j] == exch_t::both_to_left)
			{
				if(!trans_gap_sequences_right[j])
					n_gap_open_balance++;		// gap would be added at right side
			}
			else if(exch_type[j] == exch_t::both_to_right)
			{
				if(!trans_gap_sequences_left[j])
					n_gap_open_balance++;		// gap would be added at left side
			}
		}

		if(n_gap_open_balance < 0)
		{
			// Exchange columns
			for(int j = 0; j < (int) data_size; ++j)
			{
				if(exch_type[j] == exch_t::left || exch_type[j] == exch_t::both_to_left)
				{
					for(int k = 0; k < left_side-1; ++k)
						data[j]->InsertGap(right_col);
					for(int k = 0; k < left_side-1; ++k)
						data[j]->RemoveGap(left_col+1);
				}
				else if(exch_type[j] == exch_t::right || exch_type[j] == exch_t::both_to_right)
				{
					for(int k = 0; k < right_side-1; ++k)
						data[j]->RemoveGap(curr_col_id+1);
					for(int k = 0; k < right_side-1; ++k)
						data[j]->InsertGap(left_col+1);
				}
			}

			max_possible_col_id = left_col;

			r = true;
		}
	}

	delete[] trans_gap_sequences;
	delete[] raw_gap_sequences;

	return r;
}

// ****************************************************************************
// 
void CProfile::AlignProfProf(CProfile *profile1, CProfile *profile2, vector<int> *column_mapping1, vector<int> *column_mapping2)
{
	size_t prof1_width = profile1->width;
	size_t prof2_width = profile2->width;

	size_t prof1_card = profile1->data.size();
	size_t prof2_card = profile2->data.size();

	score_t gap_open      = params->gap_open;
	score_t gap_ext       = params->gap_ext;
	score_t gap_term_open = params->gap_term_open;
	score_t gap_term_ext  = params->gap_term_ext;

	CDPMatrix matrix(prof1_width + 1, prof2_width + 1);
	matrix.set_zeros();

	dp_row_t curr_row(prof2_width + 1);
	dp_row_t prev_row(prof2_width + 1);

	// Frequency of symbols at all columns of profile1 and profile2
	vector<array<int, 32>> freq1(prof1_width + 1, array<int, 32>());
	vector<pair<int, int>> col1(32);

	CProfileValues<score_t, NO_SYMBOLS> &scores1 = const_cast<CProfileValues<score_t, NO_SYMBOLS>&>(profile1->scores);
	CProfileValues<score_t, NO_SYMBOLS> &scores2 = const_cast<CProfileValues<score_t, NO_SYMBOLS>&>(profile2->scores);

	// Prepare ranges for computed cells in rows
	bool is_guided = column_mapping1 != nullptr && column_mapping2 != nullptr;
	vector<pair<int, int>> row_ranges;

	if(is_guided)
		FindRowRanges(column_mapping1, column_mapping2, row_ranges);
	else
		for(int i = 0; i <= (int) prof1_width; ++i)
			row_ranges.push_back(make_pair(0, prof2_width));

	// Precompute scores for gaps for profile2
	vector<dp_gap_costs> prof2_gaps(prof2_width+1);

	size_t n_gap_open, n_gap_ext, n_gap_term_open, n_gap_term_ext;
	       n_gap_open = n_gap_ext = n_gap_term_open = n_gap_term_ext = 0;

	size_t n_gap_prof1_start_open, n_gap_prof1_start_ext, n_gap_prof1_start_term_open, n_gap_prof1_start_term_ext,
		n_gap_prof1_cont_ext, n_gap_prof1_cont_term_ext;
	vector<dp_gap_corrections> gap_corrections(prof2_width+1);

	//scoresX.get_value(j, GAP_OPEN) returns information of how the column j of the profile X aligns with a single gap_open
	//scoresX.get_value(j, GAP_EXT) returns information of how the column j of the profile X aligns with a single gap_ext
	//etc.
	//The task is to calculate how the column j of the profile X aligns with a column of gaps of different categories
	
	for (size_t j = 0; j <= prof2_width; ++j)
	{
		prof2_gaps[j].open      = scores2.get_value(j, GAP_OPEN); 
		prof2_gaps[j].ext       = scores2.get_value(j, GAP_EXT);
		prof2_gaps[j].term_open = scores2.get_value(j, GAP_TERM_OPEN);
		prof2_gaps[j].term_ext  = scores2.get_value(j, GAP_TERM_EXT);
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
		prev_row[j].D  = -infty;
		
		prev_row[j].H = prev_row[j - 1].H + prof2_gaps[j].term_ext * prof1_card;

		prev_row[j].V  = -infty;
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
		n_gaps_prof2_to_change[j]      = profile2->counters.get_value(j, GAP_OPEN);      //the number of gaps to be changed from gap_open into gap_ext 
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
		score_t prof1_gap_open_prev = scores1.get_value(i - 1, GAP_OPEN);
		score_t prof1_gap_open_curr = scores1.get_value(i, GAP_OPEN);
		score_t prof1_gap_term_open_prev = scores1.get_value(i - 1, GAP_TERM_OPEN);
		score_t prof1_gap_term_open_curr = scores1.get_value(i, GAP_TERM_OPEN);
		score_t prof1_gap_ext_curr = scores1.get_value(i, GAP_EXT);
		score_t prof1_gap_term_ext_curr = scores1.get_value(i, GAP_TERM_EXT);

		// Boundary conditions
		curr_row[0].D  = -infty;
		curr_row[0].H  = -infty;
		matrix.set_dir_all(i, 0, direction_t::V);

		if(row_ranges[i].first)
			curr_row[row_ranges[i].first-1].D = curr_row[row_ranges[i].first-1].H = curr_row[row_ranges[i].first-1].V = -infty;

		if (i < prof1_width)
		{
			if (i == 1)
				curr_row[0].V = max(prev_row[0].D, prev_row[0].V) + prof1_gap_term_open_curr * prof2_card;
			else
				curr_row[0].V = max(prev_row[0].D, prev_row[0].V) + prof1_gap_term_ext_curr * prof2_card;

			for(int j = row_ranges[i].second+1; j <= min(row_ranges[i+1].second, (int) prof2_width); ++j)
				curr_row[j].D = curr_row[j].H = curr_row[j].V = -infty;
		}
		else
			curr_row[0].V = -infty;

		// Calculate frequency of symbols in (i-1)-th column of profile1
		col1.clear();
		size_t col1_n_non_gaps = 0;
		for(size_t k = 0; k < NO_AMINOACIDS_AND_GAPS; ++k)			
			if(profile1->counters.get_value(i, k))
			{
				size_t count = profile1->counters.get_value(i, k);
				col1.push_back(make_pair(k, count));
				if(k < NO_AMINOACIDS)
					col1_n_non_gaps += count;
			}
		size_t col1_size = col1.size();

		// Precomputing for gap correction in the profile1
		n_gap_prof1_start_open = n_gap_prof1_start_ext = n_gap_prof1_start_term_open = n_gap_prof1_start_term_ext = 0;
		DP_SolveGapsProblemWhenStarting(i , prof1_width, prof1_card, profile1, 
			n_gap_prof1_start_open, n_gap_prof1_start_ext, n_gap_prof1_start_term_open, n_gap_prof1_start_term_ext);

		DP_SolveGapsProblemWhenContinuing(i, prof1_width, prof1_card, profile1,
			n_gap_prof1_cont_ext, n_gap_prof1_cont_term_ext);

#ifndef NO_GAP_CORRECTION
		size_t n_gaps_prof1_to_change      = profile1->counters.get_value(i, GAP_OPEN);   //the number of gaps to be changed from gap_open into gap_ext
		size_t n_gaps_prof1_term_to_change = profile1->counters.get_value(i, GAP_TERM_OPEN); //the number of gaps to be changed from term_open into term_ext
#else
		size_t n_gaps_prof1_to_change = 0;
		size_t n_gaps_prof1_term_to_change = 0;
#endif

		unsigned char *matrix_cell = matrix.get_cell(i, max(1, row_ranges[i].first) - 1);
		size_t max_j = min(row_ranges[i].second, (int) prof2_width);
		for (size_t j = max(1, row_ranges[i].first); j <= max_j; ++j)
		{
			// Get current cell of matrix for faster access
			++matrix_cell;

			//...........................................................................
			// Calcualte score for D
			//...........................................................................

			//analyzing an alignment of the column j of the profile 2 with the column i of the profile 1
			score_t *scores2_column = scores2.get_column(j);

			//scores2_column contains information on how the column j aligns with a single character that could appear in the profile 1
			//the variable t is to store a cost of the alignment of the column j with a specific character (col1[k].first) 
			//of a concrete number of occurrences (col1[k].second) of the profile 1
			
			score_t t = 0;
			switch(col1_size)
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
				for(size_t k = 0; k < col1_size; ++k)
					t += col1[k].second * scores2_column[col1[k].first];
			}

			//an alignment of a column with another column, after an an alignment of two columns  
			score_t t_D = prev_row[j-1].D + t;
            
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
				score_t delta_t = n_gaps_prof1_to_change      * (scores2_column[GAP_EXT] - scores2_column[GAP_OPEN]) +
					n_gaps_prof1_term_to_change * (scores2_column[GAP_TERM_EXT] - scores2_column[GAP_TERM_OPEN]);
					
				t_H = t_H + t + delta_t;
			}
			else
				t_H = t_H + t;

			//..........................................................................
			//an alignment of a column with another column, after inserting a column of gaps before, into the profile 2
			score_t t_V = prev_row[j - 1].V;

			t_V += t + gaps_prof2_change[j] * col1_n_non_gaps;

			if(t_D > t_H && t_D > t_V)
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
			score_t gap_corr = prof2_gaps[j].open * n_gap_prof1_start_open           + prof2_gaps[j].ext * n_gap_prof1_start_ext +
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
				if(t_D > t_H)
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
				if(t_D > t_H)
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
				if(t_D > t_V)
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
				if(t_D > t_V)
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

// ****************************************************************************
// 
void CProfile::AlignSeqSeq(CProfile *profile1, CProfile *profile2)
{
	size_t prof1_width = profile1->width;
	size_t prof2_width = profile2->width;

	// These profiles contain single sequences so there are no gaps here
	vector<symbol_t> &seq1 = const_cast<vector<symbol_t>&>(profile1->data[0]->symbols);
	vector<symbol_t> &seq2 = const_cast<vector<symbol_t>&>(profile2->data[0]->symbols);

	score_t gap_open = params->gap_open;
	score_t gap_ext = params->gap_ext;
	score_t gap_term_open = params->gap_term_open;
	score_t gap_term_ext = params->gap_term_ext;

	CDPMatrix matrix(prof1_width + 1, prof2_width + 1);
	matrix.set_zeros();

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
				curr_row[j].D = t_D + params->score_matrix[seq1[i]][seq2[j]];
				matrix.set_dir_D(i, j, direction_t::D);
			}
			else if (t_H >= t_V)
			{
				curr_row[j].D = t_H + params->score_matrix[seq1[i]][seq2[j]];
				matrix.set_dir_D(i, j, direction_t::H);
				//???
			}
			else
			{
				curr_row[j].D = t_V + params->score_matrix[seq1[i]][seq2[j]];
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
//
void CProfile::ConstructProfile(CProfile *profile1, CProfile *profile2, CDPMatrix &matrix, dp_row_elem_t &last_elem)
{
	score_t gap_open	  = params->gap_open;
 	score_t gap_ext	      = params->gap_ext;
	score_t gap_term_open = params->gap_term_open;
	score_t gap_term_ext  = params->gap_term_ext;

	size_t prof1_width = profile1->width;
	size_t prof2_width = profile2->width;

	size_t prof1_size = profile1->data.size();
	size_t prof2_size = profile2->data.size();

	CProfileValues<counter_t, NO_SYMBOLS> *counters1 = &(profile1->counters);
	CProfileValues<counter_t, NO_SYMBOLS> *counters2 = &(profile2->counters);

	CProfileValues<score_t, NO_SYMBOLS> *scores1 = &(profile1->scores);
	CProfileValues<score_t, NO_SYMBOLS> *scores2 = &(profile2->scores);

	data.clear();

	size_t i = prof1_width;
	size_t j = prof2_width;

	direction_t dir, prev_dir, next_dir, init_dir;

	size_t n_gap_open, n_gap_ext, n_gap_term_open, n_gap_term_ext;
	size_t n_gap_to_transfer1, n_gap_term_to_transfer1, n_gap_to_transfer2, n_gap_term_to_transfer2;
	size_t n_gap_open_at_left1, n_gap_ext_at_left1, n_gap_term_open_at_left1, n_gap_term_ext_at_left1;
	size_t n_gap_open_at_left2, n_gap_ext_at_left2, n_gap_term_open_at_left2, n_gap_term_ext_at_left2;

    if(last_elem.D >= last_elem.H && last_elem.D >= last_elem.V)
	{
		dir = direction_t::D;
		total_score = last_elem.D;
	}
	else if(last_elem.H > last_elem.V)
	{
		dir = direction_t::H;
		total_score = last_elem.H;
	}
	else
	{
		dir = direction_t::V;
		total_score = last_elem.V;
	}
	
	prev_dir = direction_t::D;
	init_dir = dir;

	vector<direction_t> path;
	path.reserve(prof1_width + prof2_width);
	size_t path_next_index = 1;
	path.push_back(init_dir);

	// Just find the size of the new profile
	width = 0;
	while (i || j)
	{
		width++;
		if (dir == direction_t::D)
		{
			dir = matrix.get_dir_D(i--, j--);
			path.push_back(dir);
		}
		else if (dir == direction_t::H)
		{
			dir = matrix.get_dir_H(i, j--);
			path.push_back(dir);
		}
		else if (dir == direction_t::V)
		{
			dir = matrix.get_dir_V(i--, j);
			path.push_back(dir);
		}
		else
		{
			assert(0);
		}
	}

	//The last entry in the "path" is to be cut off. The variable "width" is incremented at the entrance to the loop and stores the correct path length, which is without the last entry.
	//The last entry in the "path" is made inside the loop and it means the entrance to the boundary row/column. The two indices i and j are equal to 0 .
	//Tested!
	//After making reverse() the last entry is on position 0. That means the "right" path starts at index 1
	
	reverse(path.begin(), path.end());
	
	// Construct the new profile
	i = prof1_width;
	j = prof2_width;

	prev_dir = direction_t::D; //This prevents recognition of continuing a gap sequence, when the first inserted symbol should be a gap
		
	dir = path[1];
	
	counters.resize(width + 1);
	counters.set_zeros();
	scores.resize(width + 1);
	scores.set_zeros();

	size_t new_prof_col_id = 1;
	i = 0; 
	j = 0;
	n_gap_to_transfer1 = n_gap_term_to_transfer1 = n_gap_to_transfer2 = n_gap_term_to_transfer2 = 0;
	n_gap_open_at_left1 = n_gap_ext_at_left1 = n_gap_term_open_at_left1 = n_gap_term_ext_at_left1 = 0;
	n_gap_open_at_left2 = n_gap_ext_at_left2 = n_gap_term_open_at_left2 = n_gap_term_ext_at_left2 = 0;

	path.push_back(direction_t::V); 
	
	// Traversing the DP matrix
	
	cumulate_gap_inserts = false;
	no_cumulated_gap_inserts = 0;

	for (path_next_index = 2; path_next_index < path.size(); path_next_index++)
	{
		if(dir == direction_t::D)
		{
			i++;
			j++;

			//if as a result of the insertion of a column of gaps on the left side of the current column
			//there is a need to transform gap_opens into gap_exts inside the current column. 
			//Understood as the transfer of information from the previous column
			if (n_gap_to_transfer1 || n_gap_term_to_transfer1)
			{
				counters1->add_value(i, GAP_EXT,       n_gap_to_transfer1);
				counters1->sub_value(i, GAP_OPEN,      n_gap_to_transfer1);
				counters1->add_value(i, GAP_TERM_EXT,  n_gap_term_to_transfer1);
				counters1->sub_value(i, GAP_TERM_OPEN, n_gap_term_to_transfer1);

				//"scores" should be changed together with "counters"
				score_t gap_cost = n_gap_to_transfer1 * (gap_ext - gap_open) +
								   n_gap_term_to_transfer1 * (gap_term_ext - gap_term_open);

				for (size_t sym = 0; sym < NO_AMINOACIDS; ++sym)
						scores1->add_value(i, sym, gap_cost);

				n_gap_to_transfer1 = n_gap_term_to_transfer1 = 0;
			}
			if (n_gap_to_transfer2 || n_gap_term_to_transfer2)
			{
				counters2->add_value(j, GAP_EXT, n_gap_to_transfer2);
				counters2->sub_value(j, GAP_OPEN, n_gap_to_transfer2);
				counters2->add_value(j, GAP_TERM_EXT, n_gap_term_to_transfer2);
				counters2->sub_value(j, GAP_TERM_OPEN, n_gap_term_to_transfer2);

				score_t gap_cost = n_gap_to_transfer2 * (gap_ext - gap_open) +
					n_gap_term_to_transfer2 * (gap_term_ext - gap_term_open);

				for (size_t sym = 0; sym < NO_AMINOACIDS; ++sym)
					scores2->add_value(j, sym, gap_cost);

				n_gap_to_transfer2 = n_gap_term_to_transfer2 = 0;
			}

			n_gap_open_at_left1 = n_gap_ext_at_left1 = n_gap_term_open_at_left1 = n_gap_term_ext_at_left1 = 0;
			n_gap_open_at_left2 = n_gap_ext_at_left2 = n_gap_term_open_at_left2 = n_gap_term_ext_at_left2 = 0;
			
			InsertColumn(new_prof_col_id, profile1, i);
			InsertColumn(new_prof_col_id, profile2, j);
			prev_dir = dir;
			dir = path[path_next_index];
		}
		else if(dir == direction_t::H)
		{
			next_dir = path[path_next_index];
			n_gap_open = n_gap_ext = n_gap_term_open = n_gap_term_ext = 0;

			if (prev_dir == dir)		// inside a gap
				SolveGapsProblemWhenContinuing(i, prof1_width, prof1_size, n_gap_to_transfer1, n_gap_term_to_transfer1, n_gap_open, n_gap_ext, n_gap_term_open, n_gap_term_ext, n_gap_open_at_left1, n_gap_ext_at_left1, n_gap_term_open_at_left1, n_gap_term_ext_at_left1);
			
			else //if (prev_dir != dir)							// start of a gap
				SolveGapsProblemWhenStarting(i, prof1_width, prof1_size, profile1, n_gap_to_transfer1, n_gap_term_to_transfer1, n_gap_open, n_gap_ext, n_gap_term_open, n_gap_term_ext);
			
			//the following information is to be used when extending the sequence of gaps.
			//It is then the need to trace what has been inserted in the previous step.
			n_gap_open_at_left1 = n_gap_open;
			n_gap_ext_at_left1 = n_gap_ext;
			n_gap_term_open_at_left1 = n_gap_term_open;
			n_gap_term_ext_at_left1 = n_gap_term_ext;

			prev_dir = dir;
			dir = next_dir;

			if(prev_dir == next_dir && path_next_index < path.size() - 1)
				cumulate_gap_inserts = true;
			else
				cumulate_gap_inserts = false;

			InsertGaps(new_prof_col_id, profile1, new_prof_col_id, n_gap_open, n_gap_ext, n_gap_term_open, n_gap_term_ext);

			//-----------------
			if (n_gap_to_transfer2 || n_gap_term_to_transfer2)
			{
				counters2->add_value(j + 1, GAP_EXT, n_gap_to_transfer2);
				counters2->sub_value(j + 1, GAP_OPEN, n_gap_to_transfer2);
				counters2->add_value(j + 1, GAP_TERM_EXT, n_gap_term_to_transfer2);
				counters2->sub_value(j + 1, GAP_TERM_OPEN, n_gap_term_to_transfer2);

				//"scores" should be changed together with "counters"
				score_t gap_cost = n_gap_to_transfer2 * (gap_ext - gap_open) +
					n_gap_term_to_transfer2 * (gap_term_ext - gap_term_open);

				for (size_t sym = 0; sym < NO_AMINOACIDS; ++sym)
					scores2->add_value(j + 1, sym, gap_cost);

				n_gap_to_transfer2 = n_gap_term_to_transfer2 = 0;
			}
			//-----------------
			InsertColumn(new_prof_col_id, profile2, ++j);
		}
		else if(dir == direction_t::V)
		{
			next_dir = path[path_next_index];
			n_gap_open = n_gap_ext = n_gap_term_open = n_gap_term_ext = 0;

			if (prev_dir == dir)		// inside a gap
				SolveGapsProblemWhenContinuing(j, prof2_width, prof2_size, n_gap_to_transfer2, n_gap_term_to_transfer2, n_gap_open, n_gap_ext, n_gap_term_open, n_gap_term_ext, n_gap_open_at_left2, n_gap_ext_at_left2, n_gap_term_open_at_left2, n_gap_term_ext_at_left2);

			else //if (prev_dir != dir)							// start of a gap
				SolveGapsProblemWhenStarting(j, prof2_width, prof2_size, profile2, n_gap_to_transfer2, n_gap_term_to_transfer2, n_gap_open, n_gap_ext, n_gap_term_open, n_gap_term_ext);

			//the following information is to be used when extending the sequence of gaps.
			//It is then the need to trace what has been inserted in the previous step.
			n_gap_open_at_left2 = n_gap_open;
			n_gap_ext_at_left2 = n_gap_ext;
			n_gap_term_open_at_left2 = n_gap_term_open;
			n_gap_term_ext_at_left2 = n_gap_term_ext;

			prev_dir = dir;
			dir = next_dir;

			if(prev_dir == next_dir && path_next_index < path.size() - 1)
				cumulate_gap_inserts = true;
			else
				cumulate_gap_inserts = false;

			//----------------
			if (n_gap_to_transfer1 || n_gap_term_to_transfer1)
			{
				counters1->add_value(i + 1, GAP_EXT, n_gap_to_transfer1);
				counters1->sub_value(i + 1, GAP_OPEN, n_gap_to_transfer1);
				counters1->add_value(i + 1, GAP_TERM_EXT, n_gap_term_to_transfer1);
				counters1->sub_value(i + 1, GAP_TERM_OPEN, n_gap_term_to_transfer1);

				//"scores" should be changed together with "counters"
				score_t gap_cost = n_gap_to_transfer1 * (gap_ext - gap_open) +
					n_gap_term_to_transfer1 * (gap_term_ext - gap_term_open);

				for (size_t sym = 0; sym < NO_AMINOACIDS; ++sym)
					scores1->add_value(i + 1, sym, gap_cost);

				n_gap_to_transfer1 = n_gap_term_to_transfer1 = 0;
			}
			//----------------
			
			InsertColumn(new_prof_col_id, profile1, ++i);			
			InsertGaps(new_prof_col_id, profile2, new_prof_col_id, n_gap_open, n_gap_ext, n_gap_term_open, n_gap_term_ext);
		}
		else
			assert(0);

		++new_prof_col_id;
	}

	// !!! TODO: consider moving instead of copying
	data.insert(data.end(), profile1->data.begin(), profile1->data.end());
	data.insert(data.end(), profile2->data.begin(), profile2->data.end());

	profile1->data.clear();
	profile2->data.clear();

	scores.set_value(0, GAP_OPEN     , gap_open      * data.size());
	scores.set_value(0, GAP_EXT      , gap_ext       * data.size());
	scores.set_value(0, GAP_TERM_OPEN, gap_term_open * data.size());
	scores.set_value(0, GAP_TERM_EXT , gap_term_ext  * data.size());
}

// ****************************************************************************
void CProfile::InsertGaps(size_t prof_col_id, CProfile *profile, size_t col_id, size_t n_gap_open, size_t n_gap_ext, size_t n_gap_term_open, size_t n_gap_term_ext)
{
	score_t gap_open_r	    = params->gap_open;
	score_t gap_ext_r	    = params->gap_ext;
	score_t gap_term_open_r	= params->gap_term_open;
	score_t gap_term_ext_r  = params->gap_term_ext;

	size_t num = profile->data.size();

	++no_cumulated_gap_inserts;
	if(!cumulate_gap_inserts)
	{
		for(size_t i = 0; i < num; ++i)
			profile->data[i]->InsertGaps(col_id + 1 - no_cumulated_gap_inserts, no_cumulated_gap_inserts);
		no_cumulated_gap_inserts = 0;
	}

	score_t gap_cost = n_gap_open * (gap_open_r) + 
					   n_gap_ext * (gap_ext_r) + 
					   n_gap_term_open * (gap_term_open_r) +
					   n_gap_term_ext * (gap_term_ext_r);

	counters.add_value(prof_col_id, GAP_OPEN, n_gap_open);
	counters.add_value(prof_col_id, GAP_EXT, n_gap_ext);
	counters.add_value(prof_col_id, GAP_TERM_OPEN, n_gap_term_open);
	counters.add_value(prof_col_id, GAP_TERM_EXT, n_gap_term_ext);

	counters.add_value(prof_col_id, GAP, num);				// total number of gaps in profile

	for(size_t i = 0; i < NO_AMINOACIDS; ++i)	
		scores.add_value(prof_col_id, i, gap_cost);
}

// ****************************************************************************
void CProfile::InsertColumn(size_t prof_col_id, CProfile *profile, size_t col_id)
{
	score_t gap_open	   = params->gap_open;
	score_t gap_ext		   = params->gap_ext;
	score_t gap_term_open  = params->gap_term_open;
	score_t gap_term_ext   = params->gap_term_ext;

	counters.add_column(prof_col_id, profile->counters.get_column(col_id));
	scores.add_column(prof_col_id, profile->scores.get_column(col_id));
}


// ****************************************************************************
void CProfile::SolveGapsProblemWhenContinuing(size_t source_col_id, size_t prof_width, size_t prof_size, size_t &n_gap_to_transfer, size_t &n_gap_term_to_transfer, size_t &n_gap_open, size_t &n_gap_ext, size_t &n_gap_term_open, size_t &n_gap_term_ext, size_t n_gap_open_at_left, size_t n_gap_ext_at_left, size_t n_gap_term_open_at_left, size_t n_gap_term_ext_at_left)
{
	// inside a gap
	if (source_col_id == prof_width || source_col_id == 0) 
			n_gap_term_ext += prof_size;
	else
	{
#ifndef NO_GAP_CORRECTION
			//the number of n_gap_term_ext corresponds to the number of GAP_TERM_OPEN in the column on the left 
		n_gap_term_ext += n_gap_term_open_at_left;

			//the number of n_gap_term_ext corresponds to the number of GAP_TERM_EXT in the column on the left
		n_gap_term_ext += n_gap_term_ext_at_left;

			//the number of n_gap_ext corresponds to the number of GAP_OPEN in the column on the left
		n_gap_ext = n_gap_open_at_left;

			//the number of n_gap_ext corresponds to the number of GAP_EXT in the column on the left
		n_gap_ext += n_gap_ext_at_left;

		n_gap_open = prof_size - n_gap_ext - n_gap_term_ext;

#else
		n_gap_ext = prof_size;
		n_gap_open = prof_size - n_gap_ext - n_gap_term_ext;		
#endif
	}

}

// ****************************************************************************
void CProfile::SolveGapsProblemWhenStarting(size_t source_col_id, size_t prof_width, size_t prof_size, CProfile *profile, size_t &n_gap_to_transfer, size_t &n_gap_term_to_transfer, size_t &n_gap_open, size_t &n_gap_ext, size_t &n_gap_term_open, size_t &n_gap_term_ext)
{
	if (source_col_id == 0)
	{	//the begining of a terminal gap. n_gap_term_open corresponds to the size of the profile
		n_gap_term_open += prof_size;

		//There can be the necessity of modifications in the column on the right. 
		//"open" should be turned into "ext"
		
#ifndef NO_GAP_CORRECTION
		n_gap_term_to_transfer = profile->counters.get_value(source_col_id + 1, GAP_TERM_OPEN);
#else
		n_gap_term_to_transfer = 0; 
#endif

	}
	else if (source_col_id >= prof_width)
	{
		//the begining of a terminal gap on right side of the profile
		//If - in the last column of the current profile - there were term_opens or term_exts
		//one should insert term_ext
		//Term_opens should be inserted to the rest of the sequences
		
#ifndef NO_GAP_CORRECTION
		counter_t cnt = profile->counters.get_value(source_col_id, GAP_TERM_OPEN) + profile->counters.get_value(source_col_id, GAP_TERM_EXT);
		n_gap_term_ext = cnt;

		n_gap_term_open += (prof_size - cnt);
#else
		n_gap_term_open = prof_size;
		n_gap_term_ext = 0;	
#endif
	}
	else	//start of a gap inside the profile
	{

#ifndef NO_GAP_CORRECTION
		//the number of n_gap_term_open corresponds to the number of GAP_TERM_OPEN in the column on the right, 
		n_gap_term_open += profile->counters.get_value(source_col_id + 1, GAP_TERM_OPEN);

		//GAP_TERM_OPEN from the right column should be turned into GAP_TERM_EXT
		n_gap_term_to_transfer = n_gap_term_open;

		//the number of n_gap_term_ext corresponds to GAP_TERM_OPEN in the left column, 
		n_gap_term_ext += profile->counters.get_value(source_col_id, GAP_TERM_OPEN);

		//the number of n_gap_term_ext corresponds to GAP_TERM_EXT in the left column,
		n_gap_term_ext += profile->counters.get_value(source_col_id, GAP_TERM_EXT);

		//the number of n_gap_ext corresponds GAP_OPEN in the left column,
		n_gap_ext = profile->counters.get_value(source_col_id, GAP_OPEN);

		//the number of n_gap_ext corresponds GAP_EXT in the left column,
		n_gap_ext += profile->counters.get_value(source_col_id, GAP_EXT);

		//the number of n_gap_open corresponds GAP_OPEN in the right column,
		n_gap_open = profile->counters.get_value(source_col_id + 1, GAP_OPEN);

		//GAP_OPEN from the right column should be turned into GAP_EXT
		n_gap_to_transfer += n_gap_open;

		//gap_open is to be inserted into the rest of sequences :)
		n_gap_open = prof_size - n_gap_ext - n_gap_term_open - n_gap_term_ext;

#else
		n_gap_open = prof_size;
		n_gap_term_to_transfer = 0;		
		n_gap_ext = 0;					
#endif
	
	}
}

// ****************************************************************************
void CProfile::DP_SolveGapsProblemWhenStarting(size_t source_col_id, size_t prof_width, size_t prof_size, CProfile *profile, size_t &n_gap_open, size_t &n_gap_ext, size_t &n_gap_term_open, size_t &n_gap_term_ext)
{
	//The situation below never happens in DP matrix
	//if (source_col_id == 0)
	//{	
	//	n_gap_term_open += prof_size;
	//}
	//else 
	if (source_col_id >= prof_width)
	{
		//the begining of a terminal gap on the right side of the profile
		//If - in the last column of the current profile - there were term_opens or term_exts
		//one should insert term_ext
		//Term_opens should be inserted to the rest of the sequences

#ifndef NO_GAP_CORRECTION
		counter_t *source_col = profile->counters.get_column(source_col_id);
		counter_t cnt = source_col[GAP_TERM_OPEN] + source_col[GAP_TERM_EXT];
		n_gap_term_ext = cnt;

		n_gap_term_open += (prof_size - cnt);
#else
		n_gap_term_open = prof_size;
		n_gap_term_ext = 0;				
#endif

	}
	else	//start of a gap inside the profile
	{

#ifndef NO_GAP_CORRECTION
		//the number of n_gap_term_open corresponds to the number of GAP_TERM_OPEN in the column on the right, 
		n_gap_term_open += profile->counters.get_value(source_col_id + 1, GAP_TERM_OPEN);

		counter_t *source_col = profile->counters.get_column(source_col_id);

		//the number of n_gap_term_ext corresponds to GAP_TERM_OPEN in the left column,  
		n_gap_term_ext += source_col[GAP_TERM_OPEN];

		//the number of n_gap_term_ext corresponds to GAP_TERM_EXT in the left column,
		n_gap_term_ext += source_col[GAP_TERM_EXT];

		//the number of n_gap_ext corresponds to GAP_OPEN in the left column, 
		n_gap_ext = source_col[GAP_OPEN];

		//the number of n_gap_ext corresponds to GAP_EXT in the left column,
		n_gap_ext += source_col[GAP_EXT];

		//gap_open is to be inserted into the rest of sequences :)
		n_gap_open = prof_size - n_gap_ext - n_gap_term_open - n_gap_term_ext;
#else
		n_gap_open = prof_size;
		n_gap_ext = 0;		
#endif
	}
}


// ****************************************************************************
void CProfile::DP_SolveGapsProblemWhenContinuing(size_t source_col_id, size_t prof_width, size_t prof_size, CProfile *profile, size_t &n_gap_ext, size_t &n_gap_term_ext)
{
	// inside a gap
	//if (source_col_id == prof_width) 
	//	n_gap_term_ext += prof_size;
	//else
	if (source_col_id == prof_width)
	{
		//continuation of a terminal gap on the right side of the profile 
		n_gap_term_ext = prof_size;
		n_gap_ext      = 0;
	}
	else	//continuation of a gap inside the profile 
	{

#ifndef NO_GAP_CORRECTION
		//the number of n_gap_term_ext corresponds to the number of GAP_TERM_OPEN in the column on the right,  
		n_gap_term_ext = profile->counters.get_value(source_col_id + 1, GAP_TERM_OPEN);

		counter_t *source_col = profile->counters.get_column(source_col_id);

		//the number of n_gap_term_ext corresponds to the number of GAP_TERM_OPEN in the column on the left, 
		n_gap_term_ext += source_col[GAP_TERM_OPEN];

		//the number of n_gap_term_ext corresponds to the number of GAP_TERM_EXT in the column on the left,
		n_gap_term_ext += source_col[GAP_TERM_EXT];

		n_gap_ext = prof_size - n_gap_term_ext;
#else
		n_gap_ext = prof_size;
		n_gap_term_ext = 0;			
#endif
	}
}

// ****************************************************************************
// Compute the range of cells that must be computed during the guided profile alignment according to the 
// mapping of column in the source alignment
void CProfile::FindRowRanges(vector<int> *column_mapping1, vector<int> *column_mapping2, vector<pair<int, int>> &row_ranges)
{
	int size = column_mapping1->size();
	int width = column_mapping2->size();
	int radius = params->guided_alignment_radius;

	row_ranges.resize(size+2);
	for(int i = 0; i <= size; ++i)
		row_ranges[i] = make_pair(column_mapping2->size()+1, 0);
	
	int i1 = 0;
	int i2 = 0;
	int i_res_max = max(column_mapping1->back(), column_mapping2->back());

	column_mapping1->push_back(i_res_max+1);
	column_mapping2->push_back(i_res_max+1);

	for(int i_res = 0; i_res <= i_res_max; ++i_res)
	{
		if((*column_mapping1)[i1] == i_res)
			i1++;
		if((*column_mapping2)[i2] == i_res)
			i2++;

		if(i2 - radius < 0)
			row_ranges[i1].first = 0;
		else
			row_ranges[i1].first = min(row_ranges[i1].first, i2 - radius);

		if(i2 + radius > width)
			row_ranges[i1].second = width;
		else
			row_ranges[i1].second = max(row_ranges[i1].second, i2 + radius);

		if(i1 - radius > 0)
			row_ranges[i1-radius].second = max(row_ranges[i1-radius].second, i2);
		if(i1 + radius <= size)
			row_ranges[i1+radius].first = min(row_ranges[i1+radius].first, i2);
	}

	for(int i = max(0, i1-radius); i <= size; ++i)
		row_ranges[i].second = width;

	for(int i = 0; i <= min(size, radius); ++i)
		row_ranges[i].first = 1;
}

// ****************************************************************************
void CProfile::AlignSeqProf(CProfile *profile1, CProfile *profile2, vector<int> *column_mapping1, vector<int> *column_mapping2)
{
	size_t prof1_width = profile1->width;
	size_t prof2_width = profile2->width;

	//size_t prof1_card = profile1->data.size();
	size_t prof1_card = 1;
	size_t prof2_card = profile2->data.size();
	
	// The profile1 contains a single sequence so there are no gaps here
	vector<symbol_t> &seq1 = const_cast<vector<symbol_t>&>(profile1->data[0]->symbols);

	score_t gap_open = params->gap_open;
	score_t gap_ext = params->gap_ext;
	score_t gap_term_open = params->gap_term_open;
	score_t gap_term_ext = params->gap_term_ext;

	CDPMatrix matrix(prof1_width + 1, prof2_width + 1);
	matrix.set_zeros();

	dp_row_t curr_row(prof2_width + 1);
	dp_row_t prev_row(prof2_width + 1);

	CProfileValues<score_t, NO_SYMBOLS> &scores1 = const_cast<CProfileValues<score_t, NO_SYMBOLS>&>(profile1->scores);
	CProfileValues<score_t, NO_SYMBOLS> &scores2 = const_cast<CProfileValues<score_t, NO_SYMBOLS>&>(profile2->scores);

	bool is_guided = column_mapping1 != nullptr && column_mapping2 != nullptr;
	vector<pair<int, int>> row_ranges;

	if(is_guided)
		FindRowRanges(column_mapping1, column_mapping2, row_ranges);
	else
		for(int i = 0; i <= (int) prof1_width; ++i)
			row_ranges.push_back(make_pair(0, prof2_width));

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
	}


	// Calculate matrix interior
	for (size_t i = 1; i <= prof1_width; ++i)
	{
		// Precompute scores for gaps for current and previous row of profile1
		score_t prof1_gap_open_prev = scores1.get_value(i - 1, GAP_OPEN);
		score_t prof1_gap_open_curr = scores1.get_value(i, GAP_OPEN);
		score_t prof1_gap_term_open_prev = scores1.get_value(i - 1, GAP_TERM_OPEN);
		score_t prof1_gap_term_open_curr = scores1.get_value(i, GAP_TERM_OPEN);
		score_t prof1_gap_ext_curr = scores1.get_value(i, GAP_EXT);
		score_t prof1_gap_term_ext_curr = scores1.get_value(i, GAP_TERM_EXT);

		// Boundary conditions
		curr_row[0].D = -infty;
		curr_row[0].H = -infty;
		matrix.set_dir_all(i, 0, direction_t::V);

		if(row_ranges[i].first)
			curr_row[row_ranges[i].first-1].D = curr_row[row_ranges[i].first-1].H = curr_row[row_ranges[i].first-1].V = -infty;

		if (i < prof1_width)
		{
			if (i == 1)
				curr_row[0].V = max(prev_row[0].D, prev_row[0].V) + prof1_gap_term_open_curr * prof2_card;
			else
				curr_row[0].V = max(prev_row[0].D, prev_row[0].V) + prof1_gap_term_ext_curr * prof2_card;

			for(int j = row_ranges[i].second+1; j <= min(row_ranges[i+1].second, (int) prof2_width); ++j)
				curr_row[j].D = curr_row[j].H = curr_row[j].V = -infty;
		}
		else
			curr_row[0].V = -infty;

		unsigned char *matrix_cell = matrix.get_cell(i, max(1, row_ranges[i].first) - 1);
		size_t max_j = min(row_ranges[i].second, (int) prof2_width);
		for (size_t j = max(1, row_ranges[i].first); j <= max_j; ++j)
		{
			// Get current cell of matrix for faster access
			++matrix_cell;
			//...........................................................................
			// Calcualte score for D
			//...........................................................................

			//analyzing an alignment of the column j of the profile 2 with the column i of the profile 1
			score_t *scores2_column = scores2.get_column(j);

			//an alignment of a column with another column, after an an alignment of two columns  
			score_t t = scores2_column[seq1[i]];
			score_t t_D = prev_row[j - 1].D + t;

			//..........................................................................
			//an alignment of a column with another column, after inserting a column of gaps before, into the profile 1
			score_t t_H = prev_row[j - 1].H;

			t_H = t_H + t;

			//..........................................................................
			//an alignment of a column with another column, after inserting a column of gaps before, into the profile 2 
			score_t t_V = prev_row[j - 1].V;

			//One has to check if there is a need of exchanging gap_opens/gap_term_opens (from the column of the profile 2 that is being just aligned) 
			//into gap_extensions/gap_term_extensions
			//If so, the value of the variable t has to be corrected

			t_V += t + gaps_prof2_change[j];

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

			score_t gap_corr = (i < prof1_width ? prof2_gaps[j].open : prof2_gaps[j].term_open);
			t_D = curr_row[j - 1].D + gap_corr;

			//the case when a sequence of gaps is continued (another columnn of gaps is inserted)
			t_H = curr_row[j - 1].H + (i < prof1_width ? prof2_gaps[j].ext : prof2_gaps[j].term_ext);

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

// ****************************************************************************
score_t CProfile::CalculateTotalScore(void)
{
	int i, j, len;

	score_t score = 0;

	// Calculate total score only if necessary!
	if(!params->enable_total_score_calculation)
	{
		total_score = 0;
		return total_score;
	}

	auto score_matrix = params->score_matrix;
	int64_t n_gap_open = 0;
	int64_t n_gap_term_open = 0;
	int64_t n_gap_ext = 0;
	int64_t n_gap_term_ext = 0;

	// Calculate score including only symbols (not gaps).
	// Additionally, estimates the number of gaps (all gaps are treated here as gap extensions)
	for(i = 1; i <= (int) width; ++i)
	{
		auto counters_col = counters.get_column(i);

		for(symbol_t c1 = 0; c1 < NO_AMINOACIDS; ++c1)
		{
			auto cnt_c1 = counters_col[c1];

			if (cnt_c1)
			{
				for (symbol_t c2 = c1 + 1; c2 < NO_AMINOACIDS; ++c2)
					if(counters_col[c2])
						score += score_matrix[c1][c2] * cnt_c1 * counters_col[c2];
				score += score_matrix[c1][c1] * cnt_c1 * (cnt_c1 - 1) / 2;
			}
		}

		int64_t n_gaps = counters_col[GAP_OPEN] + counters_col[GAP_EXT];
		int64_t n_gaps_term = counters_col[GAP_TERM_OPEN] + counters_col[GAP_TERM_EXT];
		int64_t n_symbols = data.size() - n_gaps - n_gaps_term;
		n_gap_ext      += n_symbols * n_gaps;
		n_gap_term_ext += n_symbols * n_gaps_term;
	}

	// Determine beginnings and ends of gaps
	size_t matrix_raw_size = ((width+1) + 1) * (width + 1) / 2;
	uint32_t *raw_gap_matrix = new uint32_t[matrix_raw_size];
	uint32_t **gap_matrix = new uint32_t*[width+1];
	uint64_t raw_offset = 0;

	for(i = 1; i <= (int) width; ++i)
	{
		gap_matrix[i] = &raw_gap_matrix[raw_offset];
		raw_offset += width + 1 - (i-1);
	}

	vector<pair<int, int>> gap_positions;
	vector<bool> gap_len_active(width+1, false);

	for(i = 0; i < (int) data.size(); ++i)
	{
		CGappedSequence *gs = data[i];
		int seq_pos = 1;
		int gap_len = 0;

		for(j = 0; j <= (int) gs->size; ++j)
		{
			gap_len = gs->n_gaps[j];
			if(gap_len)
			{
				if(!gap_len_active[gap_len])		// lazy init
				{
					gap_len_active[gap_len] = true;
					fill_n(gap_matrix[gap_len], width - gap_len + 2, 0);
				}

				if(!gap_matrix[gap_len][seq_pos]++)
					gap_positions.push_back(make_pair(gap_len, seq_pos));
				seq_pos += gap_len + 1;
			}
			else
				seq_pos++;
		}
	}

	// Calculate cumulated counters of gaps in ranges, i.e.
	// cell [i, s] contains the number of gaps included in range [i, i+s-1]
	uint32_t *raw_gap_ranges = new uint32_t[matrix_raw_size];
	uint32_t **gap_ranges = new uint32_t*[width+1];
	raw_offset = 0;

	for(i = 1; i <= (int) width; ++i)
	{
		gap_ranges[i] = &raw_gap_ranges[raw_offset];
		raw_offset += width + 1 - (i-1);
	}

	// Gaps of length 1
	if(gap_len_active[1])
		for (i = 1; i <= (int) width; ++i)
			gap_ranges[1][i] = gap_matrix[1][i];
	else
		for (i = 1; i <= (int) width; ++i)
			gap_ranges[1][i] = 0;

	// Gaps of length 2
	if(gap_len_active[2])
		for (i = 1; i <= (int) width - 1; ++i)
			gap_ranges[2][i] = gap_matrix[1][i] + gap_matrix[1][i + 1] + gap_matrix[2][i];
	else
		for (i = 1; i <= (int) width - 1; ++i)
			gap_ranges[2][i] = gap_matrix[1][i] + gap_matrix[1][i + 1];

	// Gaps longer than 2
	for (len = 3; len <= (int) width; ++len)
	{
/*		for (i = 1; i <= width - len + 1; ++i)
			gap_ranges[len][i] = gap_ranges[len - 1][i] + gap_ranges[len - 1][i + 1] - gap_ranges[len - 2][i + 1] + gap_matrix[len][i];*/
		uint32_t *dest = &gap_ranges[len][1];
		uint32_t *src1 = &gap_ranges[len-1][1];
		uint32_t *src2 = &gap_ranges[len-1][2];
		uint32_t *src3 = &gap_ranges[len-2][2];
		uint32_t *src4 = &gap_matrix[len][1];

		uint32_t *dest_end = &gap_ranges[len][width-len+1];

		if(gap_len_active[len])
			for(; dest <= dest_end; ++dest)
				*dest = *src1++ + *src2++ - *src3++ + *src4++;
		else
			for(; dest <= dest_end; ++dest)
				*dest = *src1++ + *src2++ - *src3++;
	}

	
	// Determine the exact number of gap_opens and gap_term_opens
	int64_t size = data.size();
	for (auto &x : gap_positions)
	{
		int len = x.first;
		int i = x.second;

		// No. of gaps containing the current gaps within
		int64_t n_gaps_inside = gap_ranges[width][1];
		if(len > 2)
			n_gaps_inside += gap_ranges[len - 2][i + 1];
		if(i + len - 2 > 0)
			n_gaps_inside -= gap_ranges[i + len - 2][1];
		if(i+1 <= (int) width)
			n_gaps_inside -= gap_ranges[width - i][i + 1];
		int cur_n_gaps = gap_matrix[len][i];
		n_gaps_inside -= cur_n_gaps;
		
		if(i == 1 || i+len-1 == width)
			n_gap_term_open += (size - cur_n_gaps - n_gaps_inside) * cur_n_gaps;
		else
			n_gap_open += (size - cur_n_gaps - n_gaps_inside) * cur_n_gaps;
	}
		
	n_gap_ext      -= n_gap_open;
	n_gap_term_ext -= n_gap_term_open;

	score += n_gap_ext * params->gap_ext + n_gap_open * params->gap_open +
		n_gap_term_ext * params->gap_term_ext + n_gap_term_open * params->gap_term_open;

	delete[] gap_matrix;
	delete[] raw_gap_matrix;

	delete[] gap_ranges;
	delete[] raw_gap_ranges;

	total_score = score;

	return score;
}

// ****************************************************************************
void CProfile::Swap(CProfile &profile)
{
	data.swap(profile.data);
	scores.swap(profile.scores);
	counters.swap(profile.counters);

	swap(width, profile.width);
	swap(total_score, profile.total_score);
}

