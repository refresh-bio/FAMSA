/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#include "../core/msa.h"
#include "../core/input_file.h"

#include <algorithm>
#include <set>
#include <random>
#include <map>
#include <list>
#include <stack>
#include <thread>

#include <iostream>
#include <iomanip>
#include <numeric>
#include <fstream>
#include <sstream>

#include "../core/output_file.h"
#include "../core/NewickTree.h"

#include "../libs/vectorclass.h"

using namespace std;

#define SHOW_PROGRESS

double CFAMSA::SM_MIQS[24][24] = {
//	   A     R     N     D     C     Q     E     G     H     I     L     K     M     F     P     S     T     W     Y     V     B     Z     X     *
	{  3.2, -1.3, -0.4, -0.4,  1.5, -0.2, -0.4,  0.4, -1.2, -1.3, -1.4, -0.7, -1.0, -2.3, -0.1,  0.8,  0.8, -3.6, -2.4,  0.0, -6.1, -6.1, -6.1, -6.1},	// A
	{ -1.3,  6.2, -0.1, -1.5, -2.7,  1.8, -0.7, -1.9,  0.9, -2.4, -2.5,  3.3, -1.1, -3.3, -1.1, -0.3, -0.9, -3.8, -1.9, -2.3, -6.1, -6.1, -6.1, -6.1},	// R
	{ -0.4, -0.1,  5.1,  2.6, -1.6,  0.9,  0.8,  0.2,  1.0, -3.6, -3.5,  0.7, -2.3, -3.5, -1.4,  0.9,  0.0, -4.5, -1.5, -2.6, -6.1, -6.1, -6.1, -6.1},	// N
	{ -0.4, -1.5,  2.6,  5.7, -3.7,  0.9,  2.7, -0.5,  0.3, -4.5, -4.6,  0.4, -3.3, -5.8, -0.3,  0.3, -0.2, -5.3, -3.9, -3.5, -6.1, -6.1, -6.1, -6.1},	// D
	{  1.5, -2.7, -1.6, -3.7, 11.7, -2.8, -3.2, -1.7, -1.2,  0.2, -2.3, -3.2,  0.1, -2.8, -2.8,  1.0,  0.0, -6.1, -0.7,  1.8, -6.1, -6.1, -6.1, -6.1},	// C
	{ -0.2,  1.8,  0.9,  0.9, -2.8,  3.6,  2.1, -1.6,  1.2, -2.2, -1.9,  1.7, -0.4, -2.4, -0.4,  0.4,  0.1, -5.4, -2.8, -1.8, -6.1, -6.1, -6.1, -6.1},	// Q
	{ -0.4, -0.7,  0.8,  2.7, -3.2,  2.1,  4.3, -1.3, -0.2, -3.3, -2.8,  1.1, -2.3, -4.1,  0.0,  0.4, -0.2, -5.8, -2.4, -2.3, -6.1, -6.1, -6.1, -6.1},	// E
	{  0.4, -1.9,  0.2, -0.5, -1.7, -1.6, -1.3,  7.6, -1.6, -5.4, -4.8, -1.7, -3.6, -4.6, -1.6,  0.0, -1.9, -4.8, -4.5, -3.8, -6.1, -6.1, -6.1, -6.1},	// G
	{ -1.2,  0.9,  1.0,  0.3, -1.2,  1.2, -0.2, -1.6,  7.5, -2.2, -1.9,  0.0, -2.1,  0.0, -1.5,  0.0, -0.2, -0.3,  2.1, -2.3, -6.1, -6.1, -6.1, -6.1},	// H
	{ -1.3, -2.4, -3.6, -4.5,  0.2, -2.2, -3.3, -5.4, -2.2,  4.6,  3.1, -2.3,  1.7,  0.7, -3.7, -2.8, -0.7, -0.7, -0.8,  3.3, -6.1, -6.1, -6.1, -6.1},	// I
	{ -1.4, -2.5, -3.5, -4.6, -2.3, -1.9, -2.8, -4.8, -1.9,  3.1,  4.6, -2.4,  3.2,  2.1, -2.8, -2.9, -1.6, -0.2,  0.0,  2.0, -6.1, -6.1, -6.1, -6.1},	// L
	{ -0.7,  3.3,  0.7,  0.4, -3.2,  1.7,  1.1, -1.7,  0.0, -2.3, -2.4,  3.6, -1.1, -3.7, -0.1,  0.0,  0.0, -4.0, -2.3, -2.0, -6.1, -6.1, -6.1, -6.1},	// K
	{ -1.0, -1.1, -2.3, -3.3,  0.1, -0.4, -2.3, -3.6, -2.1,  1.7,  3.2, -1.1,  5.4,  1.4, -2.8, -1.8, -0.8, -2.1, -0.9,  1.4, -6.1, -6.1, -6.1, -6.1},	// M
	{ -2.3, -3.3, -3.5, -5.8, -2.8, -2.4, -4.1, -4.6,  0.0,  0.7,  2.1, -3.7,  1.4,  7.4, -3.7, -2.6, -2.3,  4.2,  5.2, -0.3, -6.1, -6.1, -6.1, -6.1},	// F
	{ -0.1, -1.1, -1.4, -0.3, -2.8, -0.4,  0.0, -1.6, -1.5, -3.7, -2.8, -0.1, -2.8, -3.7,  8.4, -0.1, -0.5, -3.6, -4.5, -2.5, -6.1, -6.1, -6.1, -6.1},	// P
	{  0.8, -0.3,  0.9,  0.3,  1.0,  0.4,  0.4,  0.0,  0.0, -2.8, -2.9,  0.0, -1.8, -2.6, -0.1,  3.1,  1.6, -3.5, -1.5, -1.4, -6.1, -6.1, -6.1, -6.1},	// S
	{  0.8, -0.9,  0.0, -0.2,  0.0,  0.1, -0.2, -1.9, -0.2, -0.7, -1.6,  0.0, -0.8, -2.3, -0.5,  1.6,  3.8, -5.3, -2.1, -0.1, -6.1, -6.1, -6.1, -6.1},	// T
	{ -3.6, -3.8, -4.5, -5.3, -6.1, -5.4, -5.8, -4.8, -0.3, -0.7, -0.2, -4.0, -2.1,  4.2, -3.6, -3.5, -5.3, 14.8,  4.9, -3.3, -6.1, -6.1, -6.1, -6.1},	// W
	{ -2.4, -1.9, -1.5, -3.9, -0.7, -2.8, -2.4, -4.5,  2.1, -0.8,  0.0, -2.3, -0.9,  5.2, -4.5, -1.5, -2.1,  4.9,  8.3, -1.2, -6.1, -6.1, -6.1, -6.1},	// Y
	{  0.0, -2.3, -2.6, -3.5,  1.8, -1.8, -2.3, -3.8, -2.3,  3.3,  2.0, -2.0,  1.4, -0.3, -2.5, -1.4, -0.1, -3.3, -1.2,  3.5, -6.1, -6.1, -6.1, -6.1},	// V
	{ -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1},	// B
	{ -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1},	// Z
	{ -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1},	// X
	{ -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1}};	// *


// *******************************************************************
CFAMSA::CFAMSA()
{
	DetermineInstructionSet();

	n_threads = 1;
	final_profile = nullptr;

	init_sm();
}

// *******************************************************************
CFAMSA::~CFAMSA()
{
	if(final_profile)
		delete final_profile;
}

// *******************************************************************
void CFAMSA::init_sm()
{
	score_matrix.resize(NO_AMINOACIDS);

#ifdef HUGE_ALIGNEMENTS
	for(int i = 0; i < NO_AMINOACIDS; ++i)
	{
		score_vector.push_back(SM_MIQS[i][i]);
		for(int j = 0; j < NO_AMINOACIDS; ++j)
			score_matrix[i].push_back(SM_MIQS[i][j]);
	}
#else
	for(int i = 0; i < NO_AMINOACIDS; ++i)
	{
		score_vector.push_back((score_t) round(SM_MIQS[i][i] * cost_cast_factor));
		for(int j = 0; j < NO_AMINOACIDS; ++j)
			score_matrix[i].push_back((score_t) round(SM_MIQS[i][j] * cost_cast_factor));
	}
#endif
}

// *******************************************************************
bool CFAMSA::SetSequences(vector<CSequence> &_sequences)
{
	timers[3].StartTimer();

	sequences = _sequences;

	std::stable_sort(sequences.begin(), sequences.end(), [](const CSequence& a, const CSequence& b)->bool {
		return a.length > b.length || (a.length == b.length && a.data < b.data);
	});
	
	timers[3].StopTimer();

	return true;
}


// *******************************************************************
bool CFAMSA::SetSequences(vector<CSequence> &&_sequences)
{
	timers[3].StartTimer();

	sequences = std::move(_sequences);

	std::stable_sort(sequences.begin(), sequences.end(), [](const CSequence& a, const CSequence& b)->bool {
				return a.length > b.length || (a.length == b.length && a.data < b.data);
	});

	timers[3].StopTimer();

	return true;
}


// *******************************************************************
bool CFAMSA::SetParams(CParams &_params)
{
	if(sequences.empty())
		return false;

	params = _params;

	if(params.enable_gap_rescaling)
	{
		double gap_scaler = log2(sequences.size() / (double) params.scaler_log);
		if(sequences.size() < params.scaler_log)
			gap_scaler = 1.0;
		else
			gap_scaler = 1.0 + (gap_scaler / params.scaler_div);

		params.gap_ext		 = (score_t) (params.gap_ext       * gap_scaler);
		params.gap_open	 	 = (score_t) (params.gap_open      * gap_scaler);
		params.gap_term_ext  = (score_t) (params.gap_term_ext  * gap_scaler);
		params.gap_term_open = (score_t) (params.gap_term_open * gap_scaler);
	}

	params.score_matrix = score_matrix;
	params.score_vector = score_vector;

	n_threads = params.n_threads;
	// adjust automatically
	if (n_threads == 0) {
		n_threads = std::thread::hardware_concurrency();
		// if hardware_concurrency fails
		if (n_threads == 0) {
			n_threads = 8;
		}
	}

	return true;
}

// *******************************************************************
void CFAMSA::DetermineInstructionSet()
{
	int x = instrset_detect();

	if(x >= 0 && x <= 8)
		instruction_set = (instruction_set_t) x;
	else if(x < 0)
		instruction_set = instruction_set_t::none;
	else
		instruction_set = instruction_set_t::avx2;
}

#ifdef DEVELOPER_MODE
// *******************************************************************
// Compute LCS length for two sequences in the classical way - just for development
double CFAMSA::GetLCS(CSequence &seq1, CSequence &seq2)
{
	int **dp_row = new int*[2];

	for(int i = 0; i < 2; ++i)
		dp_row[i] = new int[seq2.length+1];

	fill(dp_row[0], dp_row[0]+seq2.length+1, 0);

	for(int i = 1; i <= (int) seq1.length; ++i)
	{
		int ii = i % 2;
		dp_row[ii][0] = 0;
		for(int j = 1; j <= (int) seq2.length; ++j)
			if(seq1.data[i-1] == seq2.data[j-1])
				dp_row[ii][j] = dp_row[!ii][j-1]+1;
			else
				dp_row[ii][j] = max(dp_row[ii][j-1], dp_row[!ii][j]);
	}

	return dp_row[seq1.length % 2][seq2.length];
}
#endif

// *******************************************************************
// Compute Guide Tree
bool CFAMSA::ComputeGuideTree()
{
	guide_tree.clear();
	gt_stats.clear();

	// Insert leaves (sequences)
	for(size_t i = 0; i < sequences.size(); ++i)
	{
		guide_tree.push_back(make_pair(-1, -1));
		prof_cardinalities.push_back(1);

		gt_stats.push_back(make_pair(1, 0));
	}

	// Construct guide tree using selected algorithm
	if (params.guide_tree == GT_method::single_linkage)
		SingleLinkage();
	else if (params.guide_tree == GT_method::UPGMA)
		UPGMA();
#ifdef DEVELOPER_MODE
	else if (params.guide_tree == GT_method::chained)
		GuideTreeChained();
#endif
	else if (params.guide_tree == GT_method::imported)
		if (!ImportGuideTreeFromNewick())			// file can be unexisting, so need to terminate the code here
			return false;

	// store guide tree in Newick format
	if (params.guide_tree_out_file.length() > 0) {
		ExportGuideTreeToNewick();
	}

	// Convert sequences into gapped sequences
	gapped_sequences.reserve(sequences.size());

	for(auto &p: sequences)
		gapped_sequences.push_back(CGappedSequence(p));

	return true;
}

// *******************************************************************
// Compute Alignment according to guide tree
bool CFAMSA::ComputeAlignment()
{
	final_profile = new CProfile(&params);

	final_profile->Clear();

	CProfileQueue pq(&gapped_sequences, &profiles, &guide_tree);

	params.sackin_index = pq.GetSackinIndex();

	vector<thread *> workers(n_threads, nullptr);

	uint32_t computed_prof = 0;
	mutex mtx;

	size_t ref_thr = params.thr_internal_refinement;

	for(uint32_t i = 0; i < n_threads; ++i)
		workers[i] = new thread([&]{
			CGappedSequence *gs;
			CProfile *prof1;
			CProfile *prof2;
			CProfile *prof_sol;
			size_t prof_id;

			while(pq.GetTask(prof_id, gs, prof1, prof2))
			{
				if(gs != nullptr)
				{
					prof_sol = new CProfile(*gs, &params);
				}
				else
				{
					// Internal refinement
					if (prof1->Size() + prof2->Size() > ref_thr)
					{
						if (prof1->Size() <= ref_thr && prof1->Size() > 2)
							RefineAlignment(prof1);
						if (prof2->Size() <= ref_thr && prof2->Size() > 2)
							RefineAlignment(prof2);
					}

					prof_sol = new CProfile(prof1, prof2, &params);

					delete prof1;
					delete prof2;
				}

				pq.AddSolution(prof_id, prof_sol);

				if (params.very_verbose_mode)
				{
					lock_guard<mutex> lck(mtx);
					computed_prof++;

					if (computed_prof % 100 == 0 || 
						(computed_prof % 10 == 0 && (double) computed_prof / (2 * sequences.size() - 1) > 0.95))
					{
						cout << "Computing alignment - " << fixed << setprecision(1) << 100.0 * computed_prof / (2 * sequences.size()-1) << 
							"\%    (" << computed_prof << " of " << (2 * sequences.size() - 1) << ")\r";
						fflush(stdout);
					}
				}
			}
		});

	for(auto p : workers)
	{
		p->join();
		delete p;
	}
	workers.clear();

	final_profile = profiles.begin()->second;

	return true;
}

// *******************************************************************
// Single linkage guide tree construction (SLINK algorithm) - uses native threads rather than OpenMP
void CFAMSA::SingleLinkage()
{
	int next;
	int n_seq = sequences.size();

	int prefetch_offset = 64*2;
	int prefetch_offset_2nd = 128*2;

	double indel_exp = params.indel_exp;

	vector<int> pi(n_seq + max(prefetch_offset, prefetch_offset_2nd), 0);
	vector<double> lambda(n_seq);
	vector<double> *sim_vector;

	// Temporarily resize the sequences by adding unknown symbols (just to simplify the pairwise comparison)
	size_t max_seq_len = max_element(sequences.begin(), sequences.end(), [](const CSequence &x, CSequence &y){return x.length < y.length; })->length;

	for (int i = 1; i < n_seq; ++i)
		sequences[i].data.resize(max_seq_len, UNKNOWN_SYMBOL);

	CSingleLinkageQueue slq(&sequences, n_seq, n_threads * 8);
	vector<thread *> workers(n_threads, nullptr);

	mutex mtx;

	// Calculation of similarities is made in working threads
	for (uint32_t i = 0; i < n_threads; ++i)
		workers[i] = new thread([&]{
			CLCSBP lcsbp(instruction_set);
			int row_id;
			vector<CSequence> *sequences;
			vector<double> *sim_vector;
			uint32_t lcs_lens[4];

			while (slq.GetTask(row_id, sequences, sim_vector))
			{
				for (int j = 0; j < row_id / 4; ++j)
				{
					lcsbp.GetLCSBP(&(*sequences)[row_id], &(*sequences)[j * 4 + 0], &(*sequences)[j * 4 + 1], &(*sequences)[j * 4 + 2], &(*sequences)[j * 4 + 3],
						lcs_lens[0], lcs_lens[1], lcs_lens[2], lcs_lens[3]);
					for (int k = 0; k < 4; ++k)
					{
						double indel = (*sequences)[row_id].length + (*sequences)[j * 4 + k].length - 2 * lcs_lens[k];
						//double seq_lens = (*sequences)[row_id].length + (*sequences)[j * 4 + k].length;
						(*sim_vector)[j * 4 + k] = lcs_lens[k] / pow(indel, indel_exp);
					}
				}

				if (row_id / 4 * 4 < row_id)
				{
					lcsbp.GetLCSBP(&(*sequences)[row_id],
						(row_id / 4 * 4 + 0 < row_id) ? &(*sequences)[row_id / 4 * 4 + 0] : nullptr,
						(row_id / 4 * 4 + 1 < row_id) ? &(*sequences)[row_id / 4 * 4 + 1] : nullptr,
						(row_id / 4 * 4 + 2 < row_id) ? &(*sequences)[row_id / 4 * 4 + 2] : nullptr,
						(row_id / 4 * 4 + 3 < row_id) ? &(*sequences)[row_id / 4 * 4 + 3] : nullptr,
						lcs_lens[0], lcs_lens[1], lcs_lens[2], lcs_lens[3]);
					for (int k = 0; k < 4 && row_id / 4 * 4 + k < row_id; ++k)
					{
						double indel = (*sequences)[row_id].length + (*sequences)[row_id / 4 * 4 + k].length - 2 * lcs_lens[k];
						//double seq_lens = (*sequences)[row_id].length + (*sequences)[row_id / 4 * 4 + k].length;
						(*sim_vector)[row_id / 4 * 4 + k] = lcs_lens[k] / pow(indel, indel_exp);
					}
				}

				slq.RegisterSolution(row_id);
			}
		});

	// Single linkage algorithm is here
	for (int i = 0; i < n_seq; ++i)
	{
		pi[i] = i;
		lambda[i] = -infty_double;

		if (params.very_verbose_mode && i % (100) == 0)
		{
			cout << "Computing guide tree - " << fixed << setprecision(1) << 100.0 * ((double)i * (i + 1) / 2) / ((double)n_seq * (n_seq + 1) / 2) << "\%    (" << i << " of " << n_seq << ")  \r";
			fflush(stdout);
		}

		slq.GetSolution(i, sim_vector);

		auto p_lambda = lambda.begin();
		auto p_sim_vector = (*sim_vector).begin();
		auto p_pi = pi.begin();

		for (int j = 0; j < i; ++j)
		{
			next = pi[j];

#ifdef _MSC_VER					// Visual C++
			_mm_prefetch((const char*)&(*sim_vector)[*(p_pi+prefetch_offset)], 2);
#endif
#ifdef __GNUC__
//			__builtin_prefetch((&(*sim_vector)[pi[j + prefetch_offset]]), 1, 2);
			__builtin_prefetch((&(*sim_vector)[*(p_pi+prefetch_offset)]), 1, 2);
#endif

			auto &x = (*sim_vector)[next];

			if (isgreater(*p_lambda, *p_sim_vector))
			{
				x = max(x, *p_sim_vector);
			}
			else
			{
				x = max(x, *p_lambda);
				*p_pi = i;
				*p_lambda = *p_sim_vector;
			}

			++p_pi;
			++p_lambda;
			++p_sim_vector;
		}
		
		slq.ReleaseSolution(i);

		p_pi = pi.begin();
		p_lambda = lambda.begin();
		for (int j = 0; j < i; ++j)
		{
#ifdef _MSC_VER					// Visual C++
			_mm_prefetch((const char*)&lambda[*(p_pi + prefetch_offset_2nd)], 0);
#endif
#ifdef __GNUC__
			__builtin_prefetch((&lambda[*(p_pi + prefetch_offset_2nd)]), 1, 0);
#endif

			next = *p_pi;
			if(isgreaterequal(lambda[next], *p_lambda))
				*p_pi = i;

			++p_pi;
			++p_lambda;
		}
	}

	for (auto p : workers)
	{
		p->join();
		delete p;
	}
	workers.clear();

	if (params.very_verbose_mode)
	{
		cout << "Computing guide tree - 100.0\%                                        \r";
		fflush(stdout);
	}

	vector<int> elements(n_seq - 1);
	for (int i = 0; i < n_seq - 1; ++i)
		elements[i] = i;

#ifdef DEBUG_MODE
	identity /= n_seq * (n_seq - 1) / 2.0;
#endif

	stable_sort(elements.begin(), elements.end(), [&](int x, int y){
		return lambda[x] > lambda[y];
	});

	vector<int> index(n_seq);
	for (int i = 0; i < n_seq; ++i)
		index[i] = i;

	for (int i = 0; i < n_seq - 1; ++i)
	{
		int j = elements[i];
		next = pi[j];
		guide_tree.push_back(make_pair(index[j], index[next]));
		index[next] = n_seq + i;
	}

	// Bring the sequences to the valid length
	for (int i = 1; i < n_seq; ++i)
		sequences[i].data.resize(sequences[i].length, UNKNOWN_SYMBOL);
}


// *******************************************************************
#ifdef DEVELOPER_MODE
void CFAMSA::GuideTreeChained()
{
	mt19937 rnd;

	if (sequences.size() < 2)
		return;

	vector<int> idx(sequences.size());

	for (int i = 0; i < sequences.size(); ++i)
		idx[i] = i;

	random_device rd;
	
	// Skip some number of initial values
	for (int i = 0; i < params.guide_tree_seed; ++i)
		rd();

	mt19937 g(rd());

	shuffle(idx.begin(), idx.end(), g);

	guide_tree.push_back(make_pair(idx[0], idx[1]));

	for (int i = 2; i < sequences.size(); ++i)
		guide_tree.push_back(make_pair(idx[i], guide_tree.size() - 1));
}
#endif

// *******************************************************************
bool CFAMSA::ImportGuideTreeFromNewick()
{
	// Load newick description
	ifstream newickFile;
	newickFile.open(params.guide_tree_in_file);
	if (!newickFile.good()) {
		return false;
	}

	std::stringstream ss;
	ss << newickFile.rdbuf();
	std::string description(ss.str());
	auto newend = std::remove_if(description.begin(), description.end(),
		[](char c)->bool { return c == '\r' || c == '\n';  });
	description.erase(newend, description.end());

	// Load guide tree
	NewickTree nw(params.verbose_mode);
	nw.parse(sequences, description, guide_tree);

	return true;
}

bool CFAMSA::ExportGuideTreeToNewick()
{
	// store guide tree
	string description;
	NewickTree nw(params.verbose_mode);
	nw.store(sequences, guide_tree, description);

	// Open file
	ofstream newickFile;
	newickFile.open(params.guide_tree_out_file);
	if (!newickFile.good()) {
		return false;
	}

	newickFile << description;

	return true;
}

// *******************************************************************
bool CFAMSA::ExportDistanceMatrix(float*matrix, size_t size, const std::string& fname) {
	ofstream file(fname);
	if (!file) {
		return false;
	}

	for (size_t i = 0; i < size; ++i) {
		for (size_t j = 0; j < i; ++j) {
			const size_t id = UPGMA_TriangleSubscript(i, j);
			float d = matrix[id];
			file << d << ", ";
		}
		file << std::endl;
	}
}


// *******************************************************************
void CFAMSA::RefineRandom(CProfile* profile_to_refine, vector<size_t> &dest_prof_id)
{
	for(size_t i = 0; i < profile_to_refine->data.size(); ++i)
		dest_prof_id.push_back(rnd_rfn() % 2);

	if(count(dest_prof_id.begin(), dest_prof_id.end(), 0) == 0 ||
		count(dest_prof_id.begin(), dest_prof_id.end(), 1) == 0)		// Both profiles must contain at least 1 sequence
	{
		size_t id = rnd_rfn() % dest_prof_id.size();
		dest_prof_id[id] = !dest_prof_id[id];
	}
}

// *******************************************************************
void CFAMSA::RefineMostEmptyAndFullColumn(CProfile *profile_to_refine, vector<size_t> &dest_prof_id, vector<size_t> &gap_stats, bool valid_gap_stats)
{
	size_t size = profile_to_refine->data.front()->gapped_size;
	size_t card = profile_to_refine->data.size();

	dest_prof_id.clear();

	if(!valid_gap_stats)
		profile_to_refine->GetGapStats(gap_stats);

	vector<pair<size_t, size_t>> tmp;

	for(size_t i = 1; i <= size; ++i)
	{
		int x = min(gap_stats[i], card - gap_stats[i]);
		if(x > 0)
			tmp.push_back(make_pair(i, x));
	}

	sort(tmp.begin(), tmp.end(), [](const pair<size_t, size_t> &x, const pair<size_t, size_t> &y){
		if(x.second != y.second)
			return x.second < y.second;
		else
			return x.first < y.first;
	});

	if(tmp.empty())
	{
		RefineRandom(profile_to_refine, dest_prof_id);
		return;
	}

	size_t col_id = tmp[rnd_rfn() % tmp.size()].first;

	int first_prof_id = 0;
	int second_prof_id = 1;

	if(profile_to_refine->data[0]->GetSymbol(col_id) == GAP)
		swap(first_prof_id, second_prof_id);

	for(size_t j = 0; j < card; ++j)
		if(profile_to_refine->data[j]->GetSymbol(col_id) == GAP)
			dest_prof_id.push_back(first_prof_id);
		else
			dest_prof_id.push_back(second_prof_id);
}

// *******************************************************************
// Refine alignment
#ifdef DEBUG_MODE
bool CFAMSA::RefineAlignment(string output_file_name)
#else
bool CFAMSA::RefineAlignment(CProfile *&profile_to_refine)
#endif
{
	// Restart generator
	rnd_rfn.seed(5489u);

	if (params.enable_auto_refinement && profile_to_refine->Size() > params.thr_refinement)
		return true;

	size_t n_ref = params.n_refinements;
	size_t n_seq = profile_to_refine->Size();

	vector<size_t> gap_stats;

	if(n_ref > 2*n_seq)
		n_ref = 2*n_seq;
	if(n_ref > 0 && n_ref < 100 && n_seq < 100)
		n_ref = 100;

#ifdef DEBUG_MODE
	FILE *f_stat;

	if(output_file_name != "")
	{
		vector<CGappedSequence*> result;
		GetAlignment(result);
			
		COutputFile out_file;
			
		out_file.PutSequences(result);
		out_file.SaveFile(output_file_name + to_string(0));

		f_stat = fopen((output_file_name + "_stats").c_str(), "wt");

		fprintf(f_stat, "%d  %f\n", final_profile->width, final_profile->CalculateTotalScore());
	}
#endif

	int n_ref_succ = 0;
	score_t prev_total_score = profile_to_refine->CalculateTotalScore();

	sort(profile_to_refine->data.begin(), profile_to_refine->data.end(), [](CGappedSequence *p, CGappedSequence *q){return p->id < q->id; });

	vector<size_t> dest_prof_id;
	vector<vector<size_t>> old_dest_prof_ids;

	vector<int> column_mapping1, column_mapping2;

	size_t i_ref;
	size_t i_succ_ref;
	bool valid_gap_stats = false;
#ifdef DEBUG_MODE
	int ref_upd[2] = {0};
	int hist_size[20] = {0};
#endif

	for(i_ref = i_succ_ref = 0; i_succ_ref < n_ref && i_ref < 20*n_ref; ++i_ref)
	{
		if (params.very_verbose_mode)
		{
			cout << "Computing refinement - " << fixed << setprecision(1) << 100.0 * (double) i_succ_ref / (double) n_ref << "\%    (" << i_succ_ref << " of " << n_ref << ")  \r";
			fflush(stdout);
		}

		CProfile profile1(&params), profile2(&params);

		RefineMostEmptyAndFullColumn(profile_to_refine, dest_prof_id, gap_stats, valid_gap_stats);
		valid_gap_stats = true;

		if(find(old_dest_prof_ids.begin(), old_dest_prof_ids.end(), dest_prof_id) == old_dest_prof_ids.end())
		{
			// Split into two profiles
			for(size_t i = 0; i < profile_to_refine->data.size(); ++i)
				if(dest_prof_id[i])
					profile1.AppendRawSequence(*profile_to_refine->data[i]);
				else
					profile2.AppendRawSequence(*profile_to_refine->data[i]);

			// Condense the profiles (remove empty columns)
			profile1.Condense(column_mapping1);
			profile2.Condense(column_mapping2);

			profile1.OptimizeGaps();
			profile2.OptimizeGaps();

			profile1.Size();
			profile2.Size();

#ifdef DEBUG_MODE
			int size_min = min(p1_size, p2_size);
			hist_size[min(9, size_min)]++;
#endif

			CProfile* prof = new CProfile(&params);

			prof->Align(&profile1, &profile2, &column_mapping1, &column_mapping2);
			sort(prof->data.begin(), prof->data.end(), [](CGappedSequence *p, CGappedSequence *q){return p->id < q->id; });

			if (*prof != *profile_to_refine)		// if the new profile is the same as previous do not score it
			{
				prof->CalculateTotalScore();
#ifdef DEBUG_MODE
				ref_upd[0]++;
#endif

				if (prof->total_score >= prev_total_score)
					{
					prev_total_score = prof->total_score;
					swap(profile_to_refine, prof);
					++n_ref_succ;
					old_dest_prof_ids.clear();
					valid_gap_stats = false;
#ifdef DEBUG_MODE
					ref_upd[1]++;
					hist_size[10+min(9, size_min)]++;
#endif
				}
			}

			delete prof;

			old_dest_prof_ids.push_back(dest_prof_id);
			i_succ_ref++;

#ifdef DEBUG_MODE
			if(output_file_name != "")
			{
				vector<CGappedSequence*> result;
				GetAlignment(result);

				COutputFile out_file;

				out_file.PutSequences(result);
				out_file.SaveFile(output_file_name + to_string(i_ref+1));
				fprintf(f_stat, "%d  %f  p.sizes: %5d %5d\n", final_profile->width, final_profile->total_score, p1_size, p2_size);
			}
#endif
		}
	}

#ifdef DEBUG_MODE
	if(output_file_name != "")
		fclose(f_stat);
#endif

	return true;
}

// *******************************************************************
bool CFAMSA::GetAlignment(vector<CGappedSequence*> &result)
{
	result = final_profile->data;

	sort(result.begin(), result.end(), [](CGappedSequence *p, CGappedSequence *q){return p->id < q->id;});

	return !result.empty();
}

#ifdef DEVELOPER_MODE
// *******************************************************************
bool CFAMSA::LoadRefSequences()
{
	// ***** Read input file
	CInputFile ref_file;

	if (params.verbose_mode)
		cerr << "Processing reference sequence set: " << params.ref_file_name << "\n";

	if (!ref_file.ReadFile(params.ref_file_name))
	{
		cout << "Error: no (or incorrect) input file\n";
		return false;
	}

	ref_file.StealSequences(ref_sequences);

	return true;
}

// *******************************************************************
int CFAMSA::SubTreeSize(set<int> &seq_ids)
{
	vector<pair<int, set<int>>> node_stats;
	set<int> tmp_set;
	int n_seq = sequences.size();

	// Calculate stats for single sequence nodes
	for (int i = 0; i < n_seq; ++i)
	{
		tmp_set.clear();
		if (seq_ids.count(i) > 0)
			tmp_set.insert(i);

		node_stats.push_back(make_pair(1, tmp_set));
	}

	// Calculate stats for internal nodes
	for (int i = n_seq; i < 2 * n_seq - 1; ++i)
	{
		tmp_set.clear();
		auto &x = node_stats[guide_tree[i].first];
		auto &y = node_stats[guide_tree[i].second];
		tmp_set.insert(x.second.begin(), x.second.end());
		tmp_set.insert(y.second.begin(), y.second.end());

		node_stats.push_back(make_pair(x.first + y.first, tmp_set));

		if (tmp_set.size() == ref_sequences.size())
			return node_stats.back().first;
	}

	return 0;
}

// *******************************************************************
uint64_t CFAMSA::CalculateRefSequencesSubTreeSize(double *monte_carlo_subtree_size)
{
	set<int> ref_seq_ids;
	int n_seq = sequences.size();
	int r = 0;

	if (ref_sequences.size() == 1)
		return 1;

	// Find the ids of the referential sequences in the input file
	for (int i = 0; i < n_seq; ++i)
	{
		bool is_ref = false;
		for (auto &y : ref_sequences)
			if (sequences[i].id == y.id)
				is_ref = true;

		if (is_ref)
			ref_seq_ids.insert(i);
	}

	r = SubTreeSize(ref_seq_ids);

	if (monte_carlo_subtree_size)
	{
		mt19937 mt;
		double mc_r = 0;

		for (int i = 0; i < monte_carlo_trials; ++i)
		{
			set<int> mc_seq_ids;

			while (mc_seq_ids.size() < ref_seq_ids.size())
				mc_seq_ids.insert(mt() % n_seq);
			mc_r += SubTreeSize(mc_seq_ids);
		}

		*monte_carlo_subtree_size = mc_r / (double) monte_carlo_trials;
	}

	return r;
}
#endif

// *******************************************************************
bool CFAMSA::ComputeMSA()
{
	if (params.verbose_mode)
	{
		cerr << "Params:\n";
		cerr << "  gap open cost (rescalled): " << params.gap_open << "\n";
		cerr << "  gap extension cost (rescalled): " << params.gap_ext << "\n";
		cerr << "  gap terminal open cost (rescalled): " << params.gap_term_open << "\n";
		cerr << "  gap terminal extension cost (rescalled): " << params.gap_term_ext << "\n";
		cerr << "  no. of refinements: " << params.n_refinements << "\n";
		cerr << "  gap cost scaller log-term: " << params.scaler_log << "\n";
		cerr << "  gap cost scaller div-term: " << params.scaler_div << "\n";
		cerr << "  refinement threshold: " << params.thr_refinement << "\n";
		cerr << "  enable gap rescaling: " << params.enable_gap_rescaling << "\n";
		cerr << "  enable gap optimization: " << params.enable_gap_optimization << "\n";
		cerr << "  enable total score calculation: " << params.enable_total_score_calculation << "\n";
		cerr << "  enable auto refinement: " << params.enable_auto_refinement << "\n";
		cerr << "  guided alignment radius: " << params.guided_alignment_radius << "\n";
		cerr << "  no. of threads: " << params.n_threads << "\n";
		cerr << "  guide tree method: ";
		switch (params.guide_tree)
		{
		case GT_method::chained:
			cerr << "chained\n";
			break;
		case GT_method::imported:
			cerr << "imported\n";
			cerr << "  guide tree file: " << params.guide_tree_in_file << "\n";
			break;
		case GT_method::single_linkage:
			cerr << "single linkage\n";
			break;
		case GT_method::UPGMA:
			cerr << "UPGMA\n";
			break;
		}
	}
		
#ifdef DEVELOPER_MODE
	if (params.test_ref_sequences)
		LoadRefSequences();
#endif

	if (params.very_verbose_mode)
	{
		string instr_names[] = { "None", "SSE", "SSE2", "SSE3", "SSE3S", "SSE41", "SSE42", "AVX", "AVX2" };
		cout << "Instruction set determined: " << instr_names[(int)instruction_set] << "\n";
	}

	timers[0].StartTimer();
	if (params.very_verbose_mode)
	{
		cout << "Computing guide tree...\r";
		fflush(stdout);
	}
	if(!ComputeGuideTree())
		return false;
	timers[0].StopTimer();
	if (params.very_verbose_mode)
	{
		cout << "Computing guide tree - complete                                          \n";
		fflush(stdout);
	}

#ifdef DEVELOPER_MODE
	double monte_carlo_subtree_size;
	if(params.test_ref_sequences)
		params.ref_seq_subtree_size = CalculateRefSequencesSubTreeSize(&monte_carlo_subtree_size);
#endif

	timers[1].StartTimer();
	if (params.very_verbose_mode)
	{
		cout << "Computing alignment...\r";
		fflush(stdout);
	}
	if (!ComputeAlignment())
		return false;
	if (params.very_verbose_mode)
	{
		cout << "Computing alignment - complete                                           \n";
		fflush(stdout);
	}
	timers[1].StopTimer();

	timers[2].StartTimer();
	if (params.very_verbose_mode)
	{
		cout << "Computing refinement...\r";
		fflush(stdout);
	}
	if (!RefineAlignment(final_profile))
		return false;
	if (params.very_verbose_mode)
	{
		cout << "Computing refinement - complete                                         \n";
		fflush(stdout);
	}

	timers[2].StopTimer();

	if (params.verbose_mode)
	{
		cerr << "Sequence sorting                                : " << timers[3].GetElapsedTime() << "s\n";
		cerr << "Guide tree construction (incl. similatiry calc.): " << timers[0].GetElapsedTime() << "s\n";
		cerr << "Alignment construction                          : " << timers[1].GetElapsedTime() << "s\n";
		cerr << "Iterative refinement                            : " << timers[2].GetElapsedTime() << "s\n";

		cerr << "No. of sequences : " << sequences.size() << "\n";
		cerr << "Sackin index for guide tree: " << params.sackin_index << "\n";
		cerr << "Sackin index for guide tree (normalized): " << params.sackin_index / (double)sequences.size() << "\n";

#ifdef DEVELOPER_MODE
		if (params.test_ref_sequences)
		{
			cerr << "Ref. seq. subtree size: " << params.ref_seq_subtree_size << "\n";
			cerr << "Monte Carlo subtree size: " << monte_carlo_subtree_size << "\n";
		}

#ifdef LOG_STATS
		FILE *out = fopen("execution.stats", "wt");
		fprintf(out, "[stats]\n");
		fprintf(out, "No_sequences=%d\n", sequences.size());
		fprintf(out, "Sackin_idx=%lld\n", params.sackin_index);
		fprintf(out, "Sackin_idx_norm=%.3f\n", params.sackin_index / (double) sequences.size());
		fprintf(out, "Ref_seq_subtree_size=%lld\n", params.ref_seq_subtree_size);
		fprintf(out, "Monte_carlo_subtree_size=%.1f", monte_carlo_subtree_size);
		fclose(out);
#endif
#endif
	}

	return true;
}
