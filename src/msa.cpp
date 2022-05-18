/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#include "msa.h"
#include "./tree/GuideTree.h"
#include "./tree/SingleLinkage.h"
#include "./tree/FastTree.h"
#include "./tree/MSTPrim.h"
#include "./tree/UPGMA.h"
#include "./tree/NeighborJoining.h"
#include "./tree/DistanceCalculator.h"
#include "./core/io_service.h"
#include "./utils/log.h"

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

using namespace std;

#define SHOW_PROGRESS

double CFAMSA::SM_MIQS[24][24] = {
//	   A     R     N     D     C     Q     E     G     H     I     L     K     M     F     P     S     T     W     Y     V     B     Z     X     *
	{3.2, -1.3, -0.4, -0.4,  1.5, -0.2, -0.4,  0.4, -1.2, -1.3, -1.4, -0.7, -1.0, -2.3, -0.1,  0.8,  0.8, -3.6, -2.4,  0.0, -6.1, -6.1, -6.1, -6.1},	// A
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
CFAMSA::CFAMSA(CParams& _params) 
	: params(_params), instruction_set(params.instruction_set), final_profile(nullptr)
{
	initScoreMatrix();
}

// *******************************************************************
CFAMSA::~CFAMSA()
{
	delete final_profile;
}

// *******************************************************************
void CFAMSA::initScoreMatrix()
{
	score_matrix.resize(NO_AMINOACIDS);

#ifdef HUGE_ALIGNEMENTS
	for(int i = 0; i < NO_AMINOACIDS; ++i)
	{
		score_vector.emplace_back(SM_MIQS[i][i]);
		for(int j = 0; j < NO_AMINOACIDS; ++j)
			score_matrix[i].emplace_back(SM_MIQS[i][j]);
	}
#else
	for(int i = 0; i < NO_AMINOACIDS; ++i)
	{
		score_vector.emplace_back((score_t) round(SM_MIQS[i][i] * cost_cast_factor));
		for(int j = 0; j < NO_AMINOACIDS; ++j)
			score_matrix[i].emplace_back((score_t) round(SM_MIQS[i][j] * cost_cast_factor));
	}
#endif
}

// *******************************************************************
void CFAMSA::adjustParams(int n_seqs)
{

	if ((params.gt_heuristic != GT::None) && (n_seqs < params.heuristic_threshold)) {
		params.gt_heuristic = GT::None;
	}

	if(params.enable_gap_rescaling)
	{
		double gap_scaler = log2(n_seqs / (double) params.scaler_log);
		if(n_seqs < (int) params.scaler_log)
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
}

// *******************************************************************
std::shared_ptr<AbstractTreeGenerator> CFAMSA::createTreeGenerator(const CParams& params) {
	
	std::shared_ptr<AbstractTreeGenerator> gen = nullptr;
	
	// distance only mode
	if (params.export_distances) {
		LOG_VERBOSE << "Calculating distances and storing in: " << params.output_file_name;

		std::shared_ptr<AbstractTreeGenerator> calculator;

		if (params.distance == Distance::indel_div_lcs) {
			gen = make_shared<DistanceCalculator<Distance::indel_div_lcs>>(
				params.n_threads, params.instruction_set, 
				params.output_file_name, params.generate_square_matrix, params.calculate_pid);

		}
		else if (params.distance == Distance::sqrt_indel_div_lcs) {
			gen = make_shared<DistanceCalculator<Distance::sqrt_indel_div_lcs>>(
				params.n_threads, params.instruction_set, 
				params.output_file_name, params.generate_square_matrix, params.calculate_pid);
		}
	
		return gen;
	}

	if (params.gt_method == GT::SLINK || (params.gt_heuristic != GT::None && params.gt_method == GT::MST_Prim)) {
		// single linkage (use also when MST used as a partial generator)
		if (params.distance == Distance::indel_div_lcs) {
			gen = make_shared<SingleLinkage<Distance::indel_div_lcs>>(params.n_threads, params.instruction_set);
		}
		else if (params.distance == Distance::sqrt_indel_div_lcs) {
			gen = make_shared<SingleLinkage<Distance::sqrt_indel_div_lcs>>(params.n_threads, params.instruction_set);
		}
	}
	else if (params.gt_method == GT::MST_Prim) {
		// Prim's minimum spanning tree 
		if (params.distance == Distance::indel_div_lcs) {
			gen = make_shared<MSTPrim<Distance::indel_div_lcs>>(params.n_threads, params.instruction_set);
		}
		else if (params.distance == Distance::sqrt_indel_div_lcs) {
			gen = make_shared<MSTPrim<Distance::sqrt_indel_div_lcs>>(params.n_threads, params.instruction_set);
		}
	}
	else if (params.gt_method == GT::UPGMA || params.gt_method == GT::UPGMA_modified) {
		// UPGMA
		if (params.distance == Distance::indel_div_lcs) {
			gen = make_shared<UPGMA<Distance::indel_div_lcs>>(params.n_threads, params.instruction_set,
				(params.gt_method == GT::UPGMA_modified));
		}
		else if (params.distance == Distance::sqrt_indel_div_lcs) {
			gen = make_shared<UPGMA<Distance::sqrt_indel_div_lcs>>(params.n_threads, params.instruction_set,
				(params.gt_method == GT::UPGMA_modified));
		}
	}
	else if (params.gt_method == GT::NJ) {
		// neighbour joining
		if (params.distance == Distance::indel_div_lcs) {
			gen = make_shared<NeighborJoining<Distance::indel_div_lcs>>(params.n_threads, params.instruction_set);
		}
		else if (params.distance == Distance::sqrt_indel_div_lcs) {
			gen = make_shared<NeighborJoining<Distance::sqrt_indel_div_lcs>>(params.n_threads, params.instruction_set);
		}
	} else {

		throw std::runtime_error("Error: Illegal guide tree method.");
	}

	// check if heuristic was specified
	if (params.gt_heuristic != GT::None) {

		// part-tree versus medoid-tree
		shared_ptr<IClustering> clustering = (params.gt_heuristic == GT::PartTree)
			? nullptr
			: make_shared<CLARANS>(params.cluster_fraction, params.cluster_iters);

		// verify distance measure
		if (params.distance == Distance::indel_div_lcs) {
			gen = make_shared<FastTree<Distance::indel_div_lcs>>(
				params.n_threads,
				params.instruction_set,
				dynamic_pointer_cast<IPartialGenerator>(gen),
				params.subtree_size,
				clustering,
				params.sample_size);
		}
		else if (params.distance == Distance::sqrt_indel_div_lcs) {
			gen = make_shared<FastTree<Distance::sqrt_indel_div_lcs>>(
				params.n_threads,
				params.instruction_set,
				dynamic_pointer_cast<IPartialGenerator>(gen),
				params.subtree_size,
				clustering,
				params.sample_size);
		}
	}

	return gen;
}

// *******************************************************************
// Compute Alignment according to guide tree
bool CFAMSA::ComputeAlignment(std::vector<std::pair<int,int>>& guide_tree)
{
	final_profile = new CProfile(&params);

	final_profile->Clear();

	CProfileQueue pq(&gapped_sequences, &profiles, &guide_tree, params.n_threads);

	vector<thread *> workers(params.n_threads, nullptr);

	uint32_t computed_prof = 0;
	mutex mtx;

	size_t ref_thr = params.thr_internal_refinement;

	for(uint32_t i = 0; i < params.n_threads; ++i)
		workers[i] = new thread([&]{
			CGappedSequence *gs;
			CProfile *prof1;
			CProfile *prof2;
			CProfile *prof_sol;
			size_t prof_id;
			//bool only_task;
			uint32_t no_threads;
			uint32_t no_rows_per_box;

//			while(pq.GetTask(prof_id, gs, prof1, prof2, only_task))
			while(pq.GetTask(prof_id, gs, prof1, prof2, no_threads, no_rows_per_box))
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

					prof_sol = new CProfile(prof1, prof2, &params, no_threads, no_rows_per_box);

					delete prof1;
					delete prof2;
				}

				pq.AddSolution(prof_id, prof_sol);

				if (params.very_verbose_mode)
				{
					lock_guard<mutex> lck(mtx);
					computed_prof++;

					if (computed_prof % 100 == 0 || 
						(computed_prof % 10 == 0 && (double) computed_prof / (2 * gapped_sequences.size() - 1) > 0.95))
					{
						LOG_NORMAL << "Computing alignment - " << fixed << setprecision(1) << 100.0 * computed_prof / (2 * gapped_sequences.size()-1) << 
							"%    (" << computed_prof << " of " << (2 * gapped_sequences.size() - 1) << ")\r";
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
void CFAMSA::RefineRandom(CProfile* profile_to_refine, vector<size_t> &dest_prof_id)
{
	for(size_t i = 0; i < profile_to_refine->data.size(); ++i)
		dest_prof_id.emplace_back(rnd_rfn() % 2);

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
		int x = (int) min(gap_stats[i], card - gap_stats[i]);
		if(x > 0)
			tmp.emplace_back(i, x);
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
			dest_prof_id.emplace_back(first_prof_id);
		else
			dest_prof_id.emplace_back(second_prof_id);
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
		LOG_DEBUG << "Computing refinement - " << fixed << setprecision(1) << 100.0 * (double) i_succ_ref / (double) n_ref << "%    (" << i_succ_ref << " of " << n_ref << ")  \r";
			
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

			// TODO: Enable parallelization here!
			prof->Align(&profile1, &profile2, 1, 0, &column_mapping1, &column_mapping2);
			sort(prof->data.begin(), prof->data.end(), [](CGappedSequence *p, CGappedSequence *q){return p->id < q->id; });

			if (!(*prof == *profile_to_refine))		// if the new profile is the same as previous do not score it
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

			old_dest_prof_ids.emplace_back(dest_prof_id);
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
	if (!final_profile) {
		return false;
	}		

	result = final_profile->data;

	sort(result.begin(), result.end(), [](CGappedSequence *p, CGappedSequence *q){return p->id < q->id;});

	return !result.empty();
}

#ifdef DEVELOPER_MODE
// *******************************************************************
bool CFAMSA::LoadRefSequences()
{
	// ***** Read input file	
	LOG_DEBUG << "Processing reference sequence set: " << params.ref_file_name << "\n";

	if (IOService::loadFasta(params.ref_file_name, ref_sequences) < 2)
	{
		LOG_NORMAL << "Error: no (or incorrect) input file\n";
		return false;
	}

	return true;
}
#endif

// *******************************************************************
bool CFAMSA::ComputeMSA(vector<CSequence>& sequences)
{
	adjustParams((int)sequences.size());
	
	LOG_VERBOSE 
		<< "Params:\n"
		<< "  no. of threads: " << params.n_threads << "\n"
		<< "  guide tree method: " << GT::toString(params.gt_method) << "\n"
		<< "  guide tree speeding up heuristic: " << GT::toString(params.gt_heuristic) << "\n";

	if (params.gt_method == GT::imported) {
		LOG_VERBOSE << "  guide tree file: " << params.guide_tree_in_file << "\n";
	}
		
	LOG_VERBOSE
		<< "  distance measure: " << dist2str(params.distance) << "\n\n"

		<< "Advanced params:\n"
		<< "  no. of refinements: " << params.n_refinements << "\n"
		<< "  refinement threshold: " << params.thr_refinement << "\n"
		<< "  gap open cost (rescalled): " << params.gap_open << "\n"
		<< "  gap extension cost (rescalled): " << params.gap_ext << "\n"
		<< "  gap terminal open cost (rescalled): " << params.gap_term_open << "\n"
		<< "  gap terminal extension cost (rescalled): " << params.gap_term_ext << "\n"
		<< "  gap cost scaller log-term: " << params.scaler_log << "\n"
		<< "  gap cost scaller div-term: " << params.scaler_div << "\n"
		<< "  enable gap rescaling: " << params.enable_gap_rescaling << "\n"
		<< "  enable gap optimization: " << params.enable_gap_optimization << "\n"
		<< "  enable total score calculation: " << params.enable_total_score_calculation << "\n"
		<< "  enable auto refinement: " << params.enable_auto_refinement << "\n"
		<< "  guided alignment radius: " << params.guided_alignment_radius << "\n\n";
	

#ifdef DEVELOPER_MODE
	if (params.test_ref_sequences)
		LoadRefSequences();
#endif

	string instr_names[] = { "None", "SSE", "SSE2", "SSE3", "SSE3S", "SSE41", "SSE42", "AVX", "AVX2" };
	LOG_VERBOSE << "Hardware configuration: " << endl
		<< " Number of threads: " << params.n_threads << endl
		<< " Instruction set: " << instr_names[(int)instruction_set] << endl << endl;


	GuideTree tree;

	bool goOn = true;

	// store distance matrix
	if (params.export_distances) {
		LOG_VERBOSE << "Calculating distances and storing in: " << params.output_file_name;
		
		std::shared_ptr<AbstractTreeGenerator> calculator = createTreeGenerator(params);
		tree_structure tree;

		// update sequence identifiers
		for (int i = 0; i < (int)sequences.size(); ++i) {
			sequences[i].sequence_no = i;
		}

		(*calculator)(sequences, tree);
		LOG_VERBOSE << " [OK]" << endl;
		goOn = false; // break processing at this point
	}


	if (goOn) {
		timers[TIMER_SORTING].StartTimer();
		if (params.shuffle == -1) {
			LOG_VERBOSE << "Sorting sequences...";
			std::stable_sort(sequences.begin(), sequences.end(), [](const CSequence& a, const CSequence& b)->bool {
				return a.length > b.length || 
					(a.length == b.length && std::lexicographical_compare(a.data, a.data + a.data_size, b.data, b.data + b.data_size));
			});
			LOG_VERBOSE << " [OK]" << endl;
		}
		else {
			LOG_VERBOSE << "Shuffling sequences...";
			std::mt19937 mt(params.shuffle);
			std::shuffle(sequences.begin(), sequences.end(), mt);
			LOG_VERBOSE << " [OK]" << endl;
		}

		// update sequence identifiers
		for (int i = 0; i < (int)sequences.size(); ++i) {
			sequences[i].sequence_no = i;
		}
		timers[TIMER_SORTING].StopTimer();


		timers[TIMER_TREE_BUILD].StartTimer();
		if (params.gt_method == GT::imported) {
			LOG_VERBOSE << "Importing guide tree from: " << params.guide_tree_in_file;
			tree.loadNewick(params.guide_tree_in_file, sequences);
		}
		else {
			std::shared_ptr<AbstractTreeGenerator> gen = createTreeGenerator(params);
			LOG_VERBOSE << "Computing guide tree...";
			(*gen)(sequences, tree.raw());
		}
		LOG_VERBOSE << " [OK]" << endl;
		timers[TIMER_TREE_BUILD].StopTimer();
	}

	// store guide tree in Newick format 
	if (params.export_tree) {
		timers[TIMER_TREE_STORE].StartTimer();
		LOG_VERBOSE << "Storing guide tree in: " << params.output_file_name;
		tree.saveNewick(params.output_file_name, sequences);
		LOG_VERBOSE << " [OK]" << endl;
		timers[TIMER_TREE_STORE].StopTimer();
		goOn = false; // break processing at this point
	}


#ifdef DEVELOPER_MODE
	double monte_carlo_subtree_size;
	if (params.test_ref_sequences)
		params.ref_seq_subtree_size = tree.refSequencesSubTreeSize(
			sequences,
			ref_sequences,
			&monte_carlo_subtree_size);
#endif
	int64_t sackin;
	if (params.verbose_mode || params.very_verbose_mode) {
		sackin = tree.calculateSackinIndex();
		statistics.put("guide_tree.sackin", sackin);
		statistics.put("guide_tree.sackin_norm", sackin / (double)sequences.size());
	}

	if (goOn) {
		// Convert sequences into gapped sequences
		gapped_sequences.reserve(sequences.size());
		for (auto &p : sequences)
			gapped_sequences.emplace_back(std::move(p));
		std::vector<CSequence>().swap(sequences); // clear input vector

		timers[TIMER_ALIGNMENT].StartTimer();
		LOG_VERBOSE << "Computing alignment...";
		if (!ComputeAlignment(tree.raw()))
			return false;
		LOG_VERBOSE << "[OK]" << endl;
		timers[TIMER_ALIGNMENT].StopTimer();

		timers[TIMER_REFINMENT].StartTimer();
		LOG_VERBOSE << "Computing refinement...";
		if (!RefineAlignment(final_profile))
			return false;
		LOG_VERBOSE << "[OK]" << endl;
		timers[TIMER_REFINMENT].StopTimer();

		if (final_profile->Size() != gapped_sequences.size()) {
			throw std::runtime_error("Error: incomplete guide tree - report a bug");
		}
	}

	if (params.verbose_mode || params.very_verbose_mode) {
		statistics.put("time.sort", timers[TIMER_SORTING].GetElapsedTime());
		statistics.put("time.tree_build", timers[TIMER_TREE_BUILD].GetElapsedTime());
		statistics.put("time.tree_store", timers[TIMER_TREE_STORE].GetElapsedTime());
		statistics.put("time.alignment", timers[TIMER_ALIGNMENT].GetElapsedTime());
		statistics.put("time.refinement", timers[TIMER_REFINMENT].GetElapsedTime());
	}

#ifdef DEVELOPER_MODE
	if (params.test_ref_sequences)
	{
		LOG_VERBOSE
			<< "Ref. seq. subtree size: " << params.ref_seq_subtree_size << "\n"
			<< "Monte Carlo subtree size: " << monte_carlo_subtree_size << "\n";
	}
#endif

	return true;
}
