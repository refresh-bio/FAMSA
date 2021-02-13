/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#ifndef _PARAMS_H
#define _PARAMS_H

#include <vector>
#include <math.h>
#include <string>
#include "defs.h"
#include "TreeDefs.h"

using namespace std;

struct CParams
{
	score_t gap_open;
	score_t gap_ext;
	score_t gap_term_open;
	score_t gap_term_ext;
	uint32_t n_refinements;
	uint32_t scaler_log;
	uint32_t scaler_div;
	uint32_t thr_refinement;
	uint32_t thr_internal_refinement;
	bool enable_gap_rescaling;
	bool enable_gap_optimization;
	bool enable_total_score_calculation;
	bool enable_auto_refinement;
	bool verbose_mode;
	bool very_verbose_mode;
	uint64_t ref_seq_subtree_size;
	double indel_exp;

	GT::Method gt_method;
	GT::Heuristic gt_heuristic;
	int heuristic_threshold;
	int guide_tree_seed;
	int subtree_size;
	int sample_size;
	float cluster_fraction;
	int cluster_iters;


	string guide_tree_in_file;
	bool export_distances;
	bool export_tree;

	bool test_ref_sequences;
	string ref_file_name;
	string output_file_name;

	int guided_alignment_radius;
	uint32_t n_threads;

	int64_t shuffle;

	vector<vector<score_t>> score_matrix;
	vector<score_t> score_vector;

	CParams(double _gap_open = -13.683, double _gap_ext = -1.246, double _gap_term_open = -0.619, double _gap_term_ext = -0.618, uint32_t _n_refinements = 100,
		uint32_t _scaler_log = 49, uint32_t _scaler_div = 18, uint32_t _thr_refinement = 1000,
		bool _enable_gap_rescaling = true, bool _enable_gap_optimization = true, bool _enable_total_score_calculation = true, bool _enable_auto_refinement = true,
		bool _verbose_mode = false, bool _very_verbose_mode = false) :
		guided_alignment_radius(50), n_threads(0)
	{
#ifdef HUGE_ALIGNMENTS
		gap_open      = _gap_open;
		gap_ext       = _gap_ext;
		gap_term_open = _gap_term_open;
		gap_term_ext  = _gap_term_ext;
#else
		gap_open      = (score_t) round(cost_cast_factor * _gap_open);
		gap_ext       = (score_t) round(cost_cast_factor * _gap_ext);
		gap_term_open = (score_t) round(cost_cast_factor * _gap_term_open);
		gap_term_ext  = (score_t) round(cost_cast_factor * _gap_term_ext);
#endif

		n_refinements  = _n_refinements;
		scaler_log     = _scaler_log;
		scaler_div     = _scaler_div;
		thr_refinement = _thr_refinement;
		thr_internal_refinement = 0;

		enable_gap_rescaling		   = _enable_gap_rescaling;
		enable_gap_optimization		   = _enable_gap_optimization;
		enable_total_score_calculation = _enable_total_score_calculation;
		enable_auto_refinement		   = _enable_auto_refinement;
		verbose_mode				   = _verbose_mode;
		very_verbose_mode			   = _very_verbose_mode;

		gt_method = GT::SLINK;
		gt_heuristic = GT::None;
		heuristic_threshold = 0;
		export_tree = false;
		export_distances = false;
		subtree_size = 50;
		sample_size = 500;

		cluster_fraction = 0.2;
		cluster_iters = 2;

		ref_seq_subtree_size = 0;
		test_ref_sequences = false;

		indel_exp = 1.0;
		shuffle = -1;
	};
};

#endif
