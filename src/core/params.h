/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#ifndef _PARAMS_H
#define _PARAMS_H

#include <vector>
#include <math.h>
#include <string>
#include "../core/defs.h"
#include "../tree/TreeDefs.h"


using namespace std;

class CParams
{
private:
	double gap_open_base					= 14.85;
	double gap_ext_base						= 1.25;
	double gap_term_open_base				= 0.66;
	double gap_term_ext_base				= 0.66;
	
public:	
	score_t gap_open;
	score_t gap_ext;
	score_t gap_term_open;
	score_t gap_term_ext;

	uint32_t scaler_div = 7;
	uint32_t scaler_log = 45;
	int guided_alignment_radius = 50;

	bool enable_gap_rescaling				= true;
	bool enable_gap_optimization			= true;
	bool enable_total_score_calculation		= true;
	bool enable_auto_refinement				= true;
	
	uint32_t n_refinements					= 100;
	uint32_t thr_refinement					= 1000;
	uint32_t thr_internal_refinement		= 0;
	
	GT::Method gt_method			= GT::MST_Prim;
	GT::Heuristic gt_heuristic		= GT::None;
	Distance distance				= Distance::indel_div_lcs;
	int heuristic_threshold			= 0;
	
	int guide_tree_seed				= 0;
	int subtree_size				= 100;
	int sample_size					= 2000;
	float cluster_fraction			= 0.1f;
	int cluster_iters				= 2;

	string guide_tree_in_file;
	bool export_distances				= false;
	bool export_tree					= false;
	bool generate_square_matrix			= false;
	bool calculate_pid					= false;
	
	bool test_ref_sequences				 = false;
	uint64_t ref_seq_subtree_size = 0;
	string ref_file_name;
	
	int64_t shuffle = -1;
	uint32_t n_threads					= 0;
	
	bool gzippd_output					= false;
	int gzip_level						= 7;
	
	instruction_set_t instruction_set	= instruction_set_t::none;

	bool verbose_mode = false;
	bool very_verbose_mode = false;

	string input_file_name;
	string output_file_name;

	vector<vector<score_t>> score_matrix;
	vector<score_t> score_vector;

	CParams();
	bool parse(int argc, char** argv, bool& showExpert);
	void show_usage(bool expert);
};

#endif
