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
#include <sstream>
#include <algorithm>
#include "../core/defs.h"
#include "../tree/TreeDefs.h"


using namespace std;

class Refinement {
public:
	enum Mode {ON, OFF, AUTO};
	
	static std::string toString(Mode v) {
		switch (v) {
		case ON:				return "on";
		case OFF:				return "off";
		case AUTO:				return "auto";
		default:
			throw new std::runtime_error("Error: Illegal refinment mode.");
		}

		return "Unknown";
	}

	static Mode fromString(const std::string& name) {
		if (name == "on") { return ON; }
		if (name == "off") { return OFF; }
		if (name == "auto") { return AUTO; }

		// something went wrong
		throw new std::runtime_error("Error: Illegal refinment mode.");

		return ON;
	}
};


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
	
	Refinement::Mode refinement_mode		= Refinement::AUTO;
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
	bool keepDuplicates					= false;
	
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

	bool profile_aligning = false;
	string input_file_name;
	string input_file_name_2;
	string output_file_name;

	vector<vector<score_t>> score_matrix;
	vector<score_t> score_vector;

	CParams();
	bool parse(int argc, char** argv, bool& showExpert);
	void show_usage(bool expert);

protected:
	bool findSwitch(std::vector<std::string>& params, const std::string& name) {
		auto it = find(params.begin(), params.end(), name); // verbose mode
		if (it != params.end()) {
			params.erase(it);
			return true;
		}

		return false;
	}

	template <typename T>
	bool findOption(std::vector<std::string>& params, const std::string& name, T& v) {
		auto prevToEnd = std::prev(params.end());
		auto it = find(params.begin(), prevToEnd, name); // verbose mode
		if (it != prevToEnd) {
			std::istringstream iss(*std::next(it));
			if (iss >> v) {
				params.erase(it, it + 2);
				return true;
			}
		}

		return false;
	}

	template <typename T>
	bool findOption(
		std::vector<std::string>& params, 
		const std::string& name, 
		T& v, 
		std::vector<std::string>::iterator & next) {
		
		auto prevToEnd = std::prev(params.end());
		auto it = find(params.begin(), prevToEnd, name); // verbose mode
		if (it != prevToEnd) {
			std::istringstream iss(*std::next(it));
			if (iss >> v) {
				next = params.erase(it, it + 2);
				return true;
			}
		}

		return false;
	}
};

#endif
