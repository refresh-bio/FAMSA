/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#include "params.h"
#include "../utils/log.h"

#if SIMD==SIMD_AVX1 || SIMD==SIMD_AVX2 || SIMD==SIMD_AVX512
#include "../utils/cpuid.h"
#endif

#include <thread>

//****************************************************************************
//
CParams::CParams() {
#ifdef HUGE_ALIGNMENTS
	gap_open = -gap_open_base;
	gap_ext = -gap_ext_base;
	gap_term_open = -gap_term_open_base;
	gap_term_ext = -gap_term_ext_base;
#else
	gap_open = (score_t)round(-cost_cast_factor * gap_open_base);
	gap_ext = (score_t)round(-cost_cast_factor * gap_ext_base);
	gap_term_open = (score_t)round(-cost_cast_factor * gap_term_open_base);
	gap_term_ext = (score_t)round(-cost_cast_factor * gap_term_ext_base);
#endif

	// verify instruction sets
#if SIMD==SIMD_AVX1 || SIMD==SIMD_AVX2 || SIMD==SIMD_AVX512
	if ((CPUID(1).ECX() >> 28) & 1)
		instruction_set = instruction_set_t::avx;
	if ((CPUID(7).EBX() >> 5) & 1)
		instruction_set = instruction_set_t::avx2;
#endif

#if SIMD==SIMD_NEON
	// !!! Check if NEON is supported
#endif

};


//****************************************************************************
// Show command-line parameters
void CParams::show_usage(bool expert)
{
	string bool2str[]{ "disabled", "enabled" };
	
	
	LOG_NORMAL
		<< "Usage:\n"
		<< "  famsa [options] <input_file> [<input_file_2>] <output_file>\n\n"

		<< "Positional parameters:\n"
		<< "  input_file, input_file_2 - input files in FASTA format; action depends on the number of input files:\n"
		<< "      * one input - multiple sequence alignment (input gaps, if present, are removed prior the alignment),\n"
		<< "      * two inputs - profile-profile aligment (input gaps are preserved).\n"
		<< "      First input can be replaced with STDIN string to read from standard input.\n"
		<< "  output_file - output file (pass STDOUT when writing to standard output); available outputs:\n"
		<< "      * alignment in FASTA format,\n"
		<< "      * guide tree in Newick format (-gt_export option specified),\n"
		<< "      * distance matrix in CSV format (-dist_export option specified),\n\n"

		<< "Options:\n"
		<< "  -help - print this message\n"
		<< "  -t <value> - no. of threads, 0 means all available (default: " << n_threads << ")\n"
		<< "  -v - verbose mode, show timing information (default: disabled)\n\n"

		<< "  -gt <sl | upgma | nj | import <file>> - guide tree method (default: sl):\n"
		<< "      * sl - single linkage \n"
		<< "      * upgma - UPGMA\n"
		<< "      * nj - neighbour joining\n"
		<< "      * import <file> - imported from a Newick file\n"
	//	<< "  -dist <measure> - pairwise distance measure:\n"
	//	<< "      * indel_div_lcs (default)\n"
	//	<< "      * sqrt_indel_div_lcs\n\n"

		<< "  -medoidtree - use MedoidTree heuristic for speeding up tree construction (default: disabled)\n"
		<< "  -medoid_threshold <n_seqs> - if specified, medoid trees are used only for sets with <n_seqs> or more\n"
		<< "  -gt_export - export a guide tree to output file in Newick format\n"
		<< "  -dist_export - export a distance matrix to output file in CSV format\n"
		<< "  -square_matrix - generate a square distance matrix instead of a default triangle\n"
		<< "  -pid - generate pairwise identity (the number of matching residues divided by the shorter sequence length) instead of distance\n"
		<< "  -keep-duplicates - keep duplicated sequences during alignment\n"
		<< "                     (default: disabled - duplicates are removed prior and restored after the alignment).\n\n"

		<< "  -gz - enable gzipped output (default: " << bool2str[gzippd_output] << ")\n"
		<< "  -gz-lev <value> - gzip compression level [0-9] (default: " << gzip_level << ")\n"
		<< "  -refine_mode <on | off | auto> - refinement mode (default: auto - the refinement is enabled for sets <= " << thr_refinement << " seq.)\n\n";

		
	if (expert) {
		LOG_NORMAL << "Advanced options:\n"
			<< "  -r <value> - no. of refinement iterations (default: " << n_refinements << ")\n"
			<< "  -go <value> - gap open cost (default: " << gap_open << ")\n"
			<< "  -ge <value> - gap extension cost (default: " << gap_ext << ")\n"
			<< "  -tgo <value> - terminal gap open cost (default: " << gap_term_open << ")\n"
			<< "  -tge <value> - terminal gap extenstion cost (default: " << gap_term_ext << ")\n"

			<< "  -gsd <value> - gap cost scaller div-term (default: " << scaler_div << ")\n"
			<< "  -gsl <value> - gap cost scaller log-term (default: " << scaler_log << ")\n"
			<< "  -dgr - disable gap cost rescaling (default: enabled)\n"
			<< "  -dgo - disable gap optimization (default: enabled)\n"
			<< "  -dsp - disable sum of pairs optimization during refinement (default: enabled)\n";
#ifdef DEVELOPER_MODE
		LOG_NORMAL << "  -ref <file_name> - load referential sequences (for benchmarks) and calculate the minimal subtree size containing them\n"
			<< "  -vv - very verbose mode, show timing information (default: disabled)\n";
#endif
		LOG_NORMAL << endl;
	}
	
}

// ****************************************************************************
// **** Parses parameters
bool CParams::parse(int argc, char** argv, bool& showExpert)
{
	std::vector<std::string> params(argv + 1, argv + argc);
	
	if (findSwitch(params, "-help")) {
		showExpert = true;
		return false;
	}

	if (params.size() < 2) {
		return false;
	}

	// advanced params
	findOption(params, "-go", gap_open_base);
	findOption(params, "-ge", gap_ext_base);
	findOption(params, "-tgo", gap_term_open_base);
	findOption(params, "-tge", gap_term_ext_base);
	findOption(params, "-gsd", scaler_div);
	findOption(params, "-gsl", scaler_log);
		
	enable_gap_rescaling = !findSwitch(params, "-dgr");
	enable_gap_optimization = !findSwitch(params, "-dgo");
	enable_total_score_calculation = !findSwitch(params, "-dsp");

	findOption(params, "-r", n_refinements);
	findOption(params, "-rt", thr_refinement);
	findOption(params, "-ri", thr_internal_refinement);
	
	string aux;
	if (findOption(params, "-refine_mode", aux)) {
		refinement_mode = Refinement::fromString(aux);
	}
	 
	// user params
	findOption(params, "-t", n_threads);

	auto next = params.end();
	if (findOption(params, "-gt", aux, next)) {
		gt_method = GT::fromString(aux);
		if (gt_method == GT::imported && next != params.end()) {
			guide_tree_in_file = *next;
			params.erase(next);

		}
		else if (gt_method == GT::chained && next != params.end()) {
			guide_tree_seed = stoi(*next);
			params.erase(next);
		}
	}

	if (findOption(params, "-dist", aux)) {
		distance = str2dist(aux);
	}

	if (findSwitch(params, "-parttree")) {
		gt_heuristic = GT::PartTree;
	}

	if (findSwitch(params, "-medoidtree")) {
		gt_heuristic = GT::ClusterTree;
	}

	findOption(params, "-medoid_threshold", heuristic_threshold);
	findOption(params, "-subtree_size", subtree_size);
	findOption(params, "-sample_size", sample_size);
	findOption(params, "-cluster_fraction", cluster_fraction);
	findOption(params, "-cluster_iters", cluster_iters);

	export_tree = findSwitch(params, "-gt_export");
	export_distances = findSwitch(params, "-dist_export");
	generate_square_matrix = findSwitch(params, "-square_matrix");
	calculate_pid = findSwitch(params, "-pid");
	gzippd_output = findSwitch(params, "-gz");

	int g_lev = gzip_level;
	if (findOption(params, "-gz-lev", g_lev) && (g_lev < 0 || g_lev > 12))
	{
		LOG_NORMAL << "Incorrect gzip level: " << g_lev << " was changed to default value: " << gzip_level << endl;
		g_lev = gzip_level;
	}
	gzip_level = g_lev;

	keepDuplicates = findSwitch(params, "-keep-duplicates");
		
#ifdef DEVELOPER_MODE
	findOption(params, "-shuffle", shuffle)
		
	else if (cur_par == "-ref") {
		test_ref_sequenes = true;
		ref_file_name = argv[argno++];
	}
#endif
		
	verbose_mode = findSwitch(params, "-v");
	very_verbose_mode = findSwitch(params, "-vv");
	
	// remaining parameters are positionals
	if (params.size() < 2 || params.size() > 3) {
		return false;
	}
	else if (params.size() == 2) {
		input_file_name = params[0];
		output_file_name = params[1];
	}
	else if (params.size() == 3) {
		// profile-profile aligment
		input_file_name = params[0];
		input_file_name_2 = params[1];
		output_file_name = params[2];
		profile_aligning = true;
	}

#ifdef HUGE_ALIGNMENTS
	gap_open = -gap_open_base;
	gap_ext = -gap_ext_base;
	gap_term_open = -gap_term_open_base;
	gap_term_ext = -gap_term_ext_base;
#else
	gap_open = (score_t)round(-cost_cast_factor * gap_open_base);
	gap_ext = (score_t)round(-cost_cast_factor * gap_ext_base);
	gap_term_open = (score_t)round(-cost_cast_factor * gap_term_open_base);
	gap_term_ext = (score_t)round(-cost_cast_factor * gap_term_ext_base);
#endif

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
