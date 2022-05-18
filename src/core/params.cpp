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
		<< "  famsa [options] <input_file> <output_file>\n\n"

		<< "Positional parameters:\n"
		<< "  input_file - input file in FASTA format (pass STDIN when reading from standard input)\n"
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
		<< "  -pid - generate percent identity instead of distance\n\n"

		<< "  -gz - enable gzipped output (default: " << bool2str[gzippd_output] << ")\n"
		<< "  -gz-lev <value> - gzip compression level [0-9] (default: " << gzip_level << ")\n\n";

	
	if (expert) {
		LOG_NORMAL << "Advanced options:\n"
			<< "  -r <value> - no. of refinement iterations (default: " << n_refinements << ")\n"
			<< "  -fr - force refinement (by default the refinement is disabled for sets larger than " << thr_refinement << " seq.)\n"
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
	// ugly workaround
	if (argc == 2 && string(argv[1]) == "-help") {
		showExpert = true;
		return false;
	}

	if (argc < 3)
		return false;

	int argno;

	for (argno = 1; argno + 2 < argc;)
	{
		string cur_par = argv[argno++];

		if (cur_par == "-help") {
			showExpert = true;
			return false;
		}
		else if (cur_par == "-go")
			gap_open_base = atof(argv[argno++]);
		else if (cur_par == "-ge")
			gap_ext_base = atof(argv[argno++]);
		else if (cur_par == "-tgo")
			gap_term_open_base = atof(argv[argno++]);
		else if (cur_par == "-tge")
			gap_term_ext_base = atof(argv[argno++]);
		
		else if (cur_par == "-gsd")
			scaler_div = atoi(argv[argno++]);
		else if (cur_par == "-gsl")
			scaler_log = atoi(argv[argno++]);
		else if (cur_par == "-dgr")
			enable_gap_rescaling = false;
		else if (cur_par == "-dgo")
			enable_gap_optimization = false;
		else if (cur_par == "-dsp")
			enable_total_score_calculation = false;
		else if (cur_par == "-r")
			n_refinements = atoi(argv[argno++]);
		else if (cur_par == "-fr")
			enable_auto_refinement = false;
		else if (cur_par == "-t")
			n_threads = atoi(argv[argno++]);
		else if (cur_par == "-ri")
			thr_internal_refinement = atoi(argv[argno++]);
		else if (cur_par == "-gt")
		{
			string name(argv[argno++]);
			gt_method = GT::fromString(name);
			if (gt_method == GT::imported) {
				guide_tree_in_file = argv[argno++];
			}
			else if (gt_method == GT::chained) {
				guide_tree_seed = atoi(argv[argno++]);
			}
		}
		else if (cur_par == "-dist") {
			distance = str2dist(argv[argno++]);
		}
		else if (cur_par == "-parttree") {
			gt_heuristic = GT::PartTree;
		}
		else if (cur_par == "-medoidtree") {
			gt_heuristic = GT::ClusterTree;
		}
		else if (cur_par == "-medoid_threshold") {
			heuristic_threshold = atoi(argv[argno++]);
		}
		else if (cur_par == "-subtree_size") {
			subtree_size = atoi(argv[argno++]);
		}
		else if (cur_par == "-sample_size") {
			sample_size = atoi(argv[argno++]);
		}
		else if (cur_par == "-cluster_fraction") {
			cluster_fraction = (float) atof(argv[argno++]);
		}
		else if (cur_par == "-cluster_iters") {
			cluster_iters = atoi(argv[argno++]);
		}
		else if (cur_par == "-gt_export") {
			export_tree = true;
		}
		else if (cur_par == "-dist_export") {
			export_distances = true;
		}
		else if (cur_par == "-square_matrix") {
			generate_square_matrix = true;
		}
		else if (cur_par == "-pid") {
			calculate_pid = true;
		}
		else if (cur_par == "-gz")
		{
			gzippd_output = true;
		}
		else if (cur_par == "-gz-lev")
		{
			int g_lev = atoi(argv[argno++]);
			if (g_lev < 0 || g_lev > 12)
			{
				LOG_NORMAL << "Incorrect gzip level: " << g_lev << " was changed to default value: " << gzip_level << endl;
				g_lev = gzip_level;
			}

			gzip_level = g_lev;
		}
#ifdef DEVELOPER_MODE
		else if (cur_par == "-shuffle") {
			shuffle = atoi(argv[argno++]);
		}
		else if (cur_par == "-ref") {
			test_ref_sequenes = true;
			ref_file_name = argv[argno++];
		}
#endif
		else if (cur_par == "-v")
			verbose_mode = true;
		else if (cur_par == "-vv")
			very_verbose_mode = true;
		else
		{
			LOG_NORMAL << "Unknown parameter: " << cur_par << "\n";
			return false;
		}
	}

	if (argno + 2 > argc)
	{
		LOG_NORMAL << "No file name gives\n";
		return false;
	}

	input_file_name = string(argv[argno++]);
	output_file_name = string(argv[argno++]);


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