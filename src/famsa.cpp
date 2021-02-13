/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

 */

#include <string>
#include <iostream>
#include <numeric>

#include "io_service.h"
#include "msa.h"
#include "timer.h"
#include "log.h"

// #include <R.h>
#include <Rinternals.h>
// #include <R_ext/Rdynload.h>
#undef min
#undef max

typedef struct {
    double gap_open;
    double gap_ext;
    double gap_term_open;
    double gap_term_ext;
    int gap_scaler_div;
    int gap_scaler_log;
    int guided_alignment_radius;	// internal param
    int n_refinements;
    int thr_refinement;
    int thr_internal_refinement;

    bool enable_gap_rescaling;
    bool enable_gap_optimization;
    bool enable_total_score_calculation;
    bool enable_auto_refinement;
    bool verbose_mode;
    bool very_verbose_mode;

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
    double indel_exp;

    bool test_ref_sequenes;
    string ref_file_name;

    int n_threads;
    string input_file_name;
    string output_file_name;

    int64_t shuffle;

} execution_params_t;

using namespace std;

execution_params_t execution_params;

void init_params();
bool parse_params(int argc, char **argv, bool& showExpert);
void show_usage(bool expert);
void set_famsa_params(CParams &famsa_params);

// ****************************************************************************
// Show command-line parameters
void show_usage(bool expert)
{
    init_params();		// To set default param values

    cerr
        << "Usage:\n"
        << "  famsa [options] <input_file> <output_file>\n\n"

    << "Positional parameters:\n"
    << "  input_file - input file in FASTA format (pass STDIN when reading from standard input)\n"
    << "  output_file - output file (pass STDOUT when writing to standard output); available outputs:\n"
    << "      * alignment in FASTA format,\n"
    << "      * guide tree in Newick format (-gt_export option specified),\n"
    << "      * distance matrix in CSV format (-dist_export option specified),\n\n"

    << "Options:\n"
    << "  -help - show advanced options\n"
    << "  -t <value> - no. of threads, 0 means all available (default: " << execution_params.n_threads << ")\n"
    << "  -v - verbose mode, show timing information (default: disabled)\n\n"

    << "  -gt <sl | upgma | import <file>> - guide tree method (default: sl):\n"
    << "      * sl - single linkage\n"
    << "      * upgma - UPGMA\n"
    << "      * import <file> - imported from a Newick file\n"
    //<< "  -medoidtree - use MedoidTree heuristic for speeding up tree construction (default: disabled)\n"
    //<< "  -parttree - use PartTree heuristic for speeding up tree construction (default: disabled)\n"
    << "  -gt_export - export a guide tree to output file in Newick format\n"
    << "  -dist_export - export a distance matrix to output file in CSV format\n\n";

    if (expert) {
        cerr << "Advanced options:\n"
            << "  -r <value> - no. of refinement iterations (default: " << execution_params.n_refinements << ")\n"
            << "  -fr - force refinement (by default the refinement is disabled for sets larger than " << execution_params.thr_refinement << " seq.)\n"
            << "  -go <value> - gap open cost (default: " << execution_params.gap_open << ")\n"
            << "  -ge <value> - gap extension cost (default: " << execution_params.gap_ext << ")\n"
            << "  -tgo <value> - terminal gap open cost (default: " << execution_params.gap_term_open << ")\n"
            << "  -tge <value> - terminal gap extenstion cost (default: " << execution_params.gap_term_ext << ")\n"

            << "  -gsd <value> - gap cost scaller div-term (default: " << execution_params.gap_scaler_div << ")\n"
            << "  -gsl <value> - gap cost scaller log-term (default: " << execution_params.gap_scaler_log << ")\n"
            << "  -dgr - disable gap cost rescaling (default: enabled)\n"
            << "  -dgo - disable gap optimization (default: enabled)\n"
            << "  -dsp - disable sum of pairs optimization during refinement (default: enabled)\n";
#ifdef DEVELOPER_MODE
        cerr << "  -ref <file_name> - load referential sequences (for benchmarks) and calculate the minimal subtree size containing them\n"
            << "  -vv - very verbose mode, show timing information (default: disabled)\n";
#endif
        cerr << endl;
    }
}

// ****************************************************************************
// Set default execution parameter values
void init_params()
{
    execution_params.gap_open						= 14.85;
    execution_params.gap_ext						= 1.25;
    execution_params.gap_term_open					= 0.66;
    execution_params.gap_term_ext					= 0.66;
    execution_params.gap_scaler_div					= 7;
    execution_params.gap_scaler_log					= 45;
    execution_params.guided_alignment_radius		= 50;

    execution_params.enable_gap_rescaling			= true;
    execution_params.enable_gap_optimization		= true;
    execution_params.enable_total_score_calculation = true;
    execution_params.enable_auto_refinement			= true;

    execution_params.n_refinements					= 100;
    execution_params.thr_refinement					= 1000;
    execution_params.thr_internal_refinement		= 0;

    execution_params.gt_method						= GT::SLINK;
    execution_params.gt_heuristic					= GT::None;
    execution_params.heuristic_threshold			= 0;

    execution_params.guide_tree_seed				= 0;
    execution_params.subtree_size					= 100;
    execution_params.sample_size					= 2000;
    execution_params.cluster_fraction				= 0.1;
    execution_params.cluster_iters					= 2;

    execution_params.guide_tree_in_file				= "";
    execution_params.export_tree					= false;
    execution_params.export_distances				= false;

    execution_params.test_ref_sequenes				= false;
    execution_params.ref_file_name					= "";

    execution_params.n_threads						= 0;
    execution_params.input_file_name				= "";
    execution_params.output_file_name				= "";
    execution_params.verbose_mode					= false;
    execution_params.very_verbose_mode				= false;

    execution_params.indel_exp						= 1.0;
    execution_params.shuffle						= -1;
}

// ****************************************************************************
// **** Parses parameters
bool parse_params(int argc, char **argv, bool& showExpert)
{
    // ugly workaround
    if (argc == 2 && string(argv[1]) == "-help") {
        showExpert = true;
        return false;
    }

    if(argc < 3)
        return false;

    int argno;

    for(argno = 1; argno + 2 < argc;)
    {
        string cur_par = argv[argno++];

        if (cur_par == "-help") {
            showExpert = true;
            return false;
        }
        else if	(cur_par == "-go")
            execution_params.gap_open = atof(argv[argno++]);
        else if(cur_par == "-ge")
            execution_params.gap_ext = atof(argv[argno++]);
        else if(cur_par == "-tgo")
            execution_params.gap_term_open = atof(argv[argno++]);
        else if(cur_par == "-tge")
            execution_params.gap_term_ext = atof(argv[argno++]);
        else if(cur_par == "-gsd")
            execution_params.gap_scaler_div = atoi(argv[argno++]);
        else if(cur_par == "-gsl")
            execution_params.gap_scaler_log = atoi(argv[argno++]);
        else if(cur_par == "-dgr")
            execution_params.enable_gap_rescaling = false;
        else if(cur_par == "-dgo")
            execution_params.enable_gap_optimization = false;
        else if(cur_par == "-dsp")
            execution_params.enable_total_score_calculation = false;
        else if(cur_par == "-r")
            execution_params.n_refinements = atoi(argv[argno++]);
        else if (cur_par == "-fr")
            execution_params.enable_auto_refinement = false;
        else if (cur_par == "-t")
            execution_params.n_threads = atoi(argv[argno++]);
        else if (cur_par == "-ie")
            execution_params.indel_exp = atof(argv[argno++]);
        else if (cur_par == "-ri")
            execution_params.thr_internal_refinement = atoi(argv[argno++]);
        else if (cur_par == "-gt")
        {
            string name(argv[argno++]);
            execution_params.gt_method = GT::fromString(name);
            if (execution_params.gt_method == GT::imported) {
                execution_params.guide_tree_in_file = argv[argno++];
            }
            else if (execution_params.gt_method == GT::chained) {
                execution_params.guide_tree_seed = atoi(argv[argno++]);
            }
        }
        else if (cur_par == "-parttree") {
            execution_params.gt_heuristic = GT::PartTree;
        }
        else if (cur_par == "-medoidtree") {
            execution_params.gt_heuristic = GT::ClusterTree;
        }
        else if (cur_par == "-medoid_threshold") {
            execution_params.heuristic_threshold = atoi(argv[argno++]);
        }
        else if (cur_par == "-subtree_size") {
            execution_params.subtree_size = atoi(argv[argno++]);
        }
        else if (cur_par == "-sample_size") {
            execution_params.sample_size = atoi(argv[argno++]);
        }
        else if (cur_par == "-cluster_fraction") {
            execution_params.cluster_fraction = atof(argv[argno++]);
        }
        else if (cur_par == "-cluster_iters") {
            execution_params.cluster_iters = atoi(argv[argno++]);
        }
        else if (cur_par == "-gt_export") {
            execution_params.export_tree = true;
        }
        else if (cur_par == "-dist_export") {
            execution_params.export_distances = true;
        }
#ifdef DEVELOPER_MODE
        else if (cur_par == "-shuffle") {
            execution_params.shuffle = atoi(argv[argno++]);
        }
        else if (cur_par == "-ref") {
            execution_params.test_ref_sequenes = true;
            execution_params.ref_file_name = argv[argno++];
        }
#endif
        else if (cur_par == "-v")
            execution_params.verbose_mode = true;
        else if (cur_par == "-vv")
            execution_params.very_verbose_mode = true;
        else
        {
            cout << "Unknown parameter: " << cur_par << "\n";
            return false;
        }
    }

    if(argno + 2 > argc)
    {
        cout << "No file name gives\n";
        return false;
    }

    execution_params.input_file_name  = string(argv[argno++]);
    execution_params.output_file_name = string(argv[argno++]);

    return true;
}

// ****************************************************************************
void set_famsa_params(CParams &famsa_params)
{
    famsa_params.gap_open						= round(-cost_cast_factor * execution_params.gap_open);
    famsa_params.gap_ext						= round(-cost_cast_factor * execution_params.gap_ext);
    famsa_params.gap_term_open					= round(-cost_cast_factor * execution_params.gap_term_open);
    famsa_params.gap_term_ext					= round(-cost_cast_factor * execution_params.gap_term_ext);
    famsa_params.scaler_div						= execution_params.gap_scaler_div;
    famsa_params.scaler_log						= execution_params.gap_scaler_log;
    famsa_params.guided_alignment_radius		= execution_params.guided_alignment_radius;
    famsa_params.enable_gap_optimization		= execution_params.enable_gap_optimization;
    famsa_params.enable_gap_rescaling			= execution_params.enable_gap_rescaling;
    famsa_params.enable_total_score_calculation	= execution_params.enable_total_score_calculation;
    famsa_params.enable_auto_refinement			= execution_params.enable_auto_refinement;
    famsa_params.n_refinements					= execution_params.n_refinements;
    famsa_params.thr_refinement					= execution_params.thr_refinement;
    famsa_params.thr_internal_refinement		= execution_params.thr_internal_refinement;
    famsa_params.n_threads						= execution_params.n_threads;
    famsa_params.verbose_mode					= execution_params.verbose_mode;
    famsa_params.very_verbose_mode				= execution_params.very_verbose_mode;

    famsa_params.indel_exp						= execution_params.indel_exp;

    famsa_params.gt_method						= execution_params.gt_method;
    famsa_params.gt_heuristic					= execution_params.gt_heuristic;
    famsa_params.heuristic_threshold			= execution_params.heuristic_threshold;
    famsa_params.subtree_size					= execution_params.subtree_size;
    famsa_params.sample_size					= execution_params.sample_size;
    famsa_params.cluster_fraction				= execution_params.cluster_fraction;
    famsa_params.cluster_iters					= execution_params.cluster_iters;

    famsa_params.guide_tree_in_file				= execution_params.guide_tree_in_file;
    famsa_params.export_tree					= execution_params.export_tree;
    famsa_params.export_distances				= execution_params.export_distances;
    famsa_params.output_file_name				= execution_params.output_file_name;

    famsa_params.test_ref_sequences				= execution_params.test_ref_sequenes;
    famsa_params.ref_file_name					= execution_params.ref_file_name;

    famsa_params.shuffle						= execution_params.shuffle;
}

// ****************************************************************************
extern "C" {


    int famsaCPP(int *argc, char **argv)
    {
        cerr << "FAMSA (Fast and Accurate Multiple Sequence Alignment) ver. " << FAMSA_VER << " CPU\n"
            << "  by " << FAMSA_AUTHORS << " (" << FAMSA_DATE << ")\n\n";

        bool showExpert = 0;
        init_params();
        if(!parse_params(*argc, argv, showExpert))
        {
            show_usage(showExpert);
            return 0;
        }

        CStopWatch timer;

        timer.StartTimer();

        CParams params;
        CFAMSA famsa;

        set_famsa_params(params);

        Log::getInstance(Log::LEVEL_NORMAL).enable();
        if (params.verbose_mode) {
            Log::getInstance(Log::LEVEL_VERBOSE).enable();
        }
        if (params.very_verbose_mode) {
            Log::getInstance(Log::LEVEL_VERBOSE).enable();
            Log::getInstance(Log::LEVEL_DEBUG).enable();
        }

        // ***** Read input file
        LOG_VERBOSE << "Processing: " << execution_params.input_file_name << "\n";

        vector<CGappedSequence*> result;
        vector<CSequence> sequences;

        size_t input_seq_cnt = IOService::loadFasta(execution_params.input_file_name, sequences);
        if(input_seq_cnt == 0){
            LOG_NORMAL << "Error: no (or incorrect) input file\n";
            return 1;
        } else if (input_seq_cnt == 1){
            CGappedSequence resultSeq(sequences[0]);
            result.push_back(&resultSeq);
            return IOService::saveAlignment(execution_params.output_file_name, result);
        }

        // ***** Load sequences to FAMSA
        if(!famsa.SetSequences(std::move(sequences)))
            return 1;

        if(!famsa.SetParams(params))
        {
            LOG_NORMAL << "Error: No input sequences\n";
            return 1;
        }

        if(!famsa.ComputeMSA())
        {
            LOG_NORMAL << "Some interal error occured!\n";
            return 1;
        }

        if (famsa.GetAlignment(result)) {
            LOG_VERBOSE << "Saving alignment in " << execution_params.output_file_name;
            IOService::saveAlignment(execution_params.output_file_name, result);
            LOG_VERBOSE << " [OK]" << endl;
        }

        timer.StopTimer();

        LOG_VERBOSE << "Total computation time: " << timer.GetElapsedTime() << "s\n";
        LOG_NORMAL << "Done!\n";

        return 0;
    }

    // Register the native code in R.

    R_CMethodDef cMethods[] = {
        {"famsaCPP", (DL_FUNC) &famsaCPP, 2},
        NULL
    };

    void R_init_famsa(DllInfo *info)
    {
        R_registerRoutines(info, cMethods, NULL, NULL, NULL);
    }
}
