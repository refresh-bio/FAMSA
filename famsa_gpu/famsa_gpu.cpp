/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

Version: 1.0
Date   : 2016-03-11
*/

#include <string>
#include <iostream>

#include "../core/input_file.h"
#include "../core/output_file.h"
#include "../core/msa.h"
#include "../core_gpu/gpumsa.h"
#include "../core/timer.h"

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

	bool enable_gap_rescaling;
	bool enable_gap_optimization;
	bool enable_total_score_calculation;
	bool enable_auto_refinement;
	bool verbose_mode;
	bool very_verbose_mode;

	bool use_gpu;					// internal param
	int gpu_platform;
	int gpu_device;
	int gpu_pairs_per_wave;			// internal param
	int gpu_pairs_per_task;			// internal param
	int gpu_threads_per_pair;		// internal param

	int n_threads;
	string input_file_name;
	string output_file_name;
} execution_params_t;

using namespace std;

execution_params_t execution_params;

void init_params();
bool parse_params(int argc, char **argv);
void show_usage();
void set_famsa_params(CParams &famsa_params);

// ****************************************************************************
// Show command-line parameters
void show_usage()
{
	init_params();		// To set default param values

	cout << "FAMSA (Fast and Accurate Multiple Sequence Alignment) ver. " << FAMSA_VER << " CPU and GPU\n";
	cout << "  by " << FAMSA_AUTHORS << " (" << FAMSA_DATE << ")\n";
	cout << "Usage: famsa-gpu [parameters] <input_file_name> <output_file_name>\n";
	cout << "where available parameters are:\n";
	cout << "  -go <value> - gap open cost (default: " << execution_params.gap_open << ")\n";
	cout << "  -ge <value> - gap extension cost (default: " << execution_params.gap_ext << ")\n";
	cout << "  -tgo <value> - terminal gap open cost (default: " << execution_params.gap_term_open << ")\n";
	cout << "  -tge <value> - terminal gap extenstion cost (default: " << execution_params.gap_term_ext << ")\n";
	cout << "  -gsd <value> - gap cost scaler div-term (default: " << execution_params.gap_scaler_div << ")\n";
	cout << "  -gsl <value> - gap cost scaler log-term (default: " << execution_params.gap_scaler_log << ")\n";

	cout << "  -dgr - disable gap cost rescaling (default: enabled)\n";
	cout << "  -dgo - disable gap optimization (default: enabled)\n";
	cout << "  -dsp - disable sum of pairs optimization during refinement (default: enabled)\n";

	cout << "  -r <value> - no. of refinement iterations (default: " << execution_params.n_refinements << ")\n";
	cout << "  -fr - disable auto refinement turning off (for sets larger than " << execution_params.thr_refinement << " seq.)\n";
	cout << "  -t <value> - no. of threads, 0 means all available (default: " << execution_params.n_threads << ")\n";
	cout << "  -v - verbose mode, show timing information (default: disabled)\n";

	cout << "  -gpu_p <value> - gpu platform id \n";
	cout << "  -gpu_d <value> - gpu device id \n";
	cout << "     Hint: both -gpu_p and -gpu_d must be specified to enable gpu support\n";
}

// ****************************************************************************
// Set default execution parameter values
void init_params()
{
	execution_params.gap_open = 14.85;
	execution_params.gap_ext = 1.25;
	execution_params.gap_term_open = 0.66;
	execution_params.gap_term_ext = 0.66;
	execution_params.gap_scaler_div = 7;
	execution_params.gap_scaler_log = 45;
	execution_params.guided_alignment_radius = 50;

	execution_params.enable_gap_rescaling			= true;
	execution_params.enable_gap_optimization		= true;
	execution_params.enable_total_score_calculation = true;
	execution_params.enable_auto_refinement			= true;

	execution_params.n_refinements					= 100;
	execution_params.thr_refinement					= 1000;
	execution_params.use_gpu						= false;
	execution_params.gpu_platform					= -1;
	execution_params.gpu_device						= -1;
	execution_params.gpu_pairs_per_task				= 64;		
	execution_params.gpu_pairs_per_wave				= 2 << 20;	
	execution_params.gpu_threads_per_pair			= 8;		
	execution_params.n_threads						= 1;
	execution_params.input_file_name				= "";
	execution_params.output_file_name				= "";
	execution_params.verbose_mode					= false;
	execution_params.very_verbose_mode = false;
}

// ****************************************************************************
// **** Parses parameters
bool parse_params(int argc, char **argv)
{
	if(argc < 3)
		return false;

	int argno;

	for(argno = 1; argno + 2 < argc;)
	{
		string cur_par = argv[argno++];

		if(cur_par == "-go")
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
		else if (cur_par == "-v")
			execution_params.verbose_mode = true;
		else if (cur_par == "-gpu_p")
			execution_params.gpu_platform = atoi(argv[argno++]);
		else if(cur_par == "-gpu_d")
			execution_params.gpu_device = atoi(argv[argno++]);
		else
		{
			cout << "Unknown parameter: " << cur_par << "\n";
			return false;
		}
	}

	if(execution_params.gpu_platform >= 0 && execution_params.gpu_device >= 0)
		execution_params.use_gpu = true;

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
	famsa_params.n_threads						= execution_params.n_threads;
	famsa_params.verbose_mode					= execution_params.verbose_mode;
}

// ****************************************************************************
int main(int argc, char *argv[])
{
	init_params();
	if(!parse_params(argc, argv))
	{
		show_usage();
		return 0;
	}

	std::shared_ptr<CFAMSA> ptr;
	if (execution_params.use_gpu) {
		ptr = std::make_shared<CGpuFAMSA>(execution_params.gpu_platform, execution_params.gpu_device, 
			execution_params.gpu_pairs_per_wave, execution_params.gpu_pairs_per_task, execution_params.gpu_threads_per_pair);
	} else 
		ptr = std::make_shared<CFAMSA>();

	CParams params;
	CFAMSA& famsa = *ptr;

	set_famsa_params(params);

	// ***** Read input file
	CInputFile in_file;
	COutputFile out_file;

	if(!in_file.ReadFile(execution_params.input_file_name))
	{
		cout << "Error: no (or incorrect) input file\n";
		return 0;
	}

	vector<CSequence> sequences;

	in_file.StealSequences(sequences);

	// ***** Load sequences to FAMSA
	if(!famsa.SetSequences(std::move(sequences)))
		return 0;

	if(!famsa.SetParams(params))
	{
		cout << "Error: No input sequences\n";
		return 0;
	}

	if(!famsa.ComputeMSA())
	{
		cout << "Some interal error occured!\n";
		return 0;
	}

	vector<CGappedSequence*> result;
	famsa.GetAlignment(result);

	out_file.PutSequences(std::move(result));
	out_file.SaveFile(execution_params.output_file_name);

	return 0;
}

