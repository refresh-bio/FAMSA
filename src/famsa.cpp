/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#include <string>
#include <iostream>
#include <numeric>

#include "./core/io_service.h"
#include "./core/params.h"
#include "./msa.h"
#include "./utils/timer.h"
#include "./utils/log.h"

#undef min
#undef max


// ****************************************************************************
int main(int argc, char *argv[])
{
	Log::getInstance(Log::LEVEL_NORMAL).enable();
	
	LOG_NORMAL << "FAMSA (Fast and Accurate Multiple Sequence Alignment) ver. " << FAMSA_VER << "\n"
		<< "  by " << FAMSA_AUTHORS << " (" << FAMSA_DATE << ")\n\n";
	
	bool showExpert = 0;
	CParams params;
	
	if(!params.parse(argc, argv, showExpert)) {
		// some parameters could be parsed - used default values for printing
		CParams def_params;
		def_params.show_usage(showExpert);
		return 0;
	}
	
	CStopWatch timer, timer_saving;

	timer.StartTimer();

	if (params.verbose_mode) {
		Log::getInstance(Log::LEVEL_VERBOSE).enable();
	}
	if (params.very_verbose_mode) {
		Log::getInstance(Log::LEVEL_VERBOSE).enable();
		Log::getInstance(Log::LEVEL_DEBUG).enable();
	}

	// ***** Read input file
	LOG_VERBOSE << "Processing: " << params.input_file_name << "\n";

	vector<CGappedSequence*> result;
	vector<CSequence> sequences;

	size_t input_seq_cnt = IOService::loadFasta(params.input_file_name, sequences);
    if(input_seq_cnt == 0){			
    	LOG_NORMAL << "Error: no (or incorrect) input file\n";
		return -1;
	} else if (input_seq_cnt == 1){
        CGappedSequence resultSeq(std::move(sequences[0]));
		result.push_back(&resultSeq);
		return IOService::saveAlignment(params.output_file_name, result, params.n_threads, -1);
	}

	// ***** Load sequences to FAMSA
	CFAMSA famsa(params);

	if(!famsa.ComputeMSA(sequences))
	{
		LOG_NORMAL << "Some interal error occured!\n";
		return -1;
	}

	timer_saving.StartTimer();

	if (famsa.GetAlignment(result)) {
		LOG_VERBOSE << "Saving alignment in " << params.output_file_name;		
		if(params.gzippd_output)
			IOService::saveAlignment(params.output_file_name, result, params.n_threads, params.gzip_level);
		else
			IOService::saveAlignment(params.output_file_name, result, params.n_threads, -1);
	//		IOService::saveAlignment(params.output_file_name, result);
		LOG_VERBOSE << " [OK]" << endl;		
	}

	timer_saving.StopTimer();
	timer.StopTimer();

	LOG_VERBOSE << " Alignment saving                                 : " << timer_saving.GetElapsedTime() << "s\n";

	LOG_VERBOSE << "Total computation time: " << timer.GetElapsedTime() << "s\n";
	LOG_NORMAL << "Done!\n";

	if (params.verbose_mode || params.very_verbose_mode) {
		FILE* out = fopen("famsa.stats", "at");
		fprintf(out, "time.save=%lf\n", timer_saving.GetElapsedTime());
		fclose(out);
	}

	return 0;
}

