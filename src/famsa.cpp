/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#include <string>
#include <iostream>
#include <numeric>
#include <fstream>

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

	memory_monotonic_safe mma(16 << 20, 64);
	
	vector<CGappedSequence*> result;
	vector<CSequence> sequences;

	size_t input_seq_cnt = IOService::loadFasta(params.input_file_name, sequences, &mma);
	
	bool ok = true;
	
	if (input_seq_cnt == 0) {
		// no sequences loaded - signal error
		LOG_NORMAL << "Error: no (or incorrect) input file\n";
		ok = false;
	}
	else if (input_seq_cnt == 1) {
		// single input sequence - write without alignment
		CGappedSequence resultSeq(std::move(sequences[0]));
		result.push_back(&resultSeq);
		ok = IOService::saveAlignment(params.output_file_name, result, params.n_threads, -1);
	}
	else {
		// multitple input sequences - run alignment
		CFAMSA famsa(params);
		famsa.getStatistics().put("input.n_sequences", sequences.size());

		if (famsa.ComputeMSA(sequences)) {
			timer_saving.StartTimer();

			// Save alignment if it was generated
			if (famsa.GetAlignment(result)) {

				famsa.getStatistics().put("alignment.length", result[0]->gapped_size);

				LOG_VERBOSE << "Saving alignment in " << params.output_file_name;
				if (params.gzippd_output)
					ok = IOService::saveAlignment(params.output_file_name, result, params.n_threads, params.gzip_level);
				else
					ok = IOService::saveAlignment(params.output_file_name, result, params.n_threads, -1);

				LOG_VERBOSE << " [OK]" << endl;
			}

			timer_saving.StopTimer();
			timer.StopTimer();
			LOG_NORMAL << "Done!\n";

			if (params.verbose_mode || params.very_verbose_mode) {
				famsa.getStatistics().put("time.save", timer_saving.GetElapsedTime());
				famsa.getStatistics().put("time.total", timer.GetElapsedTime());

				string stats = famsa.getStatistics().toString();

				LOG_VERBOSE << endl << endl << "Statistics:" << endl << stats << endl;

				std::ofstream ofs("famsa.stats");
				ofs << "[stats]" << endl << stats;
				ofs.close();
			}

		}
		else {
			LOG_NORMAL << "Some interal error occured!\n";
			ok = false;
		}
	}

	sequences.clear();

	return ok ? 0 : -1;
}

