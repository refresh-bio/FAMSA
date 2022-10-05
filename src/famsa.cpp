/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#include <string>
#include <iostream>
#include <numeric>
#include <fstream>
#include <stdexcept>

#include "./core/io_service.h"
#include "./core/params.h"
#include "./msa.h"
#include "./utils/timer.h"
#include "./utils/log.h"

#include "./core/version.h"

#define VAL(str) #str
#define TOSTRING(str) VAL(str)

#undef min
#undef max


// ****************************************************************************
int main(int argc, char *argv[])
{
	bool ok = true;

	try {

		Log::getInstance(Log::LEVEL_NORMAL).enable();

		LOG_NORMAL << "FAMSA (Fast and Accurate Multiple Sequence Alignment) \n" 
			<< "  version " << FAMSA_VER 
		#ifdef GIT_COMMIT
			<< "-" << TOSTRING(GIT_COMMIT) 
		#endif
			<< " (" << FAMSA_DATE << ")\n"
			<< "  " << FAMSA_AUTHORS << "\n\n";

		bool showExpert = 0;
		CParams params;

		if (!params.parse(argc, argv, showExpert)) {
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

		memory_monotonic_safe mma(16 << 20, 64);
		vector<CGappedSequence*> result;
		vector<CSequence> sequences;

		// profile - profile alignment
		if (params.profile_aligning) {
			LOG_VERBOSE << "Aligning " << params.input_file_name << " with " << params.input_file_name_2 << "\n";

			vector<CGappedSequence> profile1;
			vector<CGappedSequence> profile2;
	
			size_t size1 = IOService::loadFasta(params.input_file_name, profile1, &mma);
			size_t size2 = IOService::loadFasta(params.input_file_name_2, profile2, &mma);
			CFAMSA profile_aligner(params);

			profile_aligner.adjustParams((int)(size1 + size2));
			profile_aligner.alignProfiles(profile1, profile2);

			profile_aligner.GetAlignment(result);

			IOService::saveAlignment(params.output_file_name, result, params.n_threads, 
				params.gzippd_output ? params.gzip_level : -1);
			return 0;
		}		

		LOG_VERBOSE << "Aligning " << params.input_file_name << "\n";

		size_t input_seq_cnt = IOService::loadFasta(params.input_file_name, sequences, &mma);

		if (input_seq_cnt == 0) {
			// no sequences loaded - signal error
			throw(std::runtime_error("No (or incorrect) input file."));
		}
		else {
			// at least one input sequences - run alignment
			CFAMSA famsa(params);
			famsa.getStatistics().put("input.n_sequences", sequences.size());

			if (famsa.ComputeMSA(sequences)) {
				timer_saving.StartTimer();

				// Save alignment if it was generated
				if (famsa.GetAlignment(result)) {

					famsa.getStatistics().put("alignment.length", result[0]->gapped_size);

					LOG_VERBOSE << "Saving alignment in " << params.output_file_name;
					ok = IOService::saveAlignment(params.output_file_name, result, params.n_threads, 
						params.gzippd_output ? params.gzip_level : -1);

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
				throw(std::runtime_error("Some interal error occured"));
			}
		}

		sequences.clear();
	}
	catch (std::runtime_error& err) {
		LOG_NORMAL << endl << "[ERROR] " << err.what() << endl;
		ok = false;
	}

	return ok ? 0 : -1;
}

