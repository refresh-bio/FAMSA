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

			string inputs[2] {params.input_file_name, params.input_file_name_2};
			vector<CGappedSequence> profiles[2];
			size_t sizes[2];

			for (int i = 0; i < 2; ++i) {
				auto& profile = profiles[i];
				sizes[i] = IOService::loadFasta(inputs[i], profile, &mma);
				bool correct = (sizes[i] > 0) && std::all_of(profile.begin(), profile.end(), [&profile](const CGappedSequence& seq) {
					return seq.gapped_size == profile.front().gapped_size;
				});

				if (!correct) {
					throw std::runtime_error("Incorrect profile");
				}
			}

			CFAMSA profile_aligner(params);

			profile_aligner.adjustParams((int)(sizes[0] + sizes[1]));
			profile_aligner.alignProfiles(profiles[0], profiles[1]);

			profile_aligner.GetAlignment(result);

			IOService::saveAlignment(params.output_file_name, result, params.n_threads, 
				params.gzippd_output ? params.gzip_level : -1, params.remove_rare_columns ? params.rare_column_threshold : 1.0f);
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
						params.gzippd_output ? params.gzip_level : -1, params.remove_rare_columns ? params.rare_column_threshold : 1.0f);

					LOG_VERBOSE << " [OK]" << endl;
				}

				timer_saving.StopTimer();
				timer.StopTimer();
				LOG_NORMAL << "Done!\n";

				if (params.areStatsStored()) {
					famsa.getStatistics().put("time.save", timer_saving.GetElapsedTime());
					famsa.getStatistics().put("time.total", timer.GetElapsedTime());
					string stats = famsa.getStatistics().toString();

					if (params.verbose_mode || params.very_verbose_mode) {
						LOG_VERBOSE << endl << endl << "Statistics:" << endl << stats << endl;
					}

					if (params.stats_file_name.length()) {
						std::ofstream ofs(params.stats_file_name);
						ofs << "[stats]" << endl << stats;
						ofs.close();
					}
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

