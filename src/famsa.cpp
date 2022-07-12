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
	
	memory_monotonic_safe mma(16 << 20, 64);

	memory_monotonic_safe mma1(16 << 20, 64);

	memory_monotonic_safe mma2(16 << 20, 64);

	vector<CGappedSequence*> aligned_result;

	if (params.profile_aligning)
	  {
	    vector<CGappedSequence> first_gapped_seq_set;
	    vector<CSequence> first_seq_set;

	    int seq_no=0;

	    size_t aligned_cnt_1=IOService::loadProfile(params.input_prof1, first_seq_set, first_gapped_seq_set, seq_no, &mma1);

	    vector<CGappedSequence*> test;

	    string seq;


   	    vector<CGappedSequence> second_gapped_seq_set;
	    vector<CSequence> second_seq_set;

	    seq_no= static_cast<int>(aligned_cnt_1);

	    size_t aligned_cnt_2=IOService::loadProfile(params.input_prof2, second_seq_set, second_gapped_seq_set, seq_no, &mma2);
	    CFAMSA profile_aligner(params);

	    profile_aligner.adjustParams((int)(aligned_cnt_1+aligned_cnt_2));	  

	    CGappedSequence* gs;

	    vector<CGappedSequence>* first_gapped_sequence_pointer=&first_gapped_seq_set;
	    vector<CGappedSequence>* second_gapped_sequence_pointer=&second_gapped_seq_set;

	    CProfile prof_1=CProfile(&profile_aligner.params);

	    CProfile prof_2=CProfile(&profile_aligner.params);

	    for (int i=0; i<aligned_cnt_1;i++){
	      gs=&((*first_gapped_sequence_pointer)[i]);
	      prof_1.AppendRawSequence(*gs);
	    }
	    std::cout<<"appended the raw sequences"<<'\n';
	    prof_1.CalculateCountersScores();
	    std::cout<<"calculated counter and scores"<<'\n';

	    for (int i=0; i<aligned_cnt_2; i++){
	      gs=&((*second_gapped_sequence_pointer)[i]);
	      prof_2.AppendRawSequence(*gs);
	    }
	    std::cout<<"appended the raw sequences"<<'\n';
	    prof_2.CalculateCountersScores();
	    std::cout<<"calculated the counter and scores"<<'\n';
	    uint32_t no_threads=1;
	    uint32_t no_rows_per_box=0;

	    CProfile merged_profile=CProfile(&prof_1, &prof_2, &profile_aligner.params, no_threads, no_rows_per_box);

	    std::cout<<"merged profiles"<<'\n';

	    aligned_result=merged_profile.data;

	    sort(aligned_result.begin(), aligned_result.end(), [](CGappedSequence *p, CGappedSequence *q){return p-> id < q->id;});
	    IOService::saveAlignment(params.output_prof_align, aligned_result, params.n_threads, -1);
	    return 0;
	  }	


	LOG_VERBOSE << "Processing"<<params.input_file_name<<'\n';
	    
	
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

