/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#ifndef _MSA_H
#define _MSA_H

#include <string>
#include <vector>
#include <set>
#include <random>
#include <thread>
#include <mutex>
#include <condition_variable>

#include "../core/queues.h"
#include "../core/sequence.h"
#include "../core/profile.h"
#include "../core/params.h"
#include "../core/timer.h"

#include "../core/lcsbp.h"
#include "../core/lcsbp_classic.h"

class CFAMSA 
{
protected:
	static const int TIMER_SORTING = 0;
	static const int TIMER_TREE_BUILD = 1;
	static const int TIMER_ALIGNMENT = 2;
	static const int TIMER_REFINMENT = 3;
	static const int TIMER_TREE_STORE = 4;

	CParams params;
	static double SM_MIQS[24][24];
	vector<vector<score_t>> score_matrix;
	vector<score_t> score_vector;

	instruction_set_t instruction_set;
	vector<CSequence> sequences;
	vector<CGappedSequence> gapped_sequences;

	score_t avg_sim;

	map<size_t, CProfile*> profiles;
	CProfile *final_profile;

	mt19937 rnd_rfn;
	uint32_t n_threads;

	set<int> already_refined;

	CStopWatch timers[4];

	vector<CSequence> ref_sequences;

#ifdef DEBUG_MODE
	double estimated_identity;
#endif
	
	void DetermineInstructionSet();
	void init_sm();

#ifdef DEVELOPER_MODE
	bool LoadRefSequences();
#endif

	void RefineRandom(CProfile* profile_to_refine, vector<size_t> &dest_prof_id);
	void RefineMostEmptyAndFullColumn(CProfile *profile_to_refine, vector<size_t> &dest_prof_id, vector<size_t> &gap_stats, bool valid_gap_stats);

public:
	CFAMSA();
	~CFAMSA();

	bool SetSequences(vector<CSequence> &_sequences);
	bool SetSequences(vector<CSequence> &&_sequences);
	bool SetParams(CParams &_params);

	bool ComputeAlignment(std::vector<std::pair<int,int>>& guide_tree);
#ifdef DEBUG_MODE
	bool RefineAlignment(string output_file_name = "");
#else
	bool RefineAlignment(CProfile *&profile_to_refine);
#endif

	bool GetAlignment(vector<CGappedSequence*> &result);
	
	score_t GetScore()
	{
		if(final_profile != nullptr)
			return final_profile->total_score;
		else
			return 0;
	};

#ifdef DEBUG_MODE
	double GetEstimatedIdentity()
	{
		return estimated_identity;
	};
#endif

	bool ComputeMSA();
};


#endif