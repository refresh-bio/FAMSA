/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

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

#include "./core/queues.h"
#include "./core/sequence.h"
#include "./core/profile.h"
#include "./core/params.h"
#include "./utils/timer.h"
#include "./utils/statistics.h"

#include "./lcs/lcsbp.h"
#include "./lcs/lcsbp_classic.h"

class AbstractTreeGenerator;

class CFAMSA 
{
protected:
	static const int TIMER_SORTING = 0;
	static const int TIMER_TREE_BUILD = 1;
	static const int TIMER_ALIGNMENT = 2;
	static const int TIMER_REFINMENT = 3;
	static const int TIMER_TREE_STORE = 4;

 	static double SM_MIQS[24][24];

	CParams params;
	instruction_set_t instruction_set;
	
	vector<vector<score_t>> score_matrix;
	vector<score_t> score_vector;

	vector<CGappedSequence> gapped_sequences;

	map<size_t, CProfile*> profiles;
	CProfile *final_profile;

	mt19937 rnd_rfn;

	set<int> already_refined;

	CStopWatch timers[5];

	Statistics statistics;

#ifdef DEBUG_MODE
	double estimated_identity;
#endif
	
	void initScoreMatrix();
 
#ifdef DEVELOPER_MODE
	vector<CSequence> ref_sequences;
	bool LoadRefSequences();
#endif

	void RefineRandom(CProfile* profile_to_refine, vector<size_t> &dest_prof_id);
	void RefineMostEmptyAndFullColumn(CProfile *profile_to_refine, vector<size_t> &dest_prof_id, vector<size_t> &gap_stats, bool valid_gap_stats);

	std::shared_ptr<AbstractTreeGenerator> createTreeGenerator(const CParams& params);
	
	void sortAndExtendSequences(std::vector<CSequence>& sequences);
	void extendSequences(std::vector<CSequence>& sequences);
	void shrinkSequences(std::vector<CSequence>& sequences);
	void removeDuplicates(std::vector<CSequence*>& sorted_seqs, std::vector<int>& original2sorted);

public:
	
	
	CFAMSA(CParams& _params);
	~CFAMSA();

	CProfile* ComputeAlignment(std::vector<CGappedSequence*>& gapped_sequences, tree_structure& guide_tree);
#ifdef DEBUG_MODE
	bool RefineAlignment(string output_file_name = "");
#else
	bool RefineAlignment(CProfile *&profile_to_refine);
#endif

	bool GetAlignment(vector<CGappedSequence*> &result);
  void adjustParams(int n_seqs);
	
	score_t GetScore() { return final_profile != nullptr ? final_profile->total_score : 0;  }
	
	const Statistics& getStatistics() const { return statistics;  }
	Statistics& getStatistics() { return statistics; }

#ifdef DEBUG_MODE
	double GetEstimatedIdentity()
	{
		return estimated_identity;
	};
#endif

	bool ComputeMSA(vector<CSequence>& sequences);

	bool alignProfiles(vector<CGappedSequence>& p1, vector<CGappedSequence>& p2);
};


#endif
