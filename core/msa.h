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
#ifdef DEVELOPER_MODE
	const int monte_carlo_trials = 1000;
#endif

protected:
	CParams params;
	static double SM_MIQS[24][24];
	vector<vector<score_t>> score_matrix;
	vector<score_t> score_vector;

	instruction_set_t instruction_set;
	vector<CSequence> sequences;

	vector<CGappedSequence> gapped_sequences;

	vector<pair<int, int>> guide_tree;
	vector<pair<int, score_t>> gt_stats;
	vector<size_t> prof_cardinalities;
	score_t avg_sim;

	map<size_t, CProfile*> profiles;
	CProfile *final_profile;

	mt19937 rnd_rfn;
	uint32_t n_threads;

	CLCSBP_Classic lcsbp_classic0;

	set<int> already_refined;

	CStopWatch timers[4];

	vector<CSequence> ref_sequences;

#ifdef DEBUG_MODE
	double estimated_identity;
#endif
	
	void DetermineInstructionSet();
	void init_sm();

#ifdef DEVELOPER_MODE
	double GetLCS(CSequence &seq1, CSequence &seq2);
#endif

	// Guide tree construction method
	virtual void UPGMA();
	inline uint64_t UPGMA_TriangleSubscript(uint64_t uIndex1, uint64_t uIndex2);
	virtual void UPGMA_CalculateDistances(UPGMA_dist_t *dist_matrix);
	virtual void SingleLinkage();
#ifdef DEVELOPER_MODE
	void GuideTreeChained();
#endif
	bool ImportGuideTreeFromNewick();
	bool ExportGuideTreeToNewick();
	bool ExportDistanceMatrix(float*matrix, size_t size, const std::string& fname);

#ifdef DEVELOPER_MODE
	bool LoadRefSequences();
	uint64_t CalculateRefSequencesSubTreeSize(double *monte_carlo_subtree_size);
	int SubTreeSize(set<int> &seq_ids);
#endif

	void RefineRandom(CProfile* profile_to_refine, vector<size_t> &dest_prof_id);
	void RefineMostEmptyAndFullColumn(CProfile *profile_to_refine, vector<size_t> &dest_prof_id, vector<size_t> &gap_stats, bool valid_gap_stats);

public:
	CFAMSA();
	~CFAMSA();

	bool SetSequences(vector<CSequence> &_sequences);
	bool SetSequences(vector<CSequence> &&_sequences);
	bool SetParams(CParams &_params);

	bool ComputeGuideTree();
	bool ComputeAlignment();
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