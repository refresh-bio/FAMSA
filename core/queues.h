/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#ifndef _QUEUES_H
#define _QUEUES_H

#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <map>
#include <list>
#include <stack>

#include "../core/defs.h"
#include "../core/sequence.h"
#include "../core/profile.h"

class CProfileQueue
{
	vector<CGappedSequence> *gapped_sequences;
	map<size_t, CProfile*> *profiles;
	vector<pair<int, int>> *guide_tree;

	vector<pair<int, int>>::iterator gt_iter;

	vector<size_t> ready_profiles;
	vector<size_t> child_parent_mapping;

	vector<size_t> prof_depth;
	priority_queue<pair<int32_t, int32_t>> pq;

	queue<size_t, list<size_t>> q;

	bool eoq_flag;

	mutex mtx;
	condition_variable cv;

	size_t counter;

	uint64_t sackin_index;

public:
	CProfileQueue(vector<CGappedSequence> *_gapped_sequences, map<size_t, CProfile*> *_profiles, vector<pair<int, int>> *_guide_tree);
	~CProfileQueue();

	bool GetTask(size_t &prof_id, CGappedSequence *&gs, CProfile *&prof1, CProfile *&prof2);
	void AddSolution(size_t prof_id, CProfile *prof);

	uint64_t GetSackinIndex();
};

class CSingleLinkageQueue
{
	vector<CSequence> *sequences;

	vector<vector<double>> sim_vector_2d;
	//	vector<vector<uint32_t>> lcs_len_2d;

	vector<pair<int, bool>> ready_rows;
	stack<int, vector<int>> available_buffers;

	uint32_t lowest_uncomputed_row;
	uint32_t n_rows;
	uint32_t max_buffered_rows;

	bool eoq_flag;

	mutex mtx;
	condition_variable cv_tasks, cv_rows;

public:
	CSingleLinkageQueue(vector<CSequence> *_sequences, uint32_t _n_rows, uint32_t _max_buffered_rows);
	~CSingleLinkageQueue();

	bool GetTask(int &row_id, vector<CSequence> *&_sequences, vector<double> *&sim_vector);
	void RegisterSolution(int row_id);
	bool GetSolution(int row_id, vector<double> *&sim_vector);
	void ReleaseSolution(int row_id);
};


class CUPGMAQueue
{
	vector<CSequence> *sequences;

	vector<pair<int, bool>> ready_rows;
	stack<int, vector<int>> available_buffers;

	uint32_t lowest_uncomputed_row;
	uint32_t n_rows;
	uint32_t max_buffered_rows;
	UPGMA_dist_t *dist_matrix;

	bool eoq_flag;

	mutex mtx;
	condition_variable cv_tasks, cv_rows;

	inline long long TriangleSubscript(long long uIndex1, long long uIndex2);

public:
	CUPGMAQueue(vector<CSequence> *_sequences, uint32_t _n_rows, UPGMA_dist_t *_dist_matrix);
	~CUPGMAQueue();

	bool GetTask(int &row_id, vector<CSequence> *&_sequences, UPGMA_dist_t *&dist_row);
};

#endif
