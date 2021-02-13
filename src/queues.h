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

#include "defs.h"
#include "sequence.h"
#include "profile.h"

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

public:
	CProfileQueue(vector<CGappedSequence> *_gapped_sequences, map<size_t, CProfile*> *_profiles, vector<pair<int, int>> *_guide_tree);
	~CProfileQueue();

	bool GetTask(size_t &prof_id, CGappedSequence *&gs, CProfile *&prof1, CProfile *&prof2);
	void AddSolution(size_t prof_id, CProfile *prof);
};


#endif
