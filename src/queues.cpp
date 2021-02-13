/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#include "queues.h"
#include <iostream>

//#define PRODUCE_LOG

// *******************************************************************
// CProfileQueue
// *******************************************************************
CProfileQueue::CProfileQueue(vector<CGappedSequence> *_gapped_sequences, map<size_t, CProfile*> *_profiles, vector<pair<int, int>> *_guide_tree)
{
	gapped_sequences = _gapped_sequences;
	profiles = _profiles;
	guide_tree = _guide_tree;

	eoq_flag = false;

	prof_depth.assign(guide_tree->size(), 0);
	for (size_t i = guide_tree->size() - 1; i >= gapped_sequences->size(); --i)
	{
		prof_depth[(*guide_tree)[i].first] = prof_depth[i] + 1;
		prof_depth[(*guide_tree)[i].second] = prof_depth[i] + 1;
	}

	// Insert all sequence to profile converstions as ready to process
	vector<size_t> init_ids;
	for (size_t i = 0; i < gapped_sequences->size(); ++i)
		init_ids.push_back(i);

	for (size_t i = 0; i < gapped_sequences->size(); ++i)
		pq.push(make_pair(prof_depth[i], i));

	// Number of child profiles ready for each parent profile
	ready_profiles.assign(guide_tree->size(), 0);

	child_parent_mapping.assign(guide_tree->size(), 0);
	for (size_t i = 0; i < guide_tree->size(); ++i)
	{
		int id1 = (*guide_tree)[i].first;
		int id2 = (*guide_tree)[i].second;

		if (id1 == -1)
			continue;

		child_parent_mapping[id1] = i;
		child_parent_mapping[id2] = i;
	}

	counter = 0;
}

// *******************************************************************
CProfileQueue::~CProfileQueue()
{
	// Nothing to do
}

// *******************************************************************
bool CProfileQueue::GetTask(size_t &prof_id, CGappedSequence *&gs, CProfile *&prof1, CProfile *&prof2)
{
	unique_lock<mutex> lck(mtx);
	cv.wait(lck, [this]{return !this->pq.empty() || this->eoq_flag; });

	if (eoq_flag)
		return false;	// End of data in the profiles queue

	prof_id = pq.top().second;
	pq.pop();

	if ((*guide_tree)[prof_id].first == -1)
	{
		gs = &((*gapped_sequences)[prof_id]);
		prof1 = nullptr;
		prof2 = nullptr;
	}
	else
	{
		gs = nullptr;
		prof1 = (*profiles)[(*guide_tree)[prof_id].first];
		prof2 = (*profiles)[(*guide_tree)[prof_id].second];
	}

	return true;
}

// *******************************************************************
void CProfileQueue::AddSolution(size_t prof_id, CProfile *prof)
{
	lock_guard<mutex> lck(mtx);

	if ((*guide_tree)[prof_id].first == -1)	// Just construct profile from a sequence
		(*profiles)[prof_id] = prof;
	else
	{
		// Add new profile and remove old profiles
		(*profiles)[prof_id] = prof;

		profiles->erase((*guide_tree)[prof_id].first);
		profiles->erase((*guide_tree)[prof_id].second);
	}

	if (++ready_profiles[child_parent_mapping[prof_id]] == 2)		// Profile is ready to be computed as both child profiles are already computed
		pq.push(make_pair(prof_depth[prof_id], child_parent_mapping[prof_id]));

	if (ready_profiles[0] == 1)		// final profile was computed
		eoq_flag = true;

	cv.notify_all();
}
