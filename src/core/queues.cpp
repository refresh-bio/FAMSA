/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#include "../core/queues.h"
#include <iostream>

//#define PRODUCE_LOG

// *******************************************************************
// CProfileQueue
// *******************************************************************
CProfileQueue::CProfileQueue(vector<CGappedSequence*> *_gapped_sequences, map<size_t, CProfile*> *_profiles, vector<pair<int, int>> *_guide_tree, uint32_t _max_no_threads)
{
	gapped_sequences = _gapped_sequences;
	profiles = _profiles;
	guide_tree = _guide_tree;
	max_no_threads = _max_no_threads;
	no_working_threads = 0;

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
		init_ids.emplace_back(i);

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
}

// *******************************************************************
CProfileQueue::~CProfileQueue()
{
	// Nothing to do
}

// *******************************************************************
void CProfileQueue::CheckAlignInParallel(CProfile* prof1, CProfile* prof2, uint32_t& no_threads, uint32_t& no_rows_per_box)
{
	const uint32_t min_box_width_per_thread = 512;
	//	const uint32_t min_box_width_per_thread = 12;

	uint32_t no_available_threads = max_no_threads - no_working_threads;

	//uint32_t min_prof_width = min(prof1->width, prof2->width);
	uint32_t max_prof_width = (uint32_t) max(prof1->width, prof2->width);

	if(no_available_threads == 1)
	{
		no_threads = 1;
		no_rows_per_box = 0;

/*		cout << "max_no_threads: " + to_string(max_no_threads) +
			"  no_working_threads: " + to_string(no_working_threads) +
			"  no_threads: " + to_string(no_threads) + "\n";*/

		return;
	}

	if (max_prof_width < 2 * min_box_width_per_thread)
	{
		no_threads = 1;
		no_rows_per_box = 0;
		
/*		cout << "max_no_threads: " + to_string(max_no_threads) +
			"  no_working_threads: " + to_string(no_working_threads) +
			"  no_threads: " + to_string(no_threads) + "\n";*/

		return;
	}

	uint32_t est_no_threads = max(1u, no_available_threads / (uint32_t) (pq.size() + 1u));

	no_threads = min(est_no_threads, max_prof_width / min_box_width_per_thread);
	no_threads = max(no_threads, 1u);
	//uint32_t box_width_per_thread = max_prof_width / no_threads;

	if (no_threads > 1)
	{
//		if (box_width_per_thread < 512)
			no_rows_per_box = 4;
/*		else if (box_width_per_thread < 1024)
			no_rows_per_box = 3;
		else if (box_width_per_thread < 2048)
			no_rows_per_box = 2;
		else 
			no_rows_per_box = 1;*/
	}
	else
		no_rows_per_box = 0;

/*	cout << "max_no_threads: " + to_string(max_no_threads) +
		"  no_working_threads: " + to_string(no_working_threads) +
		"  no_threads: " + to_string(no_threads) + "\n";*/
}

// *******************************************************************
bool CProfileQueue::GetTask(size_t &prof_id, CGappedSequence *&gs, CProfile *&prof1, CProfile *&prof2, uint32_t& no_threads, uint32_t& no_rows_per_block)
{
	unique_lock<mutex> lck(mtx);
	cv.wait(lck, [this]{return !this->pq.empty() || this->eoq_flag; });

	if (eoq_flag)
		return false;	// End of data in the profiles queue

	prof_id = pq.top().second;
	pq.pop();

	if ((*guide_tree)[prof_id].first == -1)
	{
		gs = (*gapped_sequences)[prof_id];
		prof1 = nullptr;
		prof2 = nullptr;
		no_threads = 1;
		no_rows_per_block = 0;
	}
	else
	{
		gs = nullptr;
		prof1 = (*profiles)[(*guide_tree)[prof_id].first];
		prof2 = (*profiles)[(*guide_tree)[prof_id].second];

		CheckAlignInParallel(prof1, prof2, no_threads, no_rows_per_block);
	}

	no_working_threads += no_threads;
	m_reserved_threads[prof_id] = no_threads;

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

	no_working_threads -= m_reserved_threads[prof_id];
	m_reserved_threads.erase(prof_id);

	cv.notify_all();
}
