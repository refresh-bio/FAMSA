/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#include "../core/queues.h"
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

	// Calculate Sackin index for guide tree
	sackin_index = 0;
	for (int i = 0; i < gapped_sequences->size(); ++i)
		sackin_index += prof_depth[i] + 1;

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

// *******************************************************************
// Return Sackin index 
// *******************************************************************
uint64_t CProfileQueue::GetSackinIndex()
{
	return sackin_index;
}

// *******************************************************************
// CSingleLinkageQueue
// *******************************************************************
CSingleLinkageQueue::CSingleLinkageQueue(vector<CSequence> *_sequences, uint32_t _n_rows, uint32_t _max_buffered_rows)
{
	sequences = _sequences;
	n_rows = _n_rows;
	max_buffered_rows = min(n_rows, _max_buffered_rows);

	sim_vector_2d.resize(max_buffered_rows);
	for (auto &x : sim_vector_2d)
		x.resize(n_rows);

	ready_rows.resize(n_rows, make_pair(-1, false));
	
	lowest_uncomputed_row = 0;

	for (int i = 0; i < max_buffered_rows; ++i)
		available_buffers.push(i);

	eoq_flag = false;
}

// *******************************************************************
CSingleLinkageQueue::~CSingleLinkageQueue()
{	
}

// *******************************************************************
bool CSingleLinkageQueue::GetTask(int &row_id, vector<CSequence> *&_sequences, vector<double> *&sim_vector)
{
	unique_lock<mutex> lck(mtx);
	cv_tasks.wait(lck, [this]{return !this->available_buffers.empty() || this->eoq_flag; });

	if (eoq_flag)
		return false;	// End of data in the profiles queue

	row_id = lowest_uncomputed_row++;

	if (lowest_uncomputed_row >= n_rows)
		eoq_flag = true;

	_sequences = sequences;

	int buffer_row_id = available_buffers.top();
	available_buffers.pop();

	sim_vector = &sim_vector_2d[buffer_row_id];

	ready_rows[row_id].first = buffer_row_id;

#ifdef PRODUCE_LOG
	cerr << "GetTask : " << row_id << "\n";
	cerr << available_buffers.size() << "  " << lowest_uncomputed_row << "\n";
#endif

	return true;
}

// *******************************************************************
void CSingleLinkageQueue::RegisterSolution(int row_id)
{
	unique_lock<mutex> lck(mtx);

	ready_rows[row_id].second = true;

#ifdef PRODUCE_LOG
	cerr << "Registered : " << row_id << "\n";
	cerr << available_buffers.size() << "  " << lowest_uncomputed_row << "\n";
#endif

	cv_rows.notify_one();
}

// *******************************************************************
bool CSingleLinkageQueue::GetSolution(int row_id, vector<double> *&sim_vector)
{
	unique_lock<mutex> lck(mtx);
	cv_rows.wait(lck, [this, row_id]{return this->ready_rows[row_id].second; });

	int buffer_row_id = ready_rows[row_id].first;

	sim_vector = &sim_vector_2d[buffer_row_id];

#ifdef PRODUCE_LOG
	cerr << "GetSol : " << row_id << "\n";
	cerr << available_buffers.size() << "  " << lowest_uncomputed_row << "\n";
#endif

	return true;
}

// *******************************************************************
void CSingleLinkageQueue::ReleaseSolution(int row_id)
{
	unique_lock<mutex> lck(mtx);

	int buffer_row_id = ready_rows[row_id].first;

	available_buffers.push(buffer_row_id);

#ifdef PRODUCE_LOG
	cerr << "Release : " << row_id << "\n";
	cerr << available_buffers.size() << "  " << lowest_uncomputed_row << "\n";
#endif

	cv_tasks.notify_all();
}


// *******************************************************************
// Queue for UPGMA
// *******************************************************************
CUPGMAQueue::CUPGMAQueue(vector<CSequence> *_sequences, uint32_t _n_rows, UPGMA_dist_t *_dist_matrix)
{
	sequences = _sequences;
	n_rows = _n_rows;
	lowest_uncomputed_row = 0;
	eoq_flag = false;
	dist_matrix = _dist_matrix;
}

CUPGMAQueue::~CUPGMAQueue()
{
}

inline long long CUPGMAQueue::TriangleSubscript(long long uIndex1, long long uIndex2)
{
	long long v;

	if (uIndex1 >= uIndex2)
		v = uIndex2 + (uIndex1*(uIndex1 - 1)) / 2;
	else
		v = uIndex1 + (uIndex2*(uIndex2 - 1)) / 2;

	return v;
}

bool CUPGMAQueue::GetTask(int &row_id, vector<CSequence> *&_sequences, UPGMA_dist_t *&dist_row)
{
	unique_lock<mutex> lck(mtx);

	if (eoq_flag)
		return false;	// End of data in the profiles queue

	row_id = lowest_uncomputed_row++;

	if (lowest_uncomputed_row >= n_rows)
		eoq_flag = true;

	_sequences = sequences;

	dist_row = dist_matrix + TriangleSubscript(row_id, 0);

	return true;
}

