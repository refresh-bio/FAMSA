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

	uint32_t max_no_threads;
	map<size_t, uint32_t> m_reserved_threads;
	uint32_t no_working_threads;

	vector<pair<int, int>>::iterator gt_iter;

	vector<size_t> ready_profiles;
	vector<size_t> child_parent_mapping;

	vector<size_t> prof_depth;
	priority_queue<pair<int32_t, int32_t>> pq;

	queue<size_t, list<size_t>> q;

	bool eoq_flag;

	mutex mtx;
	condition_variable cv;

	void CheckAlignInParallel(CProfile *prof1, CProfile *prof2, uint32_t& no_threads, uint32_t& no_rows_per_box);

public:
	CProfileQueue(vector<CGappedSequence> *_gapped_sequences, map<size_t, CProfile*> *_profiles, vector<pair<int, int>> *_guide_tree, uint32_t _max_no_threads);
	~CProfileQueue();

	bool GetTask(size_t &prof_id, CGappedSequence *&gs, CProfile *&prof1, CProfile *&prof2, uint32_t &no_threads, uint32_t &no_rows_per_block);
	void AddSolution(size_t prof_id, CProfile *prof);
};


// *****************************************************************************************
//
class CBarrier
{
public:
	CBarrier(const CBarrier&) = delete;
	CBarrier& operator=(const CBarrier&) = delete;
	explicit CBarrier(unsigned int count) :
		m_count(count), m_generation(0),
		m_count_reset_value(count)
	{
	}
	void arrive_and_wait()
	{
		std::unique_lock< std::mutex > lock(m_mutex);
		unsigned int gen = m_generation;
		if (--m_count == 0)
		{
			m_generation++;
			m_count = m_count_reset_value;
			m_cond.notify_all();
			return;
		}
		while (gen == m_generation)
			m_cond.wait(lock);
	}
private:
	std::mutex m_mutex;
	std::condition_variable m_cond;
	unsigned int m_count;
	unsigned int m_generation;
	unsigned int m_count_reset_value;
};

#endif
