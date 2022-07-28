/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

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
	vector<CGappedSequence*> *gapped_sequences;
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
	CProfileQueue(vector<CGappedSequence*> *_gapped_sequences, map<size_t, CProfile*> *_profiles, vector<pair<int, int>> *_guide_tree, uint32_t _max_no_threads);
	~CProfileQueue();

	bool GetTask(size_t &prof_id, CGappedSequence *&gs, CProfile *&prof1, CProfile *&prof2, uint32_t &no_threads, uint32_t &no_rows_per_block);
	void AddSolution(size_t prof_id, CProfile *prof);
};

// *******************************************************************************************
// Multithreading queue with registering mechanism:
//   * The queue can report whether it is in wainitng for new data state or there will be no new data
template<typename T> class CPriorityQueue
{
	typedef map<size_t, T> queue_t;

	queue_t q;
	bool is_completed;
	int n_producers;
	uint32_t n_elements;
	size_t current_priority;

	mutable mutex mtx;								// The mutex to synchronise on
	condition_variable cv_queue_empty;

public:
	typename queue_t::iterator q_it;

	// *******************************************************************************************
	CPriorityQueue(const int _n_producers)
	{
		current_priority = 0;

		Restart(_n_producers);
	};

	// *******************************************************************************************
	~CPriorityQueue()
	{};

	// *******************************************************************************************
	void Restart(const int _n_producers)
	{
		unique_lock<mutex> lck(mtx);

		is_completed = false;
		n_producers = _n_producers;
		n_elements = 0;
		current_priority = 0;
	}

	// *******************************************************************************************
	bool IsEmpty()
	{
		lock_guard<mutex> lck(mtx);
		return n_elements == 0;
	}

	// *******************************************************************************************
	bool IsCompleted()
	{
		lock_guard<mutex> lck(mtx);

		return n_elements == 0 && n_producers == 0;
	}

	// *******************************************************************************************
	void MarkCompleted()
	{
		lock_guard<mutex> lck(mtx);
		n_producers--;

		if (!n_producers)
			cv_queue_empty.notify_all();
	}

	// *******************************************************************************************
	void Push(const T data, const size_t priority)
	{
		unique_lock<mutex> lck(mtx);

		//bool was_empty = n_elements == 0;
		q.emplace(priority, data);
		++n_elements;

		//		if (was_empty)
		cv_queue_empty.notify_all();
	}

	// *******************************************************************************************
	void Emplace(const size_t priority, T&& data)
	{
		unique_lock<mutex> lck(mtx);

		//		bool was_empty = n_elements == 0;
		q.emplace(priority, move(data));
		++n_elements;

		cv_queue_empty.notify_all();
	}

	// *******************************************************************************************
	bool Pop(T& data)
	{
		unique_lock<mutex> lck(mtx);
		cv_queue_empty.wait(lck, [this] {return (!this->q.empty() && current_priority == q.begin()->first) || !this->n_producers; });

		if (n_elements == 0)
			return false;

		data = move(q.begin()->second);

		q.erase(q.begin());
		--n_elements;
		++current_priority;

		cv_queue_empty.notify_all();

		return true;
	}

	// *******************************************************************************************
	size_t GetSize()
	{
		unique_lock<mutex> lck(mtx);

		return n_elements;
	}
};

// *******************************************************************************************
// Multithreading queue with registering mechanism:
//   * The queue can report whether it is in wainitng for new data state or there will be no new data
template<typename T> class CLimitedPriorityQueue
{
	typedef map<size_t, T> queue_t;

	queue_t q;
	bool is_completed;
	int n_producers;
	uint32_t n_elements;
	size_t current_priority;
	size_t size_limit;

	mutable mutex mtx;								// The mutex to synchronise on
	condition_variable cv_queue_empty;
	condition_variable cv_queue_full;

public:
	typename queue_t::iterator q_it;

	// *******************************************************************************************
	CLimitedPriorityQueue(const int _n_producers, const size_t _size_limit)
	{
		current_priority = 0;

		Restart(_n_producers, _size_limit);
	};

	// *******************************************************************************************
	~CLimitedPriorityQueue()
	{};

	// *******************************************************************************************
	void Restart(const int _n_producers, size_t _size_limit)
	{
		unique_lock<mutex> lck(mtx);

		is_completed = false;
		n_producers = _n_producers;
		size_limit = _size_limit;
		n_elements = 0;
		current_priority = 0;
	}

	// *******************************************************************************************
	bool IsEmpty()
	{
		lock_guard<mutex> lck(mtx);
		return n_elements == 0;
	}

	// *******************************************************************************************
	bool IsCompleted()
	{
		lock_guard<mutex> lck(mtx);

		return n_elements == 0 && n_producers == 0;
	}

	// *******************************************************************************************
	void MarkCompleted()
	{
		lock_guard<mutex> lck(mtx);
		n_producers--;

		if (!n_producers)
			cv_queue_empty.notify_all();
	}

	// *******************************************************************************************
	void Push(const T data, const size_t priority)
	{
		unique_lock<mutex> lck(mtx);
		cv_queue_full.wait(lck, [this, &priority] {return this->q.size() < this->size_limit || priority == this->current_priority; });

		//bool was_empty = n_elements == 0;
		q.emplace(priority, data);
		++n_elements;

		//		if (was_empty)
		cv_queue_empty.notify_all();
		cv_queue_full.notify_all();
	}

	// *******************************************************************************************
	void Emplace(const size_t priority, T&& data)
	{
		unique_lock<mutex> lck(mtx);
		cv_queue_full.wait(lck, [this, &priority] {return this->q.size() < this->size_limit || priority == this->current_priority; });

		//		bool was_empty = n_elements == 0;
		q.emplace(priority, move(data));
		++n_elements;

		cv_queue_empty.notify_all();
		cv_queue_full.notify_all();
	}

	// *******************************************************************************************
	bool Pop(T& data)
	{
		unique_lock<mutex> lck(mtx);
		cv_queue_empty.wait(lck, [this] {return (!this->q.empty() && current_priority == q.begin()->first) || !this->n_producers; });

		if (n_elements == 0)
			return false;

		data = move(q.begin()->second);

		q.erase(q.begin());
		--n_elements;
		++current_priority;

		cv_queue_empty.notify_all();
		cv_queue_full.notify_all();

		return true;
	}

	// *******************************************************************************************
	size_t GetSize()
	{
		unique_lock<mutex> lck(mtx);

		return n_elements;
	}
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


// ************************************************************************************
// Multithreading queue with registering mechanism:
//   * The queue can report whether it is in wainitng for new data state or there will be no new data
template<typename T> class RegisteringQueue
{
	typedef queue<T, deque<T>> queue_t;

	queue_t q;
	bool is_completed;
	int n_producers;
	int n_elements;
	int size;

	mutable mutex mtx;								// The mutex to synchronise on
	condition_variable cv_queue_empty;

public:



	// *****************************************************************************************
	//
	RegisteringQueue(int _n_producers, int size = 0)
	{
		Restart(_n_producers, size);
	};

	// *****************************************************************************************
	//
	~RegisteringQueue()
	{};

	// *****************************************************************************************
	//
	void Restart(int _n_producers, int size = 0)
	{
		unique_lock<mutex> lck(mtx);

		is_completed = false;
		n_producers = _n_producers;
		n_elements = 0;
		this->size = size;
	}

	// *****************************************************************************************
	//
	bool IsEmpty()
	{
		lock_guard<mutex> lck(mtx);
		return n_elements == 0;
	}

	// *****************************************************************************************
	//
	bool IsCompleted()
	{
		lock_guard<mutex> lck(mtx);

		return n_elements == 0 && n_producers == 0;
	}

	// *****************************************************************************************
	//
	void MarkCompleted()
	{
		lock_guard<mutex> lck(mtx);
		n_producers--;

		if (!n_producers)
			cv_queue_empty.notify_all();
	}

	// *****************************************************************************************
	//
	void Push(T data)
	{
		unique_lock<mutex> lck(mtx);

		if (size > 0) {
			cv_queue_empty.wait(lck, [this] {return this->n_elements < size; });
		}

		bool was_empty = n_elements == 0;
		q.push(data);
		++n_elements;

		if (was_empty)
			cv_queue_empty.notify_all();
	}

	// *****************************************************************************************
	//
	void PushRange(vector<T>& data)
	{
		unique_lock<mutex> lck(mtx);
		bool was_empty = n_elements == 0;

		for (auto& x : data)
			q.push(x);
		n_elements += data.size();

		if (was_empty)
			cv_queue_empty.notify_all();
	}

	// *****************************************************************************************
	//
	bool Pop(T& data)
	{
		unique_lock<mutex> lck(mtx);
		cv_queue_empty.wait(lck, [this] {return !this->q.empty() || !this->n_producers; });

		if (n_elements == 0)
			return false;

		data = q.front();
		q.pop();
		--n_elements;
		if (n_elements == 0)
			cv_queue_empty.notify_all();

		return true;
	}

	// *****************************************************************************************
	//
	uint32_t GetSize()
	{
		return n_elements;
	}
};


#endif
