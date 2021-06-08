#include "pooled_threads.h"
#include <exception>
#include <functional>
#include <thread>
#include <mutex>
#include <vector>
#include <condition_variable>

#include <iostream>

//define to create separate thread pool for each thread (than each thread is a producer of tasks)
//#define USE_THREAD_LOCAL 

namespace pooled_threads
{
	class ThreadPoolTask
	{
		int id;
		std::function<void()> fun;
		std::mutex mtx;
		std::condition_variable cv_join;
		bool done = false;
	public:
		int GetId() const
		{
			return id;
		}
		ThreadPoolTask(int id) : id(id) {}

		void SetFunction(std::function<void()>&& _f)
		{
			done = false;

			//I am not sure it it is legal, but if it is it is possible that it would work a little faster than move assignment operator
			//this->fun.~function();
			//new (&this->fun) std::function<void()>(move(_f));

			this->fun = std::move(_f);
		}
		void operator()()
		{
			fun();
		}

		void Join()
		{
			std::unique_lock<std::mutex> lck(mtx);
			cv_join.wait(lck, [this] {return done; });
		}

		void NotifyDone()
		{
			std::lock_guard<std::mutex> lck(mtx);
			done = true;
			cv_join.notify_one();
		}
	};


	class ThreadPoolTaskQueue
	{
		std::vector<ThreadPoolTask*> tasks;
		std::mutex mtx;
		std::condition_variable cv_pop;
		bool finished = false;

		unsigned long pos = 0;
	public:
		void SetMaxTasks(unsigned long m)
		{
			tasks.resize(m);
		}
		bool PopTask(ThreadPoolTask*& task)
		{
			std::unique_lock<std::mutex> lck(mtx);
			cv_pop.wait(lck, [this] {return pos || finished; });
			if (!pos)
				return false;
			task = tasks[--pos];
			return true;
		}
		void AddTask(ThreadPoolTask& task)
		{
			std::lock_guard<std::mutex> lck(mtx);
			tasks[pos++] = &task;
			if (pos == 1)
				cv_pop.notify_all();
		}

		void Finish()
		{
			std::lock_guard<std::mutex> lck(mtx);
			finished = true;
			cv_pop.notify_all();
		}
	};

	class ThreadPool
	{
		std::vector<std::thread> threads;
		std::vector<std::unique_ptr<ThreadPoolTask>> tasks;
		std::vector<int> free_tasks; //free tasks ids
		ThreadPoolTaskQueue task_queue;
		std::mutex mtx;
	public:
		void AddTask(ThreadPoolTask& task)
		{
			task_queue.AddTask(task);
		}

		ThreadPoolTask* GetFreeTask()
		{
#ifndef USE_THREAD_LOCAL
			std::lock_guard<std::mutex> lck(mtx);
#endif
			if (!free_tasks.size())
			{
				tasks.emplace_back(std::make_unique<ThreadPoolTask>(tasks.size()));
				task_queue.SetMaxTasks(static_cast<unsigned long>(tasks.size()));
				free_tasks.emplace_back(static_cast<int>(tasks.size()) - 1);
				threads.emplace_back([this]()
				{
					ThreadPoolTask* task;
					while (task_queue.PopTask(task))
					{
						(*task)();
						task->NotifyDone();
					}
				});
			}
			auto r = tasks[free_tasks.back()].get();
			free_tasks.pop_back();
			return r;
		}

		void ReturnTask(ThreadPoolTask* task)
		{
#ifndef USE_THREAD_LOCAL
			std::lock_guard<std::mutex> lck(mtx);
#endif
			free_tasks.push_back(task->GetId());
		}

		~ThreadPool()
		{
			//std::cout << "ThreadPool dtor: \n";
			//std::cout << "n threads: " << threads.size() << "\n";
			//std::cout << "n tasks: " << tasks.size() << "\n";
			task_queue.Finish();
			for (auto& th : threads)
				if (th.joinable())
					th.join();
		}
	};

	//Global Instance of thread pool
#ifdef USE_THREAD_LOCAL 
	thread_local ThreadPool thread_pool;
#else
	ThreadPool thread_pool;
#endif


	void thread::join()
	{
		if (task && _joinable)
		{
			task->Join();
			thread_pool.ReturnTask(task);
			task = nullptr;
			_joinable = false;
		}
	}

	thread::~thread()
	{
		if (task && _joinable)
			std::terminate();
	}

	void thread::Create(std::function<void()>&& f)
	{
		task = thread_pool.GetFreeTask();
		task->SetFunction(std::move(f));
		thread_pool.AddTask(*task);
	}
}