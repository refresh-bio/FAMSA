#pragma once

#include <functional>
namespace pooled_threads
{
	class ThreadPoolTask;

	class thread
	{
		bool _joinable = false;
		ThreadPoolTask* task;
		void Create(std::function<void()>&& f);
	public:
		thread(const thread&) = delete;
		thread& operator=(const thread&) = delete;

		thread(thread&& rhs) :
			_joinable(rhs._joinable), task(rhs.task)
		{
			rhs._joinable = false;
			rhs.task = nullptr;
		}
		thread& operator=(thread&& rhs)
		{
			_joinable = rhs._joinable;
			task = rhs.task;
			rhs._joinable = false;
			rhs.task = nullptr;
			return *this;
		}

		template<typename _Callable, typename... _Args>
		explicit thread(_Callable&& __f, _Args&&... __args) :
			_joinable(true)
		{
			Create(std::bind(std::forward<_Callable>(__f), std::forward<_Args>(__args)...));
		}

		bool joinable()
		{
			return _joinable;
		}

		void join();

		~thread();
	};
}