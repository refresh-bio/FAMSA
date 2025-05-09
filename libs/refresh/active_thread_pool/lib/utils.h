#pragma once

#include <atomic>
#include <cinttypes>

namespace refresh
{
	namespace utils
	{
		// *********************************************************************
		inline void noop()
		{
#if defined(_WIN32)
			_mm_pause();
#elif defined(__aarch64__)
			std::this_thread::yield();
#else
			__builtin_ia32_pause();
			//			_mm_pause();
#endif
		}
	}

	// *********************************************************************
	class spin_mutex
	{
		std::atomic_flag a_lock{};

	public:
		// *********************************************************************
		void lock()
		{
			while (a_lock.test_and_set(std::memory_order_acquire))
				while (a_lock.test(std::memory_order_relaxed))
					;
		}

		// *********************************************************************
		void lock_with_noop()
		{
			while (a_lock.test_and_set(std::memory_order_acquire))
				while (a_lock.test(std::memory_order_relaxed))
					utils::noop();
		}

		// *********************************************************************
		void unlock()
		{
			a_lock.clear();
		}
	};

    class spin_barrier 
    {
        const int initial;
        std::atomic<int> count;             
        std::atomic<int64_t> generation;    

    public:
        explicit spin_barrier(int count)
            : initial(count), count(count), generation(0)
        {
        }

        void arrive_and_wait() 
        {
            int64_t gen = generation.load(std::memory_order_acquire);

            if (count.fetch_sub(1, std::memory_order_acq_rel) == 1) 
            {
                count.store(initial, std::memory_order_release);
                generation.fetch_add(1, std::memory_order_acq_rel);
            }
            else 
            {
				while (generation.load(std::memory_order_acquire) == gen)
					utils::noop();
            }
        }

        void wait() 
		{
            int64_t gen = generation.load(std::memory_order_acquire);
            while (generation.load(std::memory_order_acquire) == gen)
				utils::noop();
        }
    };
}