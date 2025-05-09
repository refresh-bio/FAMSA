#pragma once

#include <cinttypes>
#include <thread>
#include <vector>
#include <chrono>
#include <atomic>
#include <exception>
#include <mutex>
#include <array>
#include <functional>
#include <string>
#include <future>
#include <cstdlib>
#include <cstring>
#include <bit>

#include <iostream>

#include "utils.h"

//#define REFRESH_ATP_DEBUG
//#define REFRESH_ATP_STATS

//#define REFRESH_ATP_RAW_ALLOC

//#define REFRESH_ATP_LAUNCH_PROFILE
//#define REFRESH_ATP_TASK_CALL_PROFILE

namespace refresh 
{
	// *********************************************************************
	class active_thread_pool_state
	{
		friend class active_thread_pool;
		friend class active_thread_pool_v2;

//		std::atomic_flag is_operating{};
		std::atomic<size_t> no_working{};
		std::exception_ptr exception_ptr{};

		std::mutex mtx;

		void inc()
		{
			no_working.fetch_add(1);
//			is_operating.test_and_set();
		}

		void dec()
		{
			no_working.fetch_sub(1);
/*			if (no_working.fetch_sub(1) == 1)
				is_operating.clear();

			is_operating.notify_all();*/
		}

		void set_exception(std::exception_ptr ptr)
		{
			std::lock_guard<std::mutex> lck(mtx);
			exception_ptr = ptr;
		}

	public:
		active_thread_pool_state() = default;
		
		void reserve(size_t n)	
		{
			// Do nothing here
		}
			
		std::exception_ptr get_exception()
		{
			std::lock_guard<std::mutex> lck(mtx);
			return exception_ptr;
		}

		bool was_exception()
		{
			std::lock_guard<std::mutex> lck(mtx);
			return exception_ptr != nullptr;
		}

/*		void wait()
		{
			// !!! FIXME: can work improperly
			// Można inicjalizować licznik na 1 i wtedy nigdy is_operating nie zostanie ustawiony na false zanim wszystkie wątki się odpalą i wejdzie do wait()
			is_operating.wait(true);
		}*/

		void heavy_wait()
		{
//			while (is_operating.test())
			while (no_working.load() != 0) [[likely]]
				;
		}

		void busy_wait()
		{
//			heavy_wait();
//			return;

//			while (is_operating.test())
			while (no_working.load() != 0) [[likely]]
				refresh::utils::noop();
		}

		bool is_ready()
		{
			return no_working.load() == 0;
//			return !is_operating.test();
		}
	};


	// *********************************************************************
	//template<unsigned MAX_POOL_SIZE_HARD_LIMIT = 256>
	class active_thread_pool_v2
	{
	public:
		using pool_state_t = active_thread_pool_state;

	private:
		static const size_t MAX_POOL_SIZE_HARD_LIMIT = 256;
		static const size_t NO_MASKS = (MAX_POOL_SIZE_HARD_LIMIT + 63) / 64;

		size_t busy_ticks;
		size_t max_pool_size;

		atomic<uint64_t> masks_waiting[NO_MASKS] = {};
		atomic<uint64_t> masks_sleeping[NO_MASKS] = {};
		
		atomic<size_t> tot_launches{};
		atomic<uint32_t> cur_rotation{};

		std::vector<std::thread> workers;

		//		utils::spin_mutex mtx;

		static const uint64_t thread_state_sleeping = 0;
		static const uint64_t thread_state_waiting = 1;
		static const uint64_t thread_state_code_loading = 2;
		static const uint64_t thread_state_ready_for_start = 3;
		static const uint64_t thread_state_working = 4;
		static const uint64_t thread_state_terminating = 5;

#ifdef REFRESH_ATP_STATS
		std::atomic<uint64_t> n_sleeping_use{};
		std::atomic<uint64_t> n_waiting_use{};
		std::atomic<uint64_t> n_waiting_to_sleeping{};
#endif


		// *********************************************************************
		struct worker_state_t
		{
			atomic<uint64_t> local_state{};
			pool_state_t* pool_ptr{ nullptr };
			std::function<void()> task{};

			static const size_t fields_size = sizeof(uint64_t) + sizeof(pool_state_t*) + sizeof(std::function<void()>);

			uint64_t aux[(128 - fields_size % 128) / 8 + 7]{};
			//			uint64_t aux[(8192 - fields_size % 8192) / 8] {};

			/*#ifdef REFRESH_ATP_STATS
						std::atomic<uint64_t> loc_sleeping_use{};
						std::atomic<uint64_t> loc_waiting_use{};
						std::atomic<uint64_t> loc_waiting_to_sleeping{};
			#endif*/

		public:
			// *********************************************************************
			worker_state_t() {};
		};

	public:
		void worker_body(atomic<uint64_t>& mask_waiting, atomic<uint64_t>& mask_sleeping, worker_state_t &worker_state, uint64_t worker_id, size_t busy_ticks)
		{
			auto& my_state = worker_state.local_state;
			uint64_t my_mask = 1ull << (worker_id % 64);

			while (true)
			{
				auto current_state = my_state.load();

#ifdef REFRESH_ATP_DEBUG
				std::cerr << "my state: " + std::to_string(internal_id) + " - " + std::to_string(current_state) + "\n";
#endif

				switch (current_state)
				{
				case thread_state_sleeping:
					my_state.wait(thread_state_sleeping);
					break;
				case thread_state_terminating:
					return;
				case thread_state_waiting:
				{
					for (size_t i = 0; i < busy_ticks; ++i)
					{
//						if (my_state.load(memory_order_relaxed) != thread_state_waiting)
						if (my_state.load() != thread_state_waiting) [[unlikely]]
							break;
						utils::noop();
					}

					{
						auto tmp = thread_state_waiting;
						my_state.compare_exchange_strong(tmp, thread_state_sleeping);
					}
					break;
				}
				case thread_state_code_loading:
				{
					auto tmp = thread_state_code_loading;
					my_state.compare_exchange_strong(tmp, thread_state_ready_for_start);

//					while (my_state.load(memory_order_relaxed) != thread_state_working)
					while (my_state.load() != thread_state_working)
					{
						utils::noop();
					}
//					break;
				}
					[[fallthrough]];
				case thread_state_working:
#ifdef REFRESH_ATP_TASK_CALL_PROFILE
					{
					high_resolution_clock::time_point t1, t2, t3, t4, t5;
					t1 = high_resolution_clock::now();
#endif
					try
					{
#ifdef REFRESH_ATP_TASK_CALL_PROFILE
						t2 = high_resolution_clock::now();
#endif
						worker_state.task();
#ifdef REFRESH_ATP_TASK_CALL_PROFILE
						t3 = high_resolution_clock::now();
#endif

						auto curr_pool_ptr = worker_state.pool_ptr;

						my_state.store(thread_state_waiting);
						mask_waiting |= my_mask;

						if (curr_pool_ptr)
							curr_pool_ptr->dec();

#ifdef REFRESH_ATP_TASK_CALL_PROFILE
						t4 = high_resolution_clock::now();
#endif
					}
					catch (const std::exception& e)
					{
						if (worker_state.pool_ptr)
						{
							std::cerr << "Exception: " + std::string(e.what()) + "\n";
							fflush(stderr);

							auto curr_pool_ptr = worker_state.pool_ptr;

							my_state.store(thread_state_waiting);
							mask_waiting |= my_mask;

							curr_pool_ptr->set_exception(std::current_exception());
							curr_pool_ptr->dec();

							return;
						}
					}

#ifdef REFRESH_ATP_TASK_CALL_PROFILE
					t5 = high_resolution_clock::now();

					auto d1 = 1000000 * std::chrono::duration<double>(t2 - t1).count();
					auto d2 = 1000000 * std::chrono::duration<double>(t3 - t2).count();
					auto d3 = 1000000 * std::chrono::duration<double>(t4 - t3).count();
					auto d4 = 1000000 * std::chrono::duration<double>(t5 - t4).count();

					cout << to_string(d1) + "us  " + to_string(d2) + "us  " + to_string(d3) + "us  " + to_string(d4) + "us\n";
				}
#endif

					break;
				default:
					exit(1);
				}
			}
		}

#ifndef REFRESH_ATP_RAW_ALLOC
		std::vector<worker_state_t> worker_states{ MAX_POOL_SIZE_HARD_LIMIT };
#else	
		uint64_t* worker_raw;
		worker_state_t* worker_states;
#endif

		// *********************************************************************
		void empty_loop(size_t n)
		{
			volatile size_t aux = 0;

			for (size_t i = 0; i < n; ++i)
			{
				aux = aux + 1;
				refresh::utils::noop();
			}
		}

		// *********************************************************************
		void adjust_busy_time(const std::chrono::duration<double> busy_time)
		{
			busy_ticks = 8;

			while (busy_ticks < (1 << 30))
			{
				auto t1 = std::chrono::high_resolution_clock::now();
				empty_loop(busy_ticks);
				auto t2 = std::chrono::high_resolution_clock::now();

				const std::chrono::duration<double> diff = t2 - t1;

				if (diff >= busy_time)
					break;

				if (busy_ticks & (busy_ticks - 1))
				{
					busy_ticks &= busy_ticks - 1;
					busy_ticks *= 2;
				}
				else
					busy_ticks += busy_ticks / 2;
			}
		}

		// *********************************************************************
		void init_workers(const size_t init_pool_size)
		{
			workers.resize(max_pool_size);

#ifndef REFRESH_ATP_RAW_ALLOC
//			worker_states.resize(max_pool_size);
#else
			worker_raw = (uint64_t*)std::malloc(sizeof(worker_state_t) * max_pool_size + 64);

			uint64_t* ptr;

			for (ptr = worker_raw; ((uint64_t)ptr) % 64 != 0; ptr++)
				;

			worker_states = (worker_state_t*)ptr;

			for (size_t i = 0; i < max_pool_size; ++i)
				new (worker_states + i) worker_state_t();
#endif

			for (size_t i = 0; i < max_pool_size; ++i)
			{
				if (i < init_pool_size)
				{
					masks_waiting[i / 64] |= 1ull << (i % 64);
					worker_states[i].local_state = thread_state_waiting;
				}
				else
				{
					masks_sleeping[i / 64] |= 1ull << (i % 64);
					worker_states[i].local_state = thread_state_sleeping;
				}

				worker_states[i].pool_ptr = nullptr;
				workers[i] = thread([&, i] {worker_body(masks_waiting[i / 64], masks_sleeping[i / 64], worker_states[i], i, busy_ticks); });
			}
		}

	public:
		// *********************************************************************
		active_thread_pool_v2(size_t init_pool_size, size_t max_pool_size, std::chrono::duration<double> busy_time = std::chrono::milliseconds(2)) :
			max_pool_size(max_pool_size)
		{
			adjust_busy_time(busy_time);

#ifdef REFRESH_ATP_DEBUG
			std::cout << "atp busy_ticks: " << busy_ticks << std::endl;
#endif

			init_workers(init_pool_size);
		}

		// *********************************************************************
		active_thread_pool_v2(size_t init_pool_size, size_t max_pool_size, size_t busy_ticks) :
			max_pool_size(max_pool_size),
			busy_ticks(busy_ticks)
		{
#ifdef REFRESH_ATP_DEBUG
			std::cout << "atp busy_ticks: " << busy_ticks << std::endl;
#endif

			init_workers(init_pool_size);
		}

		// *********************************************************************
		~active_thread_pool_v2()
		{
#ifdef REFRESH_ATP_DEBUG
			cout << "No. ATP launches: " << tot_launches.load() << endl;
#endif
			cancel();
		}

		// *********************************************************************
		bool is_active()	const
		{
			return true;
		}

		// *********************************************************************
		void launch(std::function<void()>&& fun, active_thread_pool_v2::pool_state_t* pool_state)
			//		void launch(const std::function<void()>& fun, active_thread_pool::pool_state_t* pool_state)
		{
			//			mtx.lock();

#ifdef REFRESH_ATP_LAUNCH_PROFILE
			auto t1 = high_resolution_clock::now();
#endif

			tot_launches++;

			auto my_rotation = cur_rotation.fetch_add(1) & 63;

			int64_t worker_id = -1;

			while (true)
			{
				// Try to find waiting worker with the same job type
				for (int i = 0; i < (max_pool_size + 63) / 64 && worker_id < 0; ++i)
				{
					uint64_t m = masks_waiting[i];
					uint64_t low_masker = ~0ull;

					while (m & low_masker)
					{
						int id = countr_zero(m & low_masker);

						uint64_t new_m = m & ~(1ull << id);

						if (masks_waiting[i].compare_exchange_strong(m, new_m)) [[likely]]
						{
							worker_id = i * 64 + id;

							if (strcmp(worker_states[worker_id].task.target_type().name(), fun.target_type().name()) != 0)
							{
								masks_waiting[i] |= 1ull << id;
								worker_id = -1;

								if (id == 63)
									low_masker = 0;
								else
									low_masker = ~0ull << (id + 1);
								continue;
							}

							uint64_t tmp = thread_state_waiting;
							if (!worker_states[worker_id].local_state.compare_exchange_strong(tmp, thread_state_code_loading))
							{
								masks_sleeping[i] |= 1ull << id;
								worker_id = -1;
								continue;
							}

							break;
						}
					}
				}

				if (worker_id >= 0)
					break;

				// Try to find any waiting worker
				for (int i = 0; i < (max_pool_size + 63) / 64 && worker_id < 0; ++i)
				{
					uint64_t m = masks_waiting[i];

					while (m)
					{
						int id = countr_zero(m);

						uint64_t new_m = m & ~(1ull << id);

						if (masks_waiting[i].compare_exchange_strong(m, new_m)) [[likely]]
						{
							worker_id = i * 64 + id;

							uint64_t tmp = thread_state_waiting;
							if (worker_states[worker_id].local_state.compare_exchange_strong(tmp, thread_state_code_loading))
								break;
							masks_sleeping[i] |= 1ull << id;
							worker_id = -1;
						}
					}
				}

				if (worker_id >= 0)
					break;

				// Try to find sleeping worker
				for (int i = 0; i < (max_pool_size + 63) / 64 && worker_id < 0; ++i)
				{
					uint64_t m = masks_sleeping[i];

					while (m)
					{
						int id = countr_zero(m);
						uint64_t new_m = m & ~(1ull << id);

						if (masks_sleeping[i].compare_exchange_strong(m, new_m)) [[likely]]
						{
							worker_id = i * 64 + id;
							worker_states[worker_id].local_state = thread_state_code_loading;
							break;
						}
					}
				}

				if (worker_id >= 0)
					break;
			}
			
#ifdef REFRESH_ATP_LAUNCH_PROFILE
			auto t2 = high_resolution_clock::now();
#endif

			pool_state->inc();
			auto& worker_state = worker_states[worker_id];
			worker_state.pool_ptr = pool_state;

#ifdef REFRESH_ATP_LAUNCH_PROFILE
			auto t3 = high_resolution_clock::now();
#endif
			worker_state.task = std::forward<std::function<void()>>(fun);
			//worker_state.task.swap(fun);

#ifdef REFRESH_ATP_LAUNCH_PROFILE
			auto t4 = high_resolution_clock::now();
#endif
			bool need_to_notify = true;
			auto tmp_state = thread_state_ready_for_start;
			if (worker_state.local_state.compare_exchange_strong(tmp_state, thread_state_working))
				need_to_notify = false;
			else
				worker_state.local_state = thread_state_working;

#ifdef REFRESH_ATP_LAUNCH_PROFILE
			auto t5 = high_resolution_clock::now();
#endif
			if(need_to_notify)
				worker_state.local_state.notify_one();

#ifdef REFRESH_ATP_LAUNCH_PROFILE
			auto t6 = high_resolution_clock::now();

			auto d1 = 1000000 * std::chrono::duration<double>(t2 - t1).count();
			auto d2 = 1000000 * std::chrono::duration<double>(t3 - t2).count();
			auto d3 = 1000000 * std::chrono::duration<double>(t4 - t3).count();
			auto d4 = 1000000 * std::chrono::duration<double>(t5 - t4).count();
			auto d5 = 1000000 * std::chrono::duration<double>(t6 - t5).count();

			cout << to_string(d1) + "us  " + to_string(d2) + "us  " + to_string(d3) + "us  " + to_string(d4) + "us   " + to_string(d5) + "us\n";
#endif

			//			mtx.unlock();
		}

		// *********************************************************************
		void cancel()
		{
			//mtx.lock();

			// Set state to terminating for all workers
#ifndef REFRESH_ATP_RAW_ALLOC
			for (size_t i = 0; i < max_pool_size; ++i)
#else
			for (size_t i = 0; i < max_pool_size; ++i)
#endif
			{
				worker_states[i].local_state = thread_state_terminating;
				worker_states[i].local_state.notify_one();

				/*#ifdef REFRESH_ATP_STATS
								n_waiting_to_sleeping += workers[i].loc_waiting_to_sleeping;
				#endif*/
			}

			for (auto& t : workers)
				t.join();
			workers.clear();

#ifndef REFRESH_ATP_RAW_ALLOC
//			worker_states.clear();
#else
			for (size_t i = 0; i < max_pool_size; ++i)
				(worker_states + i)->~worker_state_t();
			free(worker_raw);
			worker_raw = nullptr;
			worker_states = nullptr;
#endif

#ifdef REFRESH_ATP_STATS
			std::cerr << "No. waiting use        : " << n_waiting_use << std::endl;
			std::cerr << "No. sleeping use       : " << n_sleeping_use << std::endl;
			std::cerr << "No. waiting to sleeping: " << n_waiting_to_sleeping << std::endl;
#endif

			//mtx.unlock();
		}
	};


	// *********************************************************************
	class active_thread_pool
	{
	public:
		using pool_state_t = active_thread_pool_state;

	private:
		size_t busy_ticks;
		size_t max_pool_size;
		atomic<size_t> prev_worker_id{};

		atomic<size_t> tot_launches{};

//		utils::spin_mutex mtx;

		static const uint64_t thread_state_sleeping = 0;
		static const uint64_t thread_state_waiting = 1;
		static const uint64_t thread_state_code_loading = 2;
		static const uint64_t thread_state_working = 3;
		static const uint64_t thread_state_terminating = 4;

#ifdef REFRESH_ATP_STATS
		std::atomic<uint64_t> n_sleeping_use{};
		std::atomic<uint64_t> n_waiting_use{};
		std::atomic<uint64_t> n_waiting_to_sleeping{};
#endif

		// *********************************************************************
		struct worker_t
		{
			uint64_t local_state;
			pool_state_t* pool_ptr;
			std::function<void()> task{};
			std::thread internal_thread;

			static const size_t fields_size = sizeof(uint64_t) + sizeof(pool_state_t*) + sizeof(std::thread) + sizeof(std::function<void()>);

			uint64_t aux[(128 - fields_size % 128) / 8 + 7];
//			uint64_t aux[(8192 - fields_size % 8192) / 8];

/*#ifdef REFRESH_ATP_STATS
			std::atomic<uint64_t> loc_sleeping_use{};
			std::atomic<uint64_t> loc_waiting_use{};
			std::atomic<uint64_t> loc_waiting_to_sleeping{};
#endif*/

		private:
			// *********************************************************************
			void worker_body(uint64_t *global_state_ptr, size_t busy_ticks)
			{
				while (true)
				{
					auto current_state = std::atomic_ref<uint64_t>(local_state).load();

#ifdef REFRESH_ATP_DEBUG
					std::cerr << "my state: " + std::to_string(internal_id) + " - " + std::to_string(current_state) + "\n";
#endif

					switch (current_state)
					{
					case thread_state_sleeping:
						std::atomic_ref<uint64_t>(local_state).wait(thread_state_sleeping);
						break;
					case thread_state_terminating:
						return;
					case thread_state_waiting:
					{
						for (size_t i = 0; i < busy_ticks; ++i)
						{
							if (std::atomic_ref<uint64_t>(local_state).load() != thread_state_waiting)
								break;
							utils::noop();
						}

						{
							auto tmp = thread_state_waiting;
							std::atomic_ref<uint64_t>(local_state).compare_exchange_strong(tmp, thread_state_sleeping);
							tmp = thread_state_waiting;
							std::atomic_ref<uint64_t>(*global_state_ptr).compare_exchange_strong(tmp, thread_state_sleeping);

#ifdef REFRESH_ATP_STATS
#else
#endif
						}
						break;
					}
					case thread_state_code_loading:
						while (std::atomic_ref<uint64_t>(local_state).load() == thread_state_code_loading)
							utils::noop();
//						break;
						[[fallthrough]];
					case thread_state_working:
						try
						{
							task();

							auto curr_pool_ptr = pool_ptr;

							std::atomic_ref<uint64_t>(local_state).store(thread_state_waiting);
							std::atomic_ref<uint64_t>(*global_state_ptr).store(thread_state_waiting);

							if (curr_pool_ptr)
								curr_pool_ptr->dec();
						}
						catch (const std::exception &e)
						{
							if (pool_ptr)
							{
								std::cerr << "Exception: " + std::string(e.what()) + "\n";
								fflush(stderr);

								auto curr_pool_ptr = pool_ptr;

								std::atomic_ref<uint64_t>(local_state).store(thread_state_waiting);
								std::atomic_ref<uint64_t>(*global_state_ptr).store(thread_state_waiting);

								curr_pool_ptr->set_exception(std::current_exception());
								curr_pool_ptr->dec();

								return;
							}
						}

						break;
					}
				}
			}

		public:
			// *********************************************************************
			worker_t() = delete;

			worker_t(size_t internal_id, uint64_t* global_state_ptr, bool is_waiting, size_t busy_ticks) :
				local_state(is_waiting ? thread_state_waiting : thread_state_sleeping),
				pool_ptr(nullptr)
			{
				*global_state_ptr = local_state;
				internal_thread = std::thread([&, global_state_ptr, busy_ticks] {worker_body(global_state_ptr, busy_ticks); });

//				cout << sizeof(worker_t) << endl;
			}
		};

		// *********************************************************************
		alignas(64) std::vector<uint64_t> global_states;

#ifndef REFRESH_ATP_RAW_ALLOC
		std::vector<worker_t> workers;
#else	
		uint64_t* worker_raw;
		worker_t* workers;
#endif

		// *********************************************************************
		void empty_loop(size_t n)
		{
			volatile size_t aux = 0;

			for (size_t i = 0; i < n; ++i)
			{
				aux = aux + 1;
				refresh::utils::noop();
			}
		}

		// *********************************************************************
		void adjust_busy_time(const std::chrono::duration<double> busy_time)
		{
			busy_ticks = 8;

			while (busy_ticks < (1 << 30))
			{
				auto t1 = std::chrono::high_resolution_clock::now();
				empty_loop(busy_ticks);
				auto t2 = std::chrono::high_resolution_clock::now();

				const std::chrono::duration<double> diff = t2 - t1;

				if (diff >= busy_time)
					break;

				if (busy_ticks & (busy_ticks - 1))
				{
					busy_ticks &= busy_ticks - 1;
					busy_ticks *= 2;
				}
				else
					busy_ticks += busy_ticks / 2;
			}
		}
		
		// *********************************************************************
		void init_workers(const size_t init_pool_size)
		{
			global_states.reserve(max_pool_size);

#ifndef REFRESH_ATP_RAW_ALLOC
			workers.reserve(max_pool_size);

			for (size_t i = 0; i < max_pool_size; ++i)
				workers.emplace_back(i, &(global_states[i]), i < init_pool_size, busy_ticks);
#else
			worker_raw = (uint64_t*) std::malloc(sizeof(worker_t) * max_pool_size + 64);

			uint64_t *ptr;

			for (ptr = worker_raw; ((uint64_t)ptr) % 64 != 0; ptr++)
				;

			workers = (worker_t*)ptr;

			for (size_t i = 0; i < max_pool_size; ++i)
				new (workers+i) worker_t(i, &(global_states[i]), i < init_pool_size, busy_ticks);
#endif
		}

	public:
		// *********************************************************************
		active_thread_pool(size_t init_pool_size, size_t max_pool_size, std::chrono::duration<double> busy_time = std::chrono::milliseconds(2)) :
			max_pool_size(max_pool_size)
		{
			adjust_busy_time(busy_time);

#ifdef REFRESH_ATP_DEBUG
			std::cout << "atp busy_ticks: " << busy_ticks << std::endl;
#endif

			init_workers(init_pool_size);
		}

		// *********************************************************************
		active_thread_pool(size_t init_pool_size, size_t max_pool_size, size_t busy_ticks) :
			max_pool_size(max_pool_size),
			busy_ticks(busy_ticks)
		{
#ifdef REFRESH_ATP_DEBUG
			std::cout << "atp busy_ticks: " << busy_ticks << std::endl;
#endif

			init_workers(init_pool_size);
		}

		// *********************************************************************
		~active_thread_pool()
		{
			cout << "No. ATP launches: " << tot_launches.load() << endl;

			cancel();
		}

		// *********************************************************************
		bool is_active()	const
		{
			return true;
		}

		// *********************************************************************
		void launch(std::function<void()>&& fun, active_thread_pool::pool_state_t* pool_state)
//		void launch(const std::function<void()>& fun, active_thread_pool::pool_state_t* pool_state)
		{
//			mtx.lock();

			tot_launches++;

			size_t n_workers = max_pool_size;
			size_t worker_id = (prev_worker_id + 1) % n_workers;
			size_t i = 0;

			// Find any waiting thread
			for (i = 0; i < n_workers; ++i)
			{
				uint64_t curr_state = thread_state_waiting;
				
				if (std::atomic_ref<uint64_t>(global_states[worker_id]).compare_exchange_strong(curr_state, thread_state_code_loading))
				{
#ifdef REFRESH_ATP_STATS
					n_waiting_use++;
#endif
					break;
				}

				worker_id = (worker_id + 1 == n_workers) ? 0 : (worker_id + 1);
			}

			// If necessary look for sleeping thread
			if (i == n_workers)
			{
				for (i = 0; i < n_workers; ++i)
				{
					uint64_t curr_state = thread_state_sleeping;
					
					if (std::atomic_ref<uint64_t>(global_states[worker_id]).compare_exchange_strong(curr_state, thread_state_code_loading))
					{
#ifdef REFRESH_ATP_STATS
						n_sleeping_use++;
#endif
						break;
					}

					worker_id = (worker_id + 1 == n_workers) ? 0 : (worker_id + 1);
				}
			}

#ifdef REFRESH_ATP_DEBUG
			std::cerr << "worker_id: " + std::to_string(worker_id) + "\n";
#endif

			pool_state->inc();
			workers[worker_id].pool_ptr = pool_state;
//			workers[worker_id].task = std::forward<std::function<void()>>(fun);
			workers[worker_id].task.swap(fun);

			std::atomic_ref<uint64_t>(global_states[worker_id]).store(thread_state_working);
			std::atomic_ref<uint64_t>(workers[worker_id].local_state).store(thread_state_working);
			std::atomic_ref<uint64_t>(workers[worker_id].local_state).notify_one();

			prev_worker_id = worker_id;

//			mtx.unlock();
		}

		// *********************************************************************
		void cancel()
		{
			//mtx.lock();

			// Set state to terminating for all workers
#ifndef REFRESH_ATP_RAW_ALLOC
			for (size_t i = 0; i < workers.size(); ++i)
#else
			for (size_t i = 0; i < max_pool_size; ++i)
#endif
			{
				std::atomic_ref<uint64_t> thread_state(workers[i].local_state);
				thread_state.store(thread_state_terminating);
				thread_state.notify_one();
				workers[i].internal_thread.join();

/*#ifdef REFRESH_ATP_STATS
				n_waiting_to_sleeping += workers[i].loc_waiting_to_sleeping;
#endif*/
			}

#ifndef REFRESH_ATP_RAW_ALLOC
			workers.clear();
#else
			for (size_t i = 0; i < max_pool_size; ++i)
				(workers + i)->~worker_t();
			free(worker_raw);
			worker_raw = nullptr;
			workers = nullptr;
#endif

#ifdef REFRESH_ATP_STATS
			std::cerr << "No. waiting use        : " << n_waiting_use << std::endl;
			std::cerr << "No. sleeping use       : " << n_sleeping_use << std::endl;
			std::cerr << "No. waiting to sleeping: " << n_waiting_to_sleeping << std::endl;
#endif

			//mtx.unlock();
		}
	};

	// *********************************************************************
	class async_pool_state
	{
		friend class async_pool;

		std::exception_ptr exception_ptr{};

		std::mutex mtx;

		std::vector<std::future<void>> workers;

		void add_future(std::future<void>&& fut)
		{
			std::lock_guard<std::mutex> lck(mtx);
			workers.emplace_back(std::move(fut));
		}

		void set_exception(std::exception_ptr ptr)
		{
			std::lock_guard<std::mutex> lck(mtx);
			exception_ptr = ptr;
		}

	public:
		async_pool_state() = default;

		void reserve(size_t n)
		{
			std::lock_guard<std::mutex> lck(mtx);
			workers.reserve(n);						
		}

		std::exception_ptr get_exception()
		{
			std::lock_guard<std::mutex> lck(mtx);
			return exception_ptr;
		}

		bool was_exception()
		{
			std::lock_guard<std::mutex> lck(mtx);
			return exception_ptr != nullptr;
		}

		void wait()
		{
			std::lock_guard<std::mutex> lck(mtx);
			for (auto& f : workers)
			{
				try {
					f.wait();
				}
				catch (...)
				{
					set_exception(std::current_exception());
				}
			}

			workers.clear();
		}

		void heavy_wait()
		{
			wait();
		}

		void busy_wait()
		{
			wait();
		}

		bool is_ready()
		{
			std::lock_guard<std::mutex> lck(mtx);
			return workers.empty();
		}
	};

	// *********************************************************************
	class async_pool
	{
	public:
		using pool_state_t = async_pool_state;

	private:

	public:
		// *********************************************************************
		// Parameters just to give same interface at active_thread_pool
		async_pool(size_t init_pool_size = 0, size_t max_pool_size = 0, std::chrono::duration<double> busy_time = std::chrono::milliseconds(0))
		{
		}

		// *********************************************************************
		~async_pool()
		{
		}

		// *********************************************************************
		bool is_active()	const
		{
			return true;
		}

		// *********************************************************************
		void launch(std::function<void()>&& fun, async_pool::pool_state_t* pool_state)
		{
			pool_state->add_future(std::async(std::launch::async, fun));
		}

		// *********************************************************************
		void cancel()
		{
			// Do nothing
		}
	};
}
