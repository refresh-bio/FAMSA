#ifndef _MEMORY_MONOTONIC_H
#define _MEMORY_MONOTONIC_H

#include <vector>
#include <memory>
#include <algorithm>
#include <utility>
#include <mutex>
#include <cstddef>

namespace refresh {
	class memory_monotonic_base
	{
	protected:
		size_t block_size;
		size_t alignment;
		size_t total_allocated;
		size_t total_requested;
		size_t no_allocs;
		size_t no_deallocs;

		std::vector<char*> blocks;
		std::vector<char*> freezed_blocks;
		char* cur_block;
		size_t in_block_pos;

		void _release_blocks()
		{
			for (auto p : blocks)
				free(p);

			blocks.clear();
			cur_block = nullptr;
			in_block_pos = block_size;
		}

		void _allocate_block(size_t size)
		{
			cur_block = (char*)malloc(size + alignment);
			total_allocated += size + alignment;

			blocks.push_back(cur_block);

			cur_block += alignment - (uint64_t)(cur_block) % alignment;
			in_block_pos = 0;
		}

		bool _deallocation_status()
		{
			return no_allocs == no_deallocs;
		}

		void* _allocate(size_t size)
		{
			if (in_block_pos + size > block_size)
				_allocate_block(std::max<size_t>(block_size, size));

			auto p = cur_block + in_block_pos;

			in_block_pos += (size + alignment - 1) / alignment * alignment;
			++no_allocs;

			total_requested += size;

			return p;
		}

		template<typename T>
		void _deallocate(T*& p)
		{
			if (!p)
				return;

			p = nullptr;
			++no_deallocs;
		}

		void _freeze()
		{
			freezed_blocks.insert(freezed_blocks.end(), blocks.begin(), blocks.end());
			blocks.clear();
			cur_block = nullptr;
			in_block_pos = block_size;
		}

		void _release_freezed()
		{
			for (auto& p : freezed_blocks)
				free(p);

			freezed_blocks.clear();
		}

		void _release()
		{
			_release_freezed();
			_release_blocks();
		}

	public:
		memory_monotonic_base(size_t _block_size, size_t _alignment) :
			block_size(_block_size),
			alignment(_alignment),
			total_allocated(0),
			total_requested(0),
			no_allocs(0),
			no_deallocs(0),
			cur_block(nullptr),
			in_block_pos(_block_size)
		{
		}

		memory_monotonic_base() = delete;
		memory_monotonic_base(const memory_monotonic_base &x) = delete;
		memory_monotonic_base(memory_monotonic_base &&x) = delete;
		memory_monotonic_base& operator=(const memory_monotonic_base &x) = delete;
		memory_monotonic_base& operator=(const memory_monotonic_base &&x) = delete;
		
		~memory_monotonic_base()
		{
			_release_freezed();
			_release_blocks();
		}
	};

	// ****
	class memory_monotonic_unsafe : public memory_monotonic_base
	{
	public:
		memory_monotonic_unsafe(size_t _block_size = (1ull << 20), size_t _alignment = 64)
			: memory_monotonic_base(_block_size, _alignment)
		{
		}

		memory_monotonic_unsafe() = delete;
		memory_monotonic_unsafe(const memory_monotonic_unsafe& x) = delete;
		memory_monotonic_unsafe(memory_monotonic_unsafe&& x) = delete;
		memory_monotonic_unsafe& operator=(const memory_monotonic_unsafe& x) = delete;
		memory_monotonic_unsafe& operator=(const memory_monotonic_unsafe&& x) = delete;

		~memory_monotonic_unsafe()
		{}

		bool deallocation_status()
		{
			return _deallocation_status();
		}

		void* allocate(size_t size)
		{
			return _allocate(size);
		}

		template<typename T>
		void deallocate(T*& p)
		{
			_deallocate(p);
		}

		void freeze()
		{
			_freeze();
		}

		void release()
		{
			_release();
		}

		void release_freezed()
		{
			_release_freezed();
		}
	};


	// ****
	class memory_monotonic_safe : public memory_monotonic_base
	{
		std::mutex mtx;

	public:
		memory_monotonic_safe(size_t _block_size = (1ull << 20), size_t _alignment = 64)
			: memory_monotonic_base(_block_size, _alignment)
		{
		}

		memory_monotonic_safe() = delete;
		memory_monotonic_safe(const memory_monotonic_safe& x) = delete;
		memory_monotonic_safe(memory_monotonic_safe&& x) = delete;
		memory_monotonic_safe& operator=(const memory_monotonic_safe& x) = delete;
		memory_monotonic_safe& operator=(const memory_monotonic_safe&& x) = delete;

		~memory_monotonic_safe()
		{			
		}

		bool deallocation_status()
		{
			std::lock_guard<std::mutex> lck(mtx);

			return _deallocation_status();
		}

		void* allocate(size_t size)
		{
			std::lock_guard<std::mutex> lck(mtx);

			return _allocate(size);
		}

		template<typename T>
		void deallocate(T*& p)
		{
			std::lock_guard<std::mutex> lck(mtx);

			_deallocate(p);
		}

		void freeze()
		{
			std::lock_guard<std::mutex> lck(mtx);

			_freeze();
		}

		void release()
		{
			std::lock_guard<std::mutex> lck(mtx);

			_release();
		}

		void release_freezed()
		{
			std::lock_guard<std::mutex> lck(mtx);

			_release_freezed();
		}
	};
}

#endif