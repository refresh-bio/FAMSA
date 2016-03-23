#include <memory>

template <typename T>
class PlacementAllocatorData {
public:
	T* currentPtr;
	T* startPtr;
	std::size_t allocated;

	PlacementAllocatorData(T* ptr) {
		currentPtr = ptr;
		startPtr = ptr;
		allocated = 0;
	}
};


template <typename T>
class PlacementAllocator : public std::allocator<T> {
public:
	typedef size_t     size_type;
	typedef ptrdiff_t  difference_type;
	typedef T*         pointer;
	typedef const T*   const_pointer;
	typedef T&         reference;
	typedef const T&   const_reference;
	typedef T          value_type;

	template<typename U>
	struct rebind
	{
		typedef PlacementAllocator <U> other; 
	};
	
	// some variables
	bool selfManage;
	std::shared_ptr<PlacementAllocatorData<T>> data;

	// additional stuff
	PlacementAllocator<T>() : std::allocator<T>() {

	}

	PlacementAllocator<T>(shared_ptr<PlacementAllocatorData<T>> data) : std::allocator<T>() {
		this->data = data;
		selfManage = false;
	}

	template <typename U>
	PlacementAllocator<T>(const PlacementAllocator<U>& ref) {
		selfManage = true; 
	}

	T* allocate(std::size_t num) {
		T* ptr;
		if (selfManage == true) {
			ptr = new T[num];
		} else {
			ptr = data->currentPtr;
			increase(num);
		}
		
		return ptr;
	}

	void deallocate(T* p, std::size_t num) {
		if (selfManage) {
			delete p;
		} else {
			decrease(num);
		}
	}

	void increase(int num) {
		data->currentPtr += num;
		data->allocated += num;
	}

	void decrease(int num) {
		data->currentPtr -= num;
		data->allocated -= num;
	}
};
