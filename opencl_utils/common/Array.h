#pragma once
#include <vector>

template <class T>
class Array
{

public:
	size_t size() const { return width; }
	const std::vector<T>& getData() const { return v; } 
	std::vector<T>& getData() { return v; } 

	Array(int size) : width(size), height(size), v(size * size) {}
	Array() : width(0), height(0) {}
	Array(const Array<T>& ref) : width(ref.width), height(ref.height), v(ref.v) {}
	Array(Array&& rhs) : width(rhs.width), height(rhs.height), v(std::move(rhs.v)) {}

	T* operator[](const int row) { return v.data() + row * width; }
	const T* operator[](const int row) const { return v.data() + row * width; }

protected:
	int width;
	int height;
	std::vector<T> v; 
};