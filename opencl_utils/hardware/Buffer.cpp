#include <iostream>
#include <sstream>

#include "Buffer.h"

#include "../opencl_utils/common/dbgnew.h"
#include "../opencl_utils/common/mathex.h"

::size_t clex::Buffer::totalBytesAllocated = 0;
::size_t clex::Buffer::totalBytesRequested = 0;

clex::Buffer::Buffer(const OpenCL &cl, cl_mem_flags flags, ::size_t size, void* host_ptr, std::string name) 
	: cl::Buffer(*cl.context, flags, size, host_ptr, &err)
{
	this->name = name;
	this->bytesRequested = size;
	this->bytesAllocated = mathex::ceilround(size, (size_t)cl.mainDevice->info->memAddrAlign);

	totalBytesRequested += this->bytesRequested;
	totalBytesAllocated += this->bytesAllocated;
	
	if (err != CL_SUCCESS) {
		std::cout << "Error while creating " << name << " "  ; 
		clCall(err);
	}
}

clex::Buffer::~Buffer()
{
	this->totalBytesAllocated -= bytesAllocated;
	this->totalBytesRequested -= bytesRequested;
	//std::cout << "destroying buffer " << name << "..." << std::endl;
}