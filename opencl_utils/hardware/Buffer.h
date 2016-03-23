#pragma once
#include <CL/cl.hpp>

#include "../opencl_utils/common/Printable.h"
#include "OpenCl.h"

#include <sstream>

namespace clex
{

class Buffer : public cl::Buffer, public Printable
{
public:
	static ::size_t getTotalBytesAllocated() { return totalBytesAllocated; }
	static ::size_t getTotalBytesRequested() { return totalBytesRequested; }

	::size_t getBytesRequested() { return bytesRequested; }
	::size_t getBytesAllocated() { return bytesAllocated; }

	
	/// Buffer name.
	std::string name;

	cl_int getErr() const { return err; }

	Buffer(const OpenCL &cl, cl_mem_flags flags, ::size_t size, void *host_ptr = NULL, std::string name = "");

	//Buffer(const OpenCL &cl, cl_mem_flags flags, ::size_t size, const void *host_ptr = NULL, std::string name = "");

	//Buffer(const cl::Context& context, cl_mem_flags flags, ::size_t size, void *host_ptr = NULL);

	~Buffer();

	virtual std::string toString()
	{
		std::ostringstream oss;
		oss << "OpenCL " << name << ", requested[B] = " << bytesRequested << ", allocated[B]" << bytesAllocated;
		return oss.str();
	}

protected:
	
	static ::size_t totalBytesAllocated;
	static ::size_t totalBytesRequested;

	cl_int err;

	::size_t bytesRequested;
	::size_t bytesAllocated;
};

}
