#pragma once
#include <memory>
#include <exception>
#include <sstream>
#include <iostream>
#include "CL/cl.hpp"

#include "OpenCl.h"
#include "Buffer.h"

namespace clex
{

class Kernel : public cl::Kernel
{
public:
	std::string name;

	::size_t maxWorkgroupSize;
	::size_t preferredWorkgroupSizeMultiple;
	cl_ulong localMemSize;
	cl_ulong privateMemSize;

	cl_kernel getObject() { return this->object_; }

	/// <summary>
	/// </summary>
	/// <param name=openCl></param>
	/// <param name=program></param>
	/// <param name=name></param>
	/// <param name=err></param>
	/// <returns></returns>
	Kernel(std::shared_ptr<clex::OpenCL> cl, const cl::Program& program, const char* name)
		: cl::Kernel(program, name, &err), cl(cl), name(name)
	{
		updateProperties();
	}

	/// <summary>
	/// </summary>
	/// <param name=ref></param>
	/// <returns></returns>
	Kernel(const Kernel& ref) : cl::Kernel(ref), cl(ref.cl), name(ref.name)
	{
		updateProperties();
	}

	/// <summary>
	///
	/// </summary>
	template <typename T>
	cl_int setArg(cl_uint index, T value)
	{
		err = cl::Kernel::setArg(index, value);
		updateProperties();
		return clCall(err);
	}

	/// <summary>
	///
	/// </summary>
	cl_int setArg(cl_uint index, Buffer& ownBuffer)
	{
		err = cl::Kernel::setArg(index, dynamic_cast<cl::Buffer&>(ownBuffer));
		updateProperties();
		return clCall(err);
	}

	/// <summary>
	///
	/// </summary>
	cl_int setArg(cl_uint index, ::size_t size, void* argPtr)
	{
		err = cl::Kernel::setArg(index, size, argPtr);
		updateProperties();
		return clCall(err);
	}

	void updateProperties()
	{
		maxWorkgroupSize = getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(cl->mainDevice->device);
		preferredWorkgroupSizeMultiple = getWorkGroupInfo<CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE>(cl->mainDevice->device);
		localMemSize = getWorkGroupInfo<CL_KERNEL_LOCAL_MEM_SIZE>(cl->mainDevice->device);
		privateMemSize = getWorkGroupInfo<CL_KERNEL_PRIVATE_MEM_SIZE>(cl->mainDevice->device);
	}

protected:
	cl_int err;

	std::shared_ptr<clex::OpenCL> cl;
};

}