#pragma once
#include <vector>
#include <memory>
#include <CL/cl.hpp>

#include "DeviceInfo.h"

namespace clex
{

class DeviceWrapper
{
public:
	cl::Device device;

	std::vector<std::shared_ptr<cl::CommandQueue>> queues;

	std::shared_ptr<cl::CommandQueue> mainQueue;

	std::shared_ptr<DeviceInfo> info;

	DeviceWrapper(const cl::Device& ref, const cl::Context& context, bool useKernelProfiling);	

};

}