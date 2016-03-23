#pragma once
//#define __CL_ENABLE_EXCEPTIONS
#include <exception>
#include <memory>
#include <string.h>
#include <assert.h>

#include <CL/cl.hpp>
#include <CL/cl.h>

#include "OpenCLError.h"
#include "DeviceWrapper.h"
#include "EventHandler.h"

cl_int clCall(cl_int code); 

namespace clex
{

/// <summary>
/// Represents the whole OpenCL framework.
/// </summary>
class OpenCL
{
public:
	static const int CPU_DEVICE = CL_DEVICE_TYPE_CPU;
	static const int GPU_DEVICE = CL_DEVICE_TYPE_GPU;
	static const int ANY_DEVICE	= CL_DEVICE_TYPE_ALL;

	/// Device type, can be either CPU_DEVICE or GPU_DEVICE.
	const int deviceType;
	
	std::vector<std::shared_ptr<DeviceWrapper>> devices; 
	
	std::shared_ptr<DeviceWrapper> mainDevice;

	std::shared_ptr<cl::Context> context;

	static double profileTimeMsec(cl::Event & profilingEvent)
	{
		cl_ulong start = profilingEvent.getProfilingInfo<CL_PROFILING_COMMAND_START>();
		cl_ulong end = profilingEvent.getProfilingInfo<CL_PROFILING_COMMAND_END>();
		double timeMsec = (double)(end - start) / 1000000.0;
		return timeMsec;
	}

	/// <summary>
	/// Lists all OpenCL-compatible devices.
	/// </summary>
	/// <param name="deviceType">Device type.</param>
	static std::string listDevices(int deviceType);
	
	/// <summary>
	/// Creates OpenCL environment for specified device type.
	/// </summary>
	/// <param name="deviceType">Device type.</param>
	OpenCL(int deviceType, int platformNum, int deviceNum, bool kernelProfiling);

	double profileKernel(
		const cl::Kernel& kernel,
		const cl::NDRange& offset,
		const cl::NDRange& global,
		const cl::NDRange& local);

protected:
	static void errorCallback(const char* info, const void *, ::size_t, void *) {
		std::cout << "OpenCL context error:" << info << std::endl;
	}
};

}