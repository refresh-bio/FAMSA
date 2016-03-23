#include <iostream>
#include <fstream>
#include <sstream>
#include <exception>
#include "OpenCl.h"
#include "OpenCLError.h"
#include "PlatformInfo.h"
#include "DeviceInfo.h"

#include "../opencl_utils/common/dbgnew.h"
#include "../opencl_utils/common/Timer.h"
#include "../opencl_utils/common/MemoryTools.h"

using namespace cl;

cl_int clCall(cl_int code)
{
	std::string msg = clex::OpenCLError::name(code);

	if(code != CL_SUCCESS) {
		std::cout << "Error: " << msg << std::endl;
	}

	return code;
}


std::string clex::OpenCL::listDevices(int deviceType)
{
	std::ostringstream oss;
	std::vector<Platform> platforms;
	Platform::get(&platforms);

	for (int i = 0; i < platforms.size(); i++) {
		oss << "Platform " << i << ": " << PlatformInfo(platforms[i]).getName() << std::endl;
		std::vector<Device> devices;
		platforms[i].getDevices(deviceType, &devices);

		for (int j = 0; j < devices.size(); j++) {
			oss << "\tDevice " << j << ": " << DeviceInfo(devices[j]).getName() << std::endl;
		}
	}

	return oss.str();
}

clex::OpenCL::OpenCL(int deviceType, int platformNum, int deviceNum, bool kernelProfiling) :
	deviceType(deviceType)
{
	LOG_MEM("Before OpenCL initialisation");
	int code;
	std::vector<Platform> platforms;
	clCall(Platform::get(&platforms));
	LOG_MEM("Platforms get");

	std::vector<Device> tempDevices;
	clCall(platforms[platformNum].getDevices(deviceType, &tempDevices));
	LOG_MEM("Devices get");
	
	// fixme:
	std::vector<Device> finalDevices(1, tempDevices[deviceNum]);
	
	context = std::shared_ptr<Context>(new cl::Context(finalDevices, NULL, NULL, NULL, &code));
	clCall(code);
	LOG_MEM("Context created");

	this->devices.push_back(std::shared_ptr<DeviceWrapper>(new DeviceWrapper(finalDevices[0], *context, kernelProfiling)));
	this->mainDevice = devices[0];
	LOG_MEM("Device wrapper created");
}


double clex::OpenCL::profileKernel(
		const cl::Kernel& kernel,
		const cl::NDRange& offset,
		const cl::NDRange& global,
		const cl::NDRange& local)
{
	cl::Event finishEvent;

	clCall(mainDevice->mainQueue->enqueueNDRangeKernel(kernel, offset, global, local, NULL, &finishEvent));

	cl_ulong timeStart, timeEnd;
	double timeMsec;
	finishEvent.wait();
	finishEvent.getProfilingInfo(CL_PROFILING_COMMAND_START, &timeStart);
	finishEvent.getProfilingInfo(CL_PROFILING_COMMAND_END, &timeEnd);
	timeMsec = (timeEnd - timeStart) / 1000000.0;

	return timeMsec;
}
