#pragma once

#include <string>
#include <vector>
#include <CL/cl.hpp>

#include "../opencl_utils/common/Printable.h"
#include "../opencl_utils/common/Named.h"

namespace clex 
{

class PlatformInfo : public Named
{
public:
	std::string name;
	std::string version;
	std::string vendor;

	PlatformInfo(const cl::Platform &platform)
	{
		name = (std::string)platform.getInfo<CL_PLATFORM_NAME>();
		version = (std::string)platform.getInfo<CL_PLATFORM_VERSION>();
		vendor = (std::string)platform.getInfo<CL_PLATFORM_VENDOR>();
	}

	virtual std::string getName() { return name + " (" + version + ")"; }
};

}