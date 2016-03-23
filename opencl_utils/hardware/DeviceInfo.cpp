#include <sstream>
#include "DeviceInfo.h"

#include "../opencl_utils/common/dbgnew.h"

#undef min
#undef max

using namespace std;

/// <summary>
/// Fills in all device information.
/// </summary>
/// <param name="device">Reference to device being queried.</param>
clex::DeviceInfo::DeviceInfo(const cl::Device &reference)
{
	localMemTypes[CL_LOCAL]			= "local";
	localMemTypes[CL_GLOBAL]		= "global";
	cacheTypes[CL_NONE]				= "none";
	cacheTypes[CL_READ_ONLY_CACHE]	= "read only";
	cacheTypes[CL_READ_WRITE_CACHE]	= "read / write";
	
	deviceName				= reference.getInfo<CL_DEVICE_NAME>();
	vendorName				= reference.getInfo<CL_DEVICE_VENDOR>();
	deviceType				= reference.getInfo<CL_DEVICE_TYPE>();

	if (vendorName.find("NVIDIA") != std::string::npos) {
		vendor = NVidia; 
	} else if (vendorName.find("AMD") != std::string::npos || vendorName.find("Advanced") != std::string::npos) {
		vendor = AMD;
	} else {
		vendor = Other;
	}

	globalMemSize			= reference.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>();
	maxAllocSize			= reference.getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>();
	maxAllocSize_Corrected  = (vendor == NVidia) ? 2 * maxAllocSize : maxAllocSize;
	memAddrAlign			= reference.getInfo<CL_DEVICE_MEM_BASE_ADDR_ALIGN>(); 

	cacheType				= reference.getInfo<CL_DEVICE_GLOBAL_MEM_CACHE_TYPE>();
	cacheSize				= reference.getInfo<CL_DEVICE_GLOBAL_MEM_CACHE_SIZE>();
	cacheLineSize			= reference.getInfo<CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE>();

	localMemSize			= reference.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>();
	localMemType			= reference.getInfo<CL_DEVICE_LOCAL_MEM_TYPE>();
	
	maxConstantBufferSize	= reference.getInfo<CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE>();
	maxConstantArgsCount	= reference.getInfo<CL_DEVICE_MAX_CONSTANT_ARGS>();

	computeUnitsCount		= reference.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
	maxDimensionsCount		= reference.getInfo<CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS>();
	maxWorkgroupSize		= reference.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
	maxParameterSize		= reference.getInfo<CL_DEVICE_MAX_PARAMETER_SIZE>();
	
	addressBits				= reference.getInfo<CL_DEVICE_ADDRESS_BITS>();

	extensions = reference.getInfo<CL_DEVICE_EXTENSIONS>();
}

/// <summary>
/// Gets string representation of object.
/// </summary>
std::string clex::DeviceInfo::toString()
{
	ostringstream out;
	out << 
		"DEVICE_NAME = " << deviceName << endl <<
		"VENDOR_NAME = " << vendorName << endl <<
		
		"GLOBAL_MEM_SIZE = " << globalMemSize << endl <<
		"MAX_ALLOCATION_SIZE = " << maxAllocSize << endl <<
		"CACHE TYPE = " << cacheTypes[cacheType] << endl <<
		"CACHE SIZE = " << cacheSize << endl <<
		"CACHELINE SIZE = " << cacheLineSize << endl <<
		
		"LOCAL_MEM_SIZE = " << localMemSize << endl <<
		"LOCAL_MEM_TYPE = " << localMemTypes[localMemType] << endl <<
		
		"MAX_CONSTANT_BUFFER_SIZE = " << maxConstantBufferSize << endl <<
		"MAX_CONSTANT_ARGS = " << maxConstantArgsCount << endl <<
		
		"MAX_COMPUTE_UNITS = " << computeUnitsCount << endl <<
		"MAX_WORK_ITEM_DIMENSIONS = " << maxDimensionsCount << endl <<
		"MAX_WORK_GROUP_SIZE = " << maxWorkgroupSize << endl <<
		"MAX_PARAMETER_SIZE = " << maxParameterSize << endl <<
		"ADDRESS_BITS = " << addressBits << endl <<

		"EXTENSIONS = " << extensions;

	return out.str();
}