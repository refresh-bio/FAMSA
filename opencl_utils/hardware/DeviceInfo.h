#pragma once
#include <map>
#include <CL/cl.hpp>

#include "../opencl_utils/common/Printable.h"
#include "../opencl_utils/common/Named.h"
#include "../opencl_utils/common/mathex.h"

namespace clex
{

enum DeviceVendor {
	NVidia,
	AMD,
	Other
};

/// <summary>
/// Encapsulates all device information.
/// </summary>
class DeviceInfo : public Printable, public Named
{
public:
	std::string	deviceName;
	std::string	vendorName;
	DeviceVendor vendor;

	cl_device_type deviceType;

	cl_ulong globalMemSize;
	cl_ulong maxAllocSize;
	cl_ulong maxAllocSize_Corrected;

	cl_ulong localMemSize;
	cl_uint localMemType;
	cl_uint memAddrAlign;

	cl_ulong maxConstantBufferSize;
	cl_uint maxConstantArgsCount;
	cl_uint cacheType;
	cl_ulong cacheLineSize;
	cl_ulong cacheSize;

	cl_uint	maxDimensionsCount;
	cl_uint	computeUnitsCount;
	::size_t maxWorkgroupSize;
	cl_ulong maxParameterSize;

	cl_uint addressBits;

	std::string extensions;

	/// <summary>
	/// Empty constructor.
	/// </summary>
	DeviceInfo() {}

	/// <summary>
	/// Fills in all device information.
	/// </summary>
	/// <param name="device">Reference to device being queried.</param>
	DeviceInfo(const cl::Device &reference);

	/// <summary>
	/// Gets string representation of object.
	/// </summary>
	std::string toString();

	template <class T>
	T alignUp(T v) const { return mathex::ceilround(v, (T)memAddrAlign);  }

	template <class T>
	T alignDown(T v) const { return mathex::floorround(v, (T)memAddrAlign); }

	/// <summary>
	/// Gets object name.
	/// </summary>
	virtual std::string getName() { return deviceName; }

	/// <summary>
	/// Check if device supports specified extension.
	/// </summary>
	bool supportsExtension(std::string name) const { return extensions.find(name) != std::string::npos; }

private:
	std::map<cl_uint, std::string> cacheTypes;
	std::map<cl_uint, std::string> localMemTypes;
};

}