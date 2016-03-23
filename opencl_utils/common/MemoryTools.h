#pragma once
#include <cstddef>
#define LOG_MEM(x)

#if defined _WIN32
	#define WINDOWS
	#include "WindowsMemoryTools.h"
#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
	#define LINUX
	#include "UnixMemoryTools.h"
#endif


class MemoryTools
{
public:

#if defined WINDOWS
	static ::size_t processCurrentVirtual()	{ return WindowsMemoryTools::processCurrentVirtual(); }
	static ::size_t processPeakVirtual()		{ return WindowsMemoryTools::processPeakVirtual(); }
	static ::size_t systemTotalPhysical()		{ return WindowsMemoryTools::systemTotalPhysical(); }
	static ::size_t systemAvailablePhysical()	{ return WindowsMemoryTools::systemAvailablePhysical(); }
#elif defined LINUX
	static ::size_t processCurrentVirtual()	{ return UnixMemoryTools::processCurrentVirtual(); }
	static ::size_t processPeakVirtual()		{ return UnixMemoryTools::processPeakVirtual(); }
	static ::size_t systemTotalPhysical()		{ return UnixMemoryTools::systemTotalPhysical(); }
	static ::size_t systemAvailablePhysical()	{ return UnixMemoryTools::systemAvailablePhysical(); }
#else
	#error "Unable to compile MemoryTools for an unknown OS."
#endif
};