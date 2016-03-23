#include "WindowsMemoryTools.h"
#include "MemoryTools.h"

#ifdef WINDOWS
#include <Windows.h>
#include <Psapi.h>

size_t WindowsMemoryTools::processCurrentVirtual()
{
	PROCESS_MEMORY_COUNTERS counters;
	if (GetProcessMemoryInfo(GetCurrentProcess(), &counters, sizeof (counters))) {
		return counters.PagefileUsage;
	}
	
	return 0;
}

size_t WindowsMemoryTools::processPeakVirtual()
{
	PROCESS_MEMORY_COUNTERS counters;
	if (GetProcessMemoryInfo (GetCurrentProcess(), &counters, sizeof (counters))) {
		return counters.PeakPagefileUsage;
	}
	
	return 0;
}

size_t WindowsMemoryTools::systemTotalPhysical()
{
	MEMORYSTATUSEX status;
	status.dwLength = sizeof(status);
	GlobalMemoryStatusEx( &status );
	return (size_t)status.ullTotalPhys;
}

size_t WindowsMemoryTools::systemAvailablePhysical()
{
	MEMORYSTATUSEX status;
	status.dwLength = sizeof(status);
	GlobalMemoryStatusEx( &status );
	return (size_t)status.ullAvailPhys;
}

#endif