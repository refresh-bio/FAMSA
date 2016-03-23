#pragma once
#include <cstddef>

class WindowsMemoryTools
{
public:
	static size_t processCurrentVirtual();
	static size_t processPeakVirtual();
	static size_t systemTotalPhysical();
	static size_t systemAvailablePhysical();
};