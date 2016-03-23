

#if defined(_WIN32)
#include <Windows.h>
#include <Psapi.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))


#else

#endif

#include "MemoryTools.h"



size_t MemoryTools::processCurrentVirtual()
{
#if defined(__linux__)
    // Ugh, getrusage doesn't work well on Linux.  Try grabbing info
    // directly from the /proc pseudo-filesystem.  Reading from
    // /proc/self/statm gives info on your own process, as one line of
    // numbers that are: virtual mem program size, resident set size,
    // shared pages, text/code, data/stack, library, dirty pages.  The
    // mem sizes should all be multiplied by the page size.
   

#elif defined(_WINDOWS)
    // According to MSDN...
  

#else
    // No idea what platform this is
    return 0;   // Punt
#endif
}


size_t MemoryTools::processPeakVirtual()
{
#if defined(_WINDOWS)
	// According to MSDN...
	

#else
// No idea what platform this is
return 0;   // Punt
#endif
}


/**
 * Returns the size of physical memory (RAM) in bytes.
 */
size_t MemoryTools::systemTotalPhysical( )
{


#if defined(_WIN32)
	/* Windows. ------------------------------------------------- */
	/* Use new 64-bit MEMORYSTATUSEX, not old 32-bit MEMORYSTATUS */
	

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
	

#else
	return 0L;			/* Unknown OS. */
#endif
}

