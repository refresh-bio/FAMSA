#pragma once
#ifdef _DEBUG
#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define DEBUG_NEW new (_NORMAL_BLOCK, __FILE__, __LINE__)
#else
#define DEBUG_NEW new
#endif

#ifdef _DEBUG
// uncomment this line to perform new's with memory dumps
#define new DEBUG_NEW
#endif
