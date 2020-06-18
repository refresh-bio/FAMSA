/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

//#include "../libs/asmlib.h"

#include <cstring>

void mem_clear_avx(void* ptr, size_t size)
{
	memset(ptr, 0, size);
}

