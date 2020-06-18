/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#ifndef _UTILS_H
#define _UTILS_H

void mem_clear(void* ptr, size_t size);
void mem_clear_avx(void* ptr, size_t size);
void mem_clear_avx2(void* ptr, size_t size);

#endif