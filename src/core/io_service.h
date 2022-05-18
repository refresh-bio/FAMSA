/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#pragma once

#include "../core/sequence.h"

#include <string>
#include <vector>


class IOService {

public:
	static size_t loadFasta(const std::string& file_name, std::vector<CSequence>& sequences, memory_monotonic_safe* mma = nullptr);
	static bool saveAlignment(const std::string& file_name, vector<CGappedSequence*> & sequences, int no_threads, int gzip_level);
};