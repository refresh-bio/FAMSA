/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#pragma once

#include "sequence.h"

#include <string>
#include <vector>


class IOService {

public:
	static size_t loadFasta(const std::string& file_name, std::vector<CSequence>& sequences);

	static bool saveAlignment(const std::string& file_name, const vector<CGappedSequence*> & sequences);

};
