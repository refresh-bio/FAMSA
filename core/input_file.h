/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#ifndef _INPUT_FILE_H
#define _INPUT_FILE_H

#include "../core/defs.h"
#include "../core/sequence.h"

#include <string>
#include <vector>

using namespace std;

class CInputFile 
{
	vector<CSequence> sequences;

public:
	CInputFile();
	~CInputFile();

	bool ReadFile(string file_name);

	void GetSequences(vector<CSequence> &_sequences) const;
	void StealSequences(vector<CSequence> &_sequences);
};

#endif
