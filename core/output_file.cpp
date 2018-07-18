/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#include "../core/output_file.h"

#include <iostream>
#include <fstream>

// ****************************************************************************
COutputFile::COutputFile()
{
}


// ****************************************************************************
COutputFile::~COutputFile()
{
}

// ****************************************************************************
bool COutputFile::SaveFile(string file_name)
{
	ostream* out;
	ofstream outfile;

	if (file_name == "STDOUT") {
		out = &cout;
	}
	else {
		outfile.open(file_name.c_str(), ios_base::out);
		out = &outfile;
	}

	string s;
	bool is_id = false;

	string id, seq;

	for(auto &p : sequences)
	{
		*out << p->id << "\n";
		string seq = p->Decode();
		for(size_t pos = 0; pos < seq.size(); pos += 60)
			*out << seq.substr(pos, 60) << "\n";
	}

	return true;
}

// ****************************************************************************
void COutputFile::PutSequences(vector<CGappedSequence*> &_sequences)
{
	sequences = _sequences;
}

// ****************************************************************************
void COutputFile::PutSequences(vector<CGappedSequence*> &&_sequences)
{
	sequences = std::move(_sequences);
}
