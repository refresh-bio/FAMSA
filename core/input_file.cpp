/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#include "../core/input_file.h"

#include <iostream>
#include <fstream>

// *******************************************************************
CInputFile::CInputFile()
{
}

// *******************************************************************
CInputFile::~CInputFile()
{
}

// *******************************************************************
bool CInputFile::ReadFile(string file_name)
{
	istream* in;

	ifstream infile;

	if (file_name == "STDIN") {
		in = &cin;
	}
	else {
		infile.open(file_name.c_str(), ios_base::in);
		if (!infile.good())
			return false;
		in = &infile;
	}

	string s;
	string id, seq;

	while(in->good())
	{
		getline(*in, s);
		while(!s.empty() && (s[s.length()-1] == '\n' || s[s.length()-1] == '\r'))
			s.pop_back();
		if(s.empty())
			continue;

		if(s[0] == '>')
		{
			if(!id.empty() && !seq.empty())
			{
				sequences.emplace_back(id, seq);
			//	sequences.push_back(CSequence(id, seq));
				seq.clear();
			}
			id = s;
		}
		else
			seq += s;
	}

	if(!id.empty() && !seq.empty())
		sequences.push_back(CSequence(id, seq));

	return sequences.size() > 1;
}

// *******************************************************************
void CInputFile::GetSequences(vector<CSequence> &_sequences) const
{
	_sequences = sequences;
}

// *******************************************************************
void CInputFile::StealSequences(vector<CSequence> &_sequences)
{
	_sequences = std::move(sequences);
}
