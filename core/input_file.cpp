/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

Version: 1.1
Date   : 2016-06-29
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
	ifstream in;
	string s;
	bool is_id = false;

	string id, seq;

	in.open(file_name.c_str(), ios_base::in);

	if(!in.good())
		return false;

	while(in.good())
	{
		getline(in, s);
		while(!s.empty() && (s[s.length()-1] == '\n' || s[s.length()-1] == '\r'))
			s.pop_back();
		if(s.empty())
			continue;

		if(s[0] == '>')
		{
			if(!id.empty())
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

	if(!id.empty())
		sequences.push_back(CSequence(id, seq));

	return true;
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
