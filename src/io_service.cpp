/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#include "io_service.h"

#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;

// *******************************************************************
size_t IOService::loadFasta(const std::string& file_name, std::vector<CSequence>& sequences)
{
	istream* in;
	ifstream infile;

	if (file_name == "STDIN") {
		in = &cin;
	}
	else {
		infile.open(file_name.c_str(), ios_base::in);
		if (!infile.good())
			return 0;
		in = &infile;
	}

	string s;
	string id, seq;

	while (in->good())
	{
		getline(*in, s);
		while (!s.empty() && (s[s.length() - 1] == '\n' || s[s.length() - 1] == '\r'))
			s.pop_back();
		if (s.empty())
			continue;

		if (s[0] == '>')
		{
			if (!id.empty() && !seq.empty())
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

	if (!id.empty() && !seq.empty())
		sequences.push_back(CSequence(id, seq));

	return sequences.size();

}

// *******************************************************************
bool IOService::saveAlignment(const std::string& file_name, const vector<CGappedSequence*> & sequences)
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
	string id, seq;

	for (auto &p : sequences)
	{
		*out << p->id << "\n";
		string seq = p->Decode();
		for (size_t pos = 0; pos < seq.size(); pos += 60)
			*out << seq.substr(pos, 60) << "\n";
	}

	return true;
}


// *******************************************************************
/*
bool IOService::exportDistanceMatrix(const std::string& file_name, float*matrix, size_t size) {
	ofstream file(file_name);
	if (!file) {
		return false;
	}

	for (size_t i = 0; i < size; ++i) {
		for (size_t j = 0; j < i; ++j) {
			const size_t id = UPGMA_TriangleSubscript(i, j);
			float d = matrix[id];
			file << d << ", ";
		}
		file << std::endl;
	}

	return true;
}
*/
