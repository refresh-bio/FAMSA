/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#include "../core/io_service.h"

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
//	ostream* out;
//	ofstream outfile;

	string s;
	string id, seq;

/*	if (file_name == "STDOUT") {
		out = &cout;
	}
	else {
		outfile.open(file_name.c_str(), ios_base::out);
		out = &outfile;
	}

	for (auto &p : sequences)
	{
		*out << p->id << "\n";
		string seq = p->Decode();
		for (size_t pos = 0; pos < seq.size(); pos += 60)
			*out << seq.substr(pos, 60) << "\n";
	}*/

	if (file_name == "STDOUT")
	{
		for (auto& p : sequences)
		{
			cout << p->id << "\n";
			string seq = p->Decode();
			for (size_t pos = 0; pos < seq.size(); pos += 60)
				cout << seq.substr(pos, 60) << "\n";
		}
	}
	else
	{
		ofstream outfile;
		const size_t BUFFER_SIZE = 128 << 20;
		char *buffer = new char[BUFFER_SIZE];

		outfile.open(file_name.c_str(), ios_base::out);
		outfile.rdbuf()->pubsetbuf(buffer, BUFFER_SIZE);

		for (auto& p : sequences)
		{
			outfile.write(p->id.c_str(), p->id.size());
			outfile.put('\n');

			string seq = p->Decode();
			size_t seq_size = seq.size();
			auto ptr = seq.c_str();
			size_t step;

			for (size_t pos = 0; pos < seq_size; pos += step, ptr += step)
			{
				step = 60;
				if (pos + step > seq_size)
					step = seq_size - pos;

				outfile.write(ptr, step);
				outfile.put('\n');
			}
		}

		outfile.close();
		delete[] buffer;
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