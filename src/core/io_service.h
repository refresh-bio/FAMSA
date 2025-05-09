/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#pragma once
#include "../libs/refresh/compression/lib/file_wrapper.h"

#include "../core/sequence.h"
#include <algorithm>

#include <string>
#include <vector>


class IOService {

public:
//	template <class seq_t>
//	static size_t loadFasta_old(const std::string& file_name, std::vector<seq_t>& sequences, memory_monotonic_safe* mma = nullptr);

	template <class seq_t>
	static size_t loadFasta(const std::string& file_name, std::vector<seq_t>& sequences, memory_monotonic_safe* mma = nullptr);

	static bool saveAlignment(const std::string& file_name, vector<CGappedSequence*> & sequences, int no_threads, int gzip_level, float rare_column_threshold);
};

#if 0
// *******************************************************************
template <class seq_t>
size_t IOService::loadFasta_old(
	const std::string& file_name, std::vector<seq_t>& sequences, memory_monotonic_safe* mma) {

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
	int seq_no = 0;

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
				sequences.emplace_back(id, seq, seq_no++, mma);
				seq.clear();
			}
			id = s;
		}
		else {
			seq += s;
		}
	}

	if (!id.empty() && !seq.empty())
		sequences.emplace_back(id, seq, seq_no++, mma);

	return sequences.size();
}
#endif

// *******************************************************************
template <class seq_t>
size_t IOService::loadFasta(
	const std::string& file_name, std::vector<seq_t>& sequences, memory_monotonic_safe* mma) {

	refresh::stream_in_file msgz(file_name == "STDIN" ? ""s : file_name);

	if (!msgz.is_open()) {
		throw std::runtime_error("Unable to open input file " + file_name);
	}

	refresh::stream_decompression sdf(&msgz);

	string s;
	string id, seq;
	int seq_no = 0;

	while (!sdf.eof())
	{
		sdf.getline(s);

		while (!s.empty() && (s.back() == '\n' || s.back() == '\r'))
			s.pop_back();
		if (s.empty())
			continue;

		if (s[0] == '>')
		{
			if (!id.empty() && !seq.empty())
			{
				sequences.emplace_back(id, seq, seq_no++, mma);
				seq.clear();
			}
			id = s;
		}
		else {
			seq += s;
		}
	}

	if (!id.empty() && !seq.empty())
		sequences.emplace_back(id, seq, seq_no++, mma);

	return sequences.size();
}
