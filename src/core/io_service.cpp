/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#include "../core/io_service.h"
#include "../core/queues.h"
#include <libdeflate/libdeflate.h>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include <fstream>
#include <iostream>
#include <algorithm>
#include <thread>
#include <atomic>

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

#if 0
// *******************************************************************
bool IOService::saveAlignmentOld(const std::string& file_name, const vector<CGappedSequence*> & sequences)
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
bool IOService::saveGzippedAlignment(const std::string& file_name, vector<CGappedSequence*>& sequences, int no_threads, int gzip_level)
{
	string s;
	string id, seq;

	const uint32_t pack_size = 10;

	atomic<int> seq_id;

	vector<thread> v_threads;
	v_threads.reserve(no_threads);

	CPriorityQueue<vector<uint8_t>> pq(no_threads);

	for (int i = 0; i < no_threads; ++i)
		v_threads.emplace_back([pack_size, gzip_level, &seq_id, &sequences, &pq] {
		int no_seqs = static_cast<int>(sequences.size());
		string s_tmp;

		libdeflate_compressor* compressor = libdeflate_alloc_compressor(gzip_level);
		vector<uint8_t> gz_vec;

		while (true)
		{
			uint32_t id_from = seq_id.fetch_add(pack_size);
			uint32_t id_to = id_from + pack_size;

			if (id_from >= no_seqs)
				break;
			if (id_to >= no_seqs)
				id_to = no_seqs;

			s_tmp.clear();

			for (uint32_t i = id_from; i < id_to; ++i)
			{
				auto p = sequences[i];

				string seq = p->Decode();
				s_tmp.append(p->id);
				s_tmp.push_back('\n');

				size_t seq_size = seq.size();
				auto ptr = seq.c_str();
				size_t step;

				for (size_t pos = 0; pos < seq_size; pos += step, ptr += step)
				{
					step = 60;
					if (pos + step > seq_size)
						step = seq_size - pos;

					s_tmp.append(ptr, step);
					s_tmp.push_back('\n');
				}

				// Remove CGappedSequence* here to save memory
				p->Clear();
			}

			size_t need_alloc = libdeflate_gzip_compress_bound(compressor, s_tmp.size());
			if (gz_vec.size() < need_alloc)
				gz_vec.resize(need_alloc);

			size_t gzipped_size = libdeflate_gzip_compress(compressor, s_tmp.data(), s_tmp.size(), gz_vec.data(), gz_vec.size());

			vector<uint8_t> v(gz_vec.begin(), gz_vec.begin() + gzipped_size);
			//v_gz_sequences[id_from / pack_size].assign(gz_vec.begin(), gz_vec.begin() + gzipped_size);
			pq.Emplace(id_from / pack_size, move(v));
		}

		pq.MarkCompleted();

		libdeflate_free_compressor(compressor);

		});

	if (file_name == "STDOUT")
	{
#ifdef _WIN32
		_setmode(_fileno(stdout), _O_BINARY);
#endif
		vector<uint8_t> dat;

		while (!pq.IsCompleted())
		{
			if (!pq.Pop(dat))
				continue;

			fwrite(dat.data(), 1, dat.size(), stdout);
		}
	}
	else
	{
		ofstream outfile;
		const size_t BUFFER_SIZE = 128 << 20;
		char* buffer = new char[BUFFER_SIZE];

		outfile.open(file_name.c_str(), ios_base::out | ios_base::binary);
		outfile.rdbuf()->pubsetbuf(buffer, BUFFER_SIZE);

		vector<uint8_t> dat;

		while (!pq.IsCompleted())
		{
			if (!pq.Pop(dat))
				continue;

			outfile.write((char*)dat.data(), dat.size());
		}

		outfile.close();
		delete[] buffer;
	}

	for (auto& t : v_threads)
		t.join();

	return true;
}
#endif

// *******************************************************************
bool IOService::saveAlignment(const std::string& file_name, vector<CGappedSequence*>& sequences, int no_threads, int gzip_level)
{
	string s;
	string id, seq;

	uint32_t pack_size = gzip_level < 0 ? 5 : 10;

	atomic<int> seq_id {0};

	vector<thread> v_threads;
	v_threads.reserve(no_threads);

	CLimitedPriorityQueue<vector<uint8_t>> lpq(no_threads, 5 * no_threads);

	for (int i = 0; i < no_threads; ++i)
		v_threads.emplace_back([pack_size, gzip_level, &seq_id, &sequences, &lpq] {
		int no_seqs = static_cast<int>(sequences.size());
		string s_tmp;

		libdeflate_compressor* compressor = gzip_level >= 0 ? libdeflate_alloc_compressor(gzip_level) : nullptr;
		vector<uint8_t> gz_vec, raw_vec;

		while (true)
		{
			uint32_t id_from = seq_id.fetch_add(pack_size);
			uint32_t id_to = id_from + pack_size;

			if (id_from >= no_seqs)
				break;
			if (id_to >= no_seqs)
				id_to = no_seqs;

			s_tmp.clear();

			for (uint32_t i = id_from; i < id_to; ++i)
			{
				auto p = sequences[i];

				string seq = p->Decode();
				s_tmp.append(p->id);
				s_tmp.push_back('\n');

				size_t seq_size = seq.size();
				auto ptr = seq.c_str();
				size_t step;

				for (size_t pos = 0; pos < seq_size; pos += step, ptr += step)
				{
					step = 60;
					if (pos + step > seq_size)
						step = seq_size - pos;

					s_tmp.append(ptr, step);
					s_tmp.push_back('\n');
				}

				// Remove CGappedSequence* here to save memory
				p->Clear();
			}

			if (gzip_level >= 0)
			{
				size_t need_alloc = libdeflate_gzip_compress_bound(compressor, s_tmp.size());
				if (gz_vec.size() < need_alloc)
					gz_vec.resize(need_alloc);

				size_t gzipped_size = libdeflate_gzip_compress(compressor, s_tmp.data(), s_tmp.size(), gz_vec.data(), gz_vec.size());

				vector<uint8_t> v(gz_vec.begin(), gz_vec.begin() + gzipped_size);
				//v_gz_sequences[id_from / pack_size].assign(gz_vec.begin(), gz_vec.begin() + gzipped_size);
				lpq.Emplace(id_from / pack_size, move(v));
			}
			else
			{
				raw_vec.assign(s_tmp.begin(), s_tmp.end());
				lpq.Emplace(id_from / pack_size, move(raw_vec));
			}
		}

		lpq.MarkCompleted();

		if(gzip_level >= 0)
			libdeflate_free_compressor(compressor);

			});

	if (file_name == "STDOUT")
	{
#ifdef _WIN32
		_setmode(_fileno(stdout), _O_BINARY);
#endif
		vector<uint8_t> dat;

		while (!lpq.IsCompleted())
		{
			if (!lpq.Pop(dat))
				continue;

			fwrite(dat.data(), 1, dat.size(), stdout);
		}
	}
	else
	{
		ofstream outfile;
		const size_t BUFFER_SIZE = 128 << 20;
		char* buffer = new char[BUFFER_SIZE];

		outfile.open(file_name.c_str(), ios_base::out | ios_base::binary);
		outfile.rdbuf()->pubsetbuf(buffer, BUFFER_SIZE);

		vector<uint8_t> dat;

		while (!lpq.IsCompleted())
		{
			if (!lpq.Pop(dat))
				continue;

			outfile.write((char*)dat.data(), dat.size());
		}

		outfile.close();
		delete[] buffer;
	}

	for (auto& t : v_threads)
		t.join();

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