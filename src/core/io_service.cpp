
/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

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
bool IOService::saveAlignment(const std::string& file_name, vector<CGappedSequence*>& sequences, int no_threads, int gzip_level)
{
	string s;
	string id, seq;

	int pack_size = gzip_level < 0 ? 5 : 10;
	int clear_pack_size = 100;

	atomic<int> seq_id {0};
	vector<thread> v_threads;
	v_threads.reserve(no_threads);

	for (int i = 0; i < no_threads; ++i)
		v_threads.emplace_back([clear_pack_size, &seq_id, &sequences] {
		int no_seqs = static_cast<int>(sequences.size());

		while (true)
		{
			int32_t id_from = seq_id.fetch_add(clear_pack_size);
			int32_t id_to = id_from + clear_pack_size;

			if (id_from >= no_seqs)
				break;
			if (id_to >= no_seqs)
				id_to = no_seqs;

			for (int i = id_from; i < id_to; ++i)
				sequences[i]->ClearDPS();
		}
		});

	for (auto& t : v_threads)
		t.join();

	v_threads.clear();

	seq_id = 0;

	CLimitedPriorityQueue<vector<uint8_t>> lpq(no_threads, 5 * no_threads);

	for (int i = 0; i < no_threads; ++i)
		v_threads.emplace_back([pack_size, gzip_level, &seq_id, &sequences, &lpq] {
		int no_seqs = static_cast<int>(sequences.size());
		string s_tmp;

		libdeflate_compressor* compressor = gzip_level >= 0 ? libdeflate_alloc_compressor(gzip_level) : nullptr;
		vector<uint8_t> gz_vec, raw_vec;

		while (true)
		{
			int id_from = seq_id.fetch_add(pack_size);
			int id_to = id_from + pack_size;

			if (id_from >= no_seqs)
				break;
			if (id_to >= no_seqs)
				id_to = no_seqs;

			s_tmp.clear();

			for (int i = id_from; i < id_to; ++i)
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

				// Clear internal data here to save memory
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
