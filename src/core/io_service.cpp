
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
bool IOService::saveAlignment(const std::string& file_name, vector<CGappedSequence*>& sequences, int no_threads, int gzip_level, float rare_column_threshold)
{
	string s;
	string id, seq;

	int pack_size = gzip_level < 0 ? 5 : 10;
	int clear_pack_size = 100;

	atomic<int> seq_id {0};
	vector<thread> v_threads;
	v_threads.reserve(no_threads);

	// *** Clear DPS data
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

	vector<bool> accepted_columns(sequences[0]->gapped_size, true);

	// *** Find rare columns
	if (rare_column_threshold < 1.0)
	{
		seq_id = 0;
		int rare_column_pack_size = 32;

		vector<uint32_t> global_column_freq(sequences[0]->gapped_size, 0);
		mutex mtx;

		for (int i = 0; i < no_threads; ++i)
			v_threads.emplace_back([rare_column_pack_size, &seq_id, &sequences, &mtx, &global_column_freq, rare_column_threshold] {

			vector<uint32_t> local_column_freq(sequences[0]->gapped_size, 0);

			int no_seqs = static_cast<int>(sequences.size());
			while (true)
			{
				int32_t id_from = seq_id.fetch_add(rare_column_pack_size);
				int32_t id_to = id_from + rare_column_pack_size;
				if (id_from >= no_seqs)
					break;
				if (id_to >= no_seqs)
					id_to = no_seqs;

				for (int i = id_from; i < id_to; ++i)
				{
					auto p = sequences[i];
					auto& symbols = p->symbols;
					auto& n_gaps = p->n_gaps;
					int pos = n_gaps[0];

					for (uint32_t j = 1; j <= p->size; ++j)
					{
						local_column_freq[pos]++;
						pos += n_gaps[j] + 1;
					}
				}
			}

			lock_guard<mutex> lck(mtx);
			for (int i = 0; i < sequences[0]->gapped_size; ++i)
				global_column_freq[i] += local_column_freq[i];
				});

		for (auto& t : v_threads)
			t.join();
		v_threads.clear();

		uint32_t no_seqs = static_cast<uint32_t>(sequences.size());
		uint32_t min_no_symbols = static_cast<uint32_t>(no_seqs * rare_column_threshold);

		for (size_t i = 0; i < global_column_freq.size(); ++i)
			accepted_columns[i] = global_column_freq[i] >= min_no_symbols;
	}

	// *** Save sequences
	seq_id = 0;
	CLimitedPriorityQueue<vector<uint8_t>> lpq(no_threads, 5 * no_threads);

	for (int i = 0; i < no_threads; ++i)
		v_threads.emplace_back([pack_size, gzip_level, rare_column_threshold, &seq_id, &sequences, &lpq, &accepted_columns] {
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

			string seq_densified;

			for (int i = id_from; i < id_to; ++i)
			{
				auto p = sequences[i];

				string seq = p->Decode();
				s_tmp.append(p->id);
				s_tmp.push_back('\n');

				if (rare_column_threshold < 1.0)
				{
					swap(seq, seq_densified);
					seq.clear();

					seq.reserve(seq_densified.size());
					
					for (size_t j = 0; j < seq_densified.size(); ++j)
						if (accepted_columns[j])
							seq.push_back(seq_densified[j]);
				}

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
