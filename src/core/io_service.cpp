
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
bool IOService::saveAlignment(const std::string& file_name, vector<CGappedSequence*>& sequences, const std::string& format, seq_t seq_type, int no_threads, int gzip_level)
{
	if (format == "fasta"){
		return saveFASTAFormat(file_name, sequences, no_threads, gzip_level);
	} else if (format == "clustal"){
		return saveCLUSTALFormat(file_name, sequences, seq_type);
	} else {
		return saveFASTAFormat(file_name, sequences, no_threads, gzip_level);
	}		
}

// *******************************************************************
bool IOService::saveFASTAFormat(const std::string& file_name, vector<CGappedSequence*>& sequences, int no_threads, int gzip_level)
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


// BEGIN CLUSTAL OUTPUT SPECIFIC CODE
// *******************************************************************
FILE* openFile(const std::string& file_name, const char* mode = "r")
{
	FILE* file_ptr = std::fopen(file_name.c_str(), mode);
	if (!file_ptr) {
		std::cerr << "Error: could not open file \"" 
				  << file_name << "\" with mode \"" << mode << "\"\n";
		std::exit(EXIT_FAILURE);
	}
	return file_ptr;
}

// *******************************************************************
size_t utf8len(const char *s)
{
	size_t len = 0;
	for (; *s; ++s) if ((*s & 0xC0) != 0x80) ++len;
	return len;
}

// This was taken and modified from:
// https://github.com/GSLBiotech/clustal-omega/blob/d21fab82d380638c568c9427ed39cb42dd87d93b/src/squid/clustal.c#L191

// *******************************************************************
bool IOService::saveCLUSTALFormat(const std::string& file_name, vector<CGappedSequence*>& sequences, seq_t seq_type)
{
	int    idx;			/* counter for sequences         */
	int    len;			/* tmp variable for name lengths */
	int    namelen; /* maximum name length used      */
	int    pos;			/* position counter              */
	char  *buf;    	/* buffer for writing seq        */
	int    cpl = 60; //msa->alen<iWrap ? msa->alen+10 : (iWrap > 0 ? iWrap : 60);		/* char per line (< 64)          */

	/* consensus line stuff */
	int subpos;
	char first;
	int bail;
	int strong_bins[9];
	int weak_bins[11];
	/*int cons;*/
	int bin;
	int nseq = sequences.size(); // msa->nseq
	bool bResno = false;
	FILE *fp = openFile(file_name, "w");
	int *piResCnt = NULL;

	if (1 == bResno) {
		if (NULL == (piResCnt = (int *)malloc(nseq * sizeof(int)))) {
			printf("%s:%s:%d: could not malloc %d int for residue count", 
			__FUNCTION__, __FILE__, __LINE__, nseq);
			return false;
		} else {
			memset(piResCnt, 0, nseq * sizeof(int));
		}
	} /* do print residue numbers */

	if (NULL == (buf = (char *)malloc(cpl+20))) {
		printf("%s:%s:%d: could not malloc %d char for buffer",
		__FUNCTION__, __FILE__, __LINE__, cpl+20);
		return false;
	} else {
		memset(buf, 0, cpl+20);
	}

	/* calculate max namelen used */
	namelen = 0;
	for (idx = 0; idx < nseq; idx++)
	{
		/*if ((len = strlen(msa->sqname[idx])) > namelen) */ /* strlen() gives problems for unicode, FS, -> 290 */
		if ((len = utf8len(sequences[idx]->id.c_str())) > namelen)
		namelen = len; 
	}

	fprintf(fp, "FAMSA multiple sequence alignment\n");

	/*****************************************************
	 * Write the sequences
	 *****************************************************/

	fprintf(fp, "\n");	/* original had two blank lines */

	map<string, string> alignments;
	int alen = 0; // length of global alignment

	string firstseq;
	for (idx = 0; idx < nseq; idx++)
	{
		auto p = sequences[idx];
		string seq = p->Decode();
		string id = p->id;
		if (idx == 0){
			alen = seq.size();
			firstseq = seq;
		} 
		alignments[id] = seq;

		// cout << id << " -- " << alignments[id] << endl;
	}
	
	for (pos = 0; pos < alen; pos += cpl)
	{
		fprintf(fp, "\n");	/* Blank line between sequence blocks */
		for (idx = 0; idx < nseq; idx++)
		{
			auto p = sequences[idx];
			string id = p->id;
			const char *sqname = id.c_str();// msa->sqname[idx]
			const char *aseq = alignments[id].c_str(); // msa->aseq
			strncpy(buf, aseq + pos, cpl);

			buf[cpl] = '\0';
			if (1 == bResno) {
				char *pc = NULL;
				for (pc = buf; *pc != '\0'; pc++){
					if ( ( (*pc >= 'a') && (*pc <= 'z') ) || ( (*pc >= 'A') && (*pc <= 'Z') ) ) {
						piResCnt[idx]++;
					}
				}
				/* printf("%*s") gives problems for unicode, FS, -> 290 */
				/*fprintf(fp, "%-*s\t%s\t%d\n", namelen+5, msa->sqname[idx], buf, piResCnt[idx]);*/
				fprintf(fp, "%s%*s %s\t%d\n", sqname, (int)(namelen+5-utf8len(sqname)), "", buf, piResCnt[idx]);
			} else {
				/* printf("%*s") gives problems for unicode, FS, -> 290 */
				/*fprintf(fp, "%-*s\t%s\n", namelen+5, msa->sqname[idx], buf);*/
				fprintf(fp, "%s%*s %s\n", sqname, (int)(namelen+5-utf8len(sqname)), "", buf); 	
			}
		}
		/* do consensus dots */

		/* print namelen+5 spaces */
		for(subpos = 0; subpos <= namelen+5; subpos++)
			fprintf(fp, " ");

		const char *firstaseq = firstseq.c_str(); // msa->aseq
		for(subpos = pos; subpos < min(pos + cpl, alen); subpos++)
		{
			/* see if 100% conservation */
			first = firstaseq[subpos];
			bail = 0;
			for (idx = 1; idx < nseq; idx++)
			{
				auto p = sequences[idx];
				string id = p->id;
				const char *aseq = alignments[id].c_str(); // msa->aseq
				if(toupper(aseq[subpos]) != toupper(first)) { /* toupper makes consensus case-insensitive, FS, r290 -> */
					bail = 1;
					break;
				}
			}
			if(!bail)
				fprintf(fp, "*");
			else {
				/* if not then check strong */
				for(bin = 0; bin < 9; bin++)
					strong_bins[bin] = 0; /* clear the bins */

				for (idx = 0; (seq_type == seq_t::AA) && (idx < nseq); idx++)
				{ /* do this only for amino acids, no strong/weak consensus for nucleotide, FS, r290 -> */
					auto p = sequences[idx];
					string id = p->id;
					const char *aseq = alignments[id].c_str(); // msa->aseq
				
					switch(toupper(aseq[subpos])) /* toupper makes consensus case-insensitive, FS, r290 -> */
					{
						case 'S': strong_bins[0]++; break;
						case 'T': strong_bins[0]++; break;
						case 'A': strong_bins[0]++; break;
						case 'N': strong_bins[1]++; strong_bins[2]++; strong_bins[3]++; break;
						case 'E': strong_bins[1]++; strong_bins[3]++; break;
						case 'Q': strong_bins[1]++; strong_bins[2]++; strong_bins[3]++; strong_bins[4]++; break;
						case 'K': strong_bins[1]++; strong_bins[2]++; strong_bins[4]++; break;
						case 'D': strong_bins[3]++; break;
						case 'R': strong_bins[4]++; break;
						case 'H': strong_bins[2]++; strong_bins[4]++; strong_bins[7]++; break; /* added bin-2 (NHQK), FS 2016-07-14 */
						case 'M': strong_bins[5]++; strong_bins[6]++; break;
						case 'I': strong_bins[5]++; strong_bins[6]++; break;
						case 'L': strong_bins[5]++; strong_bins[6]++; break;
						case 'V': strong_bins[5]++; break;
						case 'F': strong_bins[6]++; strong_bins[8]++; break;
						case 'Y': strong_bins[7]++; strong_bins[8]++; break;
						case 'W': strong_bins[8]++; break;
					}
				}
				bail = 0;
				for(bin = 0; bin < 9; bin++)
					if(strong_bins[bin] == nseq) {
						bail = 1;
						break;
					}

				if(bail)
					fprintf(fp, ":");
				else {
					/* check weak */
					for(bin = 0; bin < 11; bin++)
						weak_bins[bin] = 0; /* clear the bins */

					for(idx = 0; (seq_type == seq_t::AA) && (idx < nseq); idx++)
					{ /* do this only for amino acids, no strong/weak consensus for nucleotide, FS, r290 -> */
						auto p = sequences[idx];
						string id = p->id;
						const char *aseq = alignments[id].c_str(); // msa->aseq
					
						switch(toupper(aseq[subpos])) /* toupper makes consensus case-insensitive, FS, r290 -> */
						{
							case 'C': weak_bins[0]++; break;
							case 'S': weak_bins[0]++; weak_bins[2]++; weak_bins[3]++; weak_bins[4]++; weak_bins[5]++; weak_bins[6]++; break;
							case 'A': weak_bins[0]++; weak_bins[1]++; weak_bins[2]++; weak_bins[4]++; break;
							case 'T': weak_bins[1]++; weak_bins[3]++; weak_bins[4]++; break;
							case 'V': weak_bins[1]++; weak_bins[9]++; break;
							case 'G': weak_bins[2]++; weak_bins[5]++; break; /* Added bin-5, FS 2016-07-14 */
							case 'N': weak_bins[3]++; weak_bins[5]++; weak_bins[6]++; weak_bins[7]++; weak_bins[8]++; break;
							case 'K': weak_bins[3]++; weak_bins[6]++; weak_bins[7]++; weak_bins[8]++; break;
							case 'D': weak_bins[5]++; weak_bins[6]++; weak_bins[7]++; break;
							case 'E': weak_bins[6]++; weak_bins[7]++; weak_bins[8]++; break;
							case 'Q': weak_bins[6]++; weak_bins[7]++; weak_bins[8]++; break;
							case 'H': weak_bins[7]++; weak_bins[8]++; weak_bins[10]++; break;
							case 'R': weak_bins[8]++; break;
							case 'F': weak_bins[9]++; weak_bins[10]++; break;
							case 'L': weak_bins[9]++; break;
							case 'I': weak_bins[9]++; break;
							case 'M': weak_bins[9]++; break;
							case 'Y': weak_bins[10]++; break;
						}
					}
					bail = 0;
					for(bin = 0; bin < 11; bin++)
						if(weak_bins[bin] == nseq) {
							bail = 1;
							break;
						}
					if(bail)
						fprintf(fp, ".");
					else
						fprintf(fp, " ");
				}
			}
		}
		fprintf(fp,"\n");
	} // end position for loop

	free(piResCnt); piResCnt = NULL;
	return true;
}
