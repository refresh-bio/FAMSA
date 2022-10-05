/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#ifndef _SEQUENCE_H
#define _SEQUENCE_H

#include "../core/defs.h"
#include "../utils/memory_monotonic.h"
#include "../utils/array.h"
#include <string>
#include <vector>
#include <array>


using namespace std;
using namespace refresh;

// *******************************************************************
class CSequence
{
	static char mapping_table[25];

public:
	uint32_t length;
	uint32_t data_size;
	symbol_t* data = nullptr;
	bit_vec_t *p_bit_masks;
	uint32_t p_bv_len;

	const int original_no;
	int sequence_no;
	string id;

	memory_monotonic_safe *mma;

	vector<bool> uppercase;
	vector<std::pair<int, char>> extra_symbols;
	
public:
	CSequence() = delete;
	CSequence(const string& _id, const string& seq, int sequence_no = -1, memory_monotonic_safe *mma = nullptr);
	
	// sequences are not copyable
	CSequence(const CSequence& x) noexcept = delete;
	CSequence& operator=(const CSequence& x) noexcept = delete;

	CSequence(CSequence&& x) noexcept;
	CSequence& operator=(CSequence&& x) noexcept = delete; 
	
	~CSequence();
	
	void DataResize(uint32_t new_size, symbol_t new_symbol);

	void ComputeBitMasks();
	void ReleaseBitMasks();
	//string DecodeSequence();

	memory_monotonic_safe* get_mma()
	{
		return mma;
	}
};

// *******************************************************************
struct CSequenceView
{
	uint32_t length;
	uint32_t padding1;
	symbol_t* data;
};

// *******************************************************************
class CGappedSequence
{
	static char mapping_table[25];

	void RecalculateDPS();
	void InitialiseDPS();

	memory_monotonic_safe* mma = nullptr;

public:
	symbol_t* symbols = nullptr;
    size_t size;
	size_t symbols_size;
	size_t gapped_size;

	size_t dps_size;
	size_t dps_size_div2;
	int original_no;
	int sequence_no;

    vector<uint32_t> n_gaps;
    vector<uint32_t> dps;        // dynamic position statistics (DSP) for the sequence

	string id;
	vector<bool> uppercase;
	vector<std::pair<int, char>> extra_symbols;

	CGappedSequence() = delete;
	CGappedSequence(const string& _id, const string& seq, int seq_no=-1, memory_monotonic_safe *mma=nullptr);
    CGappedSequence(CSequence &&_sequence);
    CGappedSequence(const CGappedSequence &_gapped_sequence);
    CGappedSequence(CGappedSequence &&_gapped_sequence) noexcept;
	~CGappedSequence();

	bool operator==(const CGappedSequence &gs) const;

    void InsertGap(uint32_t pos);
    void InsertGaps(uint32_t pos, uint32_t n);
	void InsertGapsVector(const vector<pair<uint32_t, uint32_t>> &v_gaps);

	void RemoveGap(size_t pos);
    void RemoveGaps(size_t pos, uint32_t n);
	symbol_t GetSymbol(size_t pos);

   // void DecodeRaw(symbol_t *seq);
    string Decode();
	uint32_t NoSymbols();

	void InsertFront(symbol_t new_symbol);

	void Clear();
	void ClearDPS();
};
#endif
