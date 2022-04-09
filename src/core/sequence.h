/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#ifndef _SEQUENCE_H
#define _SEQUENCE_H

#include "../core/defs.h"
#include "../utils/array.h"
#include <string>
#include <vector>
#include <array>


using namespace std;

class CSequence
{
	static char mapping_table[25];

public:
	int sequence_no = -1;

	string id;
	vector<symbol_t> data;
	vector<bool> uppercase;
	uint32_t length;

#ifdef __APPLE__
	uint16_t hist[NO_SYMBOLS];
#else
	uint16_t alignas(32) hist[NO_SYMBOLS];
#endif

	Array<bit_vec_t> bit_masks;
	
public:
	CSequence();
	CSequence(const string& _id, const string& seq);
	CSequence(const CSequence& x);
	CSequence(CSequence&& x) noexcept;
	~CSequence();

	CSequence& operator=(const CSequence& x) noexcept;

	void ComputeBitMasks();
	void ReleaseBitMasks();
	string DecodeSequence();
	void PrepareHistogram();
};

#define COMPACT_GAPPED_SEQUENCE

#ifndef COMPACT_GAPPED_SEQUENCE
class CGappedSequence 
{
	static char mapping_table[25];

public:
	string id;
    vector<symbol_t> symbols;
	vector<bool> uppercase;
    size_t size;
	size_t dps_size;

    vector<int32_t> n_gaps;
    vector<int32_t> dps;        // dynamic position statistics (DSP) for the sequence
    size_t gapped_size;

    CGappedSequence(CSequence &_sequence);
    CGappedSequence(const CGappedSequence &_gapped_sequence);
    ~CGappedSequence();

	bool operator==(const CGappedSequence &gs);

	void InitialiseDPS();

    void InsertGap(size_t pos);
    void InsertGaps(size_t pos, uint32_t n);
	void InsertGapsVector(const vector<pair<uint32_t, uint32_t>> &v_gaps);

	void RemoveGap(size_t pos);
    void RemoveGaps(size_t pos, uint32_t n);
	symbol_t GetSymbol(size_t pos);

    void DecodeRaw(symbol_t *seq);
    string Decode();
	uint32_t NoSymbols();

	void AddGuard() { symbols.insert(symbols.begin(), GUARD); };

	void Clear();
};
#else

class CGappedSequence 
{
	static char mapping_table[25];

	void RecalculateDPS();
	void InitialiseDPS();
	void FillDPS();

public:
	string id;
    vector<symbol_t> symbols;
	vector<bool> uppercase;
    size_t size;
	size_t dps_size;
	size_t dps_size_div2;

    vector<int32_t> n_gaps;
    vector<int32_t> dps;        // dynamic position statistics (DSP) for the sequence

    size_t gapped_size;

    CGappedSequence(CSequence &_sequence);
    CGappedSequence(const CGappedSequence &_gapped_sequence);
    CGappedSequence(CGappedSequence &&_gapped_sequence);
    ~CGappedSequence();

	bool operator==(const CGappedSequence &gs);

    void InsertGap(size_t pos);
    void InsertGaps(size_t pos, uint32_t n);
	void InsertGapsVector(const vector<pair<uint32_t, uint32_t>> &v_gaps);

	void RemoveGap(size_t pos);
    void RemoveGaps(size_t pos, uint32_t n);
	symbol_t GetSymbol(size_t pos);

    void DecodeRaw(symbol_t *seq);
    string Decode();
	uint32_t NoSymbols();

	void AddGuard() { symbols.insert(symbols.begin(), GUARD); };

	void Clear();
};
#endif

#endif
