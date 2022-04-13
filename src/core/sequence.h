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
	int sequence_no;
	uint32_t length;
	string id;
	vector<symbol_t> data;
	vector<bool> uppercase;
	
/*
#ifdef __APPLE__
	uint16_t hist[NO_SYMBOLS];
#else
	uint16_t alignas(32) hist[NO_SYMBOLS];
#endif
*/
	uint16_t* hist = nullptr;
	Array<bit_vec_t> bit_masks;
	
public:
	//CSequence();
	CSequence(const string& _id, const string& seq);
	
	// sequences are not copyable
	CSequence(const CSequence& x) noexcept = delete;
	CSequence& operator=(const CSequence& x) noexcept = delete;

	CSequence(CSequence&& x) noexcept = default;
	CSequence& operator=(CSequence&& x) noexcept = default; 
	
	~CSequence();

	

	void ComputeBitMasks();
	void ReleaseBitMasks();
	string DecodeSequence();
	void PrepareHistogram();
};

class CGappedSequence 
{
	static char mapping_table[25];

	void RecalculateDPS();
	void InitialiseDPS();

public:
	string id;
    vector<symbol_t> symbols;
	vector<bool> uppercase;
    size_t size;
	size_t gapped_size;

	size_t dps_size;
	size_t dps_size_div2;

    vector<int32_t> n_gaps;
    vector<int32_t> dps;        // dynamic position statistics (DSP) for the sequence

    CGappedSequence(CSequence &&_sequence);
    CGappedSequence(const CGappedSequence &_gapped_sequence);
    CGappedSequence(CGappedSequence &&_gapped_sequence) noexcept;
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
	void ClearDPS();
};
#endif