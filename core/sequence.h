/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#ifndef _SEQUENCE_H
#define _SEQUENCE_H

#include "../core/defs.h"
#include "../core/params.h"
#include <string>
#include <vector>

using namespace std;

class CSequence
{
	static char mapping_table[25];

public:
	string id;
	vector<symbol_t> data;
	vector<bool> uppercase;
	size_t length;

	vector<vector<bit_vec_t>> bit_masks;

	void ComputeBitMasks();
	void ComputeBitMasks64();

public:
	CSequence(const string& _id, const string& seq);
	CSequence(string _id, vector<symbol_t> &symbol_seq);
	CSequence(const CSequence &_sequence);

	CSequence(CSequence&& rhs);
	CSequence& operator=(CSequence&& rhs);
	CSequence& operator=(const CSequence& rhs);

	~CSequence();

	string DecodeSequence();
};

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
	bool operator!=(const CGappedSequence &gs);

	void InitialiseDPS();

    void InsertGap(size_t pos);
    void InsertGaps(size_t pos, uint32_t n);
    void RemoveGap(size_t pos);
    void RemoveGaps(size_t pos, uint32_t n);
	symbol_t GetSymbol(size_t pos);

    void DecodeRaw(symbol_t *seq);
    string Decode();
};

#endif