/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#include "sequence.h"

#include <algorithm>
#include <assert.h>

using namespace std;

char CSequence::mapping_table[25] = "ARNDCQEGHILKMFPSTWYVBZX*";
char CGappedSequence::mapping_table[25] = "ARNDCQEGHILKMFPSTWYVBZX*";


// *******************************************************************
CSequence::CSequence(const string& _id, const string& seq) : id(_id), bit_masks(nullptr)
{
	length = seq.length();
	data.resize(length);
	uppercase.resize(length);

	for(size_t i = 0; i < length; ++i)
	{
		char c = seq[i];
		if(c > 'Z')
		{
			c -= 32;
			uppercase[i] = false;
		}
		else
			uppercase[i] = true;

		char *q = find(mapping_table, mapping_table+25, c);
		if(q == mapping_table+25)
			data[i] = (symbol_t) UNKNOWN_SYMBOL;
		else
			data[i] = (symbol_t) (q - mapping_table);
	}
}


// *******************************************************************
CSequence::~CSequence()
{
	ReleaseBitMasks();
}


// *******************************************************************
string CSequence::DecodeSequence()
{
	string s;

	s.reserve(data.size()+1);
	for(auto &p : data)
		if(p == GUARD)
			continue;
		else if(p == GAP)
			s += "-";
		else
			s += mapping_table[p];

	return s;
}

// *******************************************************************
// Compute Bit-mask Vectors for bit-parallel computation of LCS for the sequences
void CSequence::ComputeBitMasks()
{
	size_t bv_len = (data.size() + bv_size - 1) / bv_size;

	/*
	bit_masks.clear();
	bit_masks.resize(NO_SYMBOLS);

	for (size_t i = 0; i < NO_SYMBOLS; ++i)
		bit_masks[i].resize(bv_len, (bit_vec_t) 0);
	*/

	bit_masks = new Array<bit_vec_t>(bv_len, NO_SYMBOLS, (bit_vec_t)0);

	for(size_t i = 0; i < length; ++i)
		if(data[i] >= 0 && data[i] < NO_VALID_AMINOACIDS)
			(*bit_masks)[data[i]][i / bv_size] |= ((bit_vec_t) 1) << (i % bv_size);
}

// *******************************************************************
void CSequence::ReleaseBitMasks()
{
	delete bit_masks;
	bit_masks = nullptr;
}

// *******************************************************************
//
// *******************************************************************

// *******************************************************************
CGappedSequence::CGappedSequence(CSequence &_sequence)
{
	id = _sequence.id;
    symbols = _sequence.data;
    size    = symbols.size();
	uppercase = _sequence.uppercase;

    gapped_size = size;
    n_gaps.resize(size+1, 0);

	InitialiseDPS();
}

// *******************************************************************
CGappedSequence::CGappedSequence(const CGappedSequence &_gapped_sequence)
{
	id          = _gapped_sequence.id;

    size        = _gapped_sequence.size;
    gapped_size = _gapped_sequence.gapped_size;
	dps_size    = _gapped_sequence.dps_size;

    symbols     = _gapped_sequence.symbols;
    n_gaps      = _gapped_sequence.n_gaps;
    dps         = _gapped_sequence.dps;
	uppercase   = _gapped_sequence.uppercase;
}

// *******************************************************************
CGappedSequence::~CGappedSequence()
{
}

// *******************************************************************
bool CGappedSequence::operator==(const CGappedSequence &gs)
{
	return id == gs.id &&
		gapped_size == gs.gapped_size &&
		symbols == gs.symbols &&
		n_gaps == gs.n_gaps;
}

// *******************************************************************
bool CGappedSequence::operator!=(const CGappedSequence &gs)
{
	return !(*this == gs);
}

// *******************************************************************
void CGappedSequence::InitialiseDPS()
{
	// Round size of n_gap to the nearest not smaller power of 2
	dps_size = size+1;
	if((dps_size & (dps_size - 1)) != 0)
	{
		while((dps_size & (dps_size - 1)) != 0)
			dps_size &= dps_size - 1;
		dps_size <<= 1;
	}

	dps.resize(2 * dps_size, 0);

	for(size_t i = 0; i <= size; ++i)
		dps[dps_size+i] = n_gaps[i] + 1;

	for(size_t n = dps_size / 2; n; n /= 2)
		for(size_t i = 0; i < n; ++i)
			dps[n + i] = dps[2*n + 2*i] + dps[2*n + 2*i + 1];
}

// *******************************************************************
string CGappedSequence::Decode()
{
    string s;

    s.reserve(gapped_size);

	// Starting gaps
    for(int j = 0; j < n_gaps[0]; ++j)
        s.push_back('-');
    for(int i = 1; i <= size; ++i)
	{
		char symbol = mapping_table[symbols[i]];

		if(!uppercase[i-1])
			symbol += 32;		// change to lowercase

        s.push_back(symbol);
        for(int j = 0; j < n_gaps[i]; ++j)
            s.push_back('-');
	}

	return std::move(s);
}

// *******************************************************************
void CGappedSequence::DecodeRaw(symbol_t *seq)
{
	uint32_t seq_pos = 1;

	seq[0] = GUARD;

	// Starting gaps
    for(int j = 0; j < n_gaps[0]; ++j)
		seq[seq_pos++] = GAP;
    for(int i = 1; i <= size; ++i)
	{
		seq[seq_pos++] = symbols[i];

		for(int j = 0; j < n_gaps[i]; ++j)
            seq[seq_pos++] = GAP;
	}
}

// *******************************************************************
void CGappedSequence::InsertGap(size_t pos)
{
 	// Look for the place to insert the gap
	size_t x = 1;

	while(x < dps_size)
	{
		if(dps[2 * x] >= pos)
			x = 2 * x;
		else
		{
			pos -= dps[2 * x];
			x = 2 * x + 1;
		}
	}

	// Increment the no. of gaps
	++n_gaps[x - dps_size];

	// Update DSP
	for(; x; x /= 2)
		++dps[x];

    ++gapped_size;
}

// *******************************************************************
void CGappedSequence::InsertGaps(size_t pos, uint32_t n)
{
 	// Look for the place to insert the gap
	size_t x = 1;

	while(x < dps_size)
	{
		if(dps[2 * x] >= pos)
			x = 2 * x;
		else
		{
			pos -= dps[2 * x];
			x = 2 * x + 1;
		}
	}

	// Increment the no. of gaps
	n_gaps[x - dps_size] += n;

	// Update DSP
	for(; x; x /= 2)
		dps[x] += n;

    gapped_size += n;
}

// *******************************************************************
void CGappedSequence::RemoveGap(size_t pos)
{
 	// Look for the place to remove the gap
	size_t x = 1;

	while(x < dps_size)
	{
		if(dps[2 * x] >= pos)
			x = 2 * x;
		else
		{
			pos -= dps[2 * x];
			x = 2 * x + 1;
		}
	}

	// Decrement the no. of gaps
	--n_gaps[x - dps_size];

	// Update DSP
	for(; x; x /= 2)
		--dps[x];

    --gapped_size;
}

// *******************************************************************
void CGappedSequence::RemoveGaps(size_t pos, uint32_t n)
{
 	// Look for the place to remove the gap
	size_t x = 1;

	while(x < dps_size)
	{
		if(dps[2 * x] >= pos)
			x = 2 * x;
		else
		{
			pos -= dps[2 * x];
			x = 2 * x + 1;
		}
	}

	// Decrement the no. of gaps
	n_gaps[x - dps_size] -= n;

	// Update DSP
	for(; x; x /= 2)
		dps[x] -= n;

    gapped_size -= n;
}

// *******************************************************************
symbol_t CGappedSequence::GetSymbol(size_t pos)
{
 	// Look for the place to remove the gap
	size_t x = 1;

	while(x < dps_size)
	{
		if(dps[2 * x] >= pos)
			x = 2 * x;
		else
		{
			pos -= dps[2 * x];
			x = 2 * x + 1;
		}
	}

	if(pos == dps[x])
		return symbols[x - dps_size + 1];
	else
		return GAP;
}
