/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#include "../core/sequence.h"
#include "../utils/utils.h"

#include <algorithm>
#include <assert.h>

using namespace std;

char CSequence::mapping_table[25] = "ARNDCQEGHILKMFPSTWYVBZX*";
char CGappedSequence::mapping_table[25] = "ARNDCQEGHILKMFPSTWYVBZX*";


// *******************************************************************
CSequence::CSequence(const string& _id, const string& seq) 
	: 
	sequence_no(-1), 
	length((uint32_t)seq.length()), 
	id(_id),
	data(seq.length()), 
	uppercase(seq.length())
{
	for(uint32_t i = 0; i < length; ++i)
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
	delete[] hist;
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
	uint32_t bv_len = (uint32_t) ((data.size() + bv_size - 1) / bv_size);

	/*
	bit_masks.clear();
	bit_masks.resize(NO_SYMBOLS);

	for (size_t i = 0; i < NO_SYMBOLS; ++i)
		bit_masks[i].resize(bv_len, (bit_vec_t) 0);
	*/

	bit_masks.resize(bv_len, (uint32_t) NO_SYMBOLS, (bit_vec_t)0);

	for(size_t i = 0; i < length; ++i)
		if(data[i] >= 0 && data[i] < NO_VALID_AMINOACIDS)
			bit_masks[data[i]][i / bv_size] |= ((bit_vec_t) 1) << (i % bv_size);
}

// *******************************************************************
void CSequence::ReleaseBitMasks() 
{
	bit_masks.clear();
}

// *******************************************************************
void CSequence::PrepareHistogram()
{
	hist = new uint16_t[NO_AMINOACIDS];

/*	int p1 = length / 3;
	int p2 = 2 * length / 3;*/
	
//	fill_n(i_hist, NO_AMINOACIDS, 0u);
	fill_n(hist, NO_AMINOACIDS, 0u);
//	fill_n(hist, NO_SYMBOLS, 0u);

/*	fill_n(hist0, NO_SYMBOLS, 0u);
	fill_n(hist1, NO_SYMBOLS, 0u);
	fill_n(hist2, NO_SYMBOLS, 0u);
	fill_n(hist01, NO_SYMBOLS, 0u);
	fill_n(hist12, NO_SYMBOLS, 0u);
	fill_n(hist012, NO_SYMBOLS, 0u);*/

	for (size_t i = 0; i < length; ++i)
	{
//		++i_hist[data[i]];
		++hist[data[i]];

/*		if (i < p1)
		{
			++hist0[data[i]];
			++hist01[data[i]];
			++hist012[data[i]];
		}
		else if (i < p2)
		{
			++hist1[data[i]];
			++hist01[data[i]];
			++hist012[data[i]];
		}
		else
		{
			++hist2[data[i]];
			++hist12[data[i]];
			++hist012[data[i]];
		}*/
	}
}


// *******************************************************************
//
// *******************************************************************

// *******************************************************************
CGappedSequence::CGappedSequence(CSequence&& _sequence)
{
	id = std::move(_sequence.id);
	symbols = std::move(_sequence.data);
	uppercase = std::move(_sequence.uppercase);
	_sequence.bit_masks.clear();
	
	size = symbols.size();
	gapped_size = size;
	n_gaps.resize(size + 1, 0);

	InitialiseDPS();
}

// *******************************************************************
CGappedSequence::CGappedSequence(const CGappedSequence& _gapped_sequence)
{
	id = _gapped_sequence.id;

	size = _gapped_sequence.size;
	gapped_size = _gapped_sequence.gapped_size;
	dps_size = _gapped_sequence.dps_size;
	dps_size_div2 = _gapped_sequence.dps_size_div2;

	symbols = _gapped_sequence.symbols;
	n_gaps = _gapped_sequence.n_gaps;
	dps = _gapped_sequence.dps;
	uppercase = _gapped_sequence.uppercase;
}

// *******************************************************************
CGappedSequence::CGappedSequence(CGappedSequence&& _gapped_sequence) noexcept
{
	id = move(_gapped_sequence.id);

	size = _gapped_sequence.size;
	gapped_size = _gapped_sequence.gapped_size;
	dps_size = _gapped_sequence.dps_size;
	dps_size_div2 = _gapped_sequence.dps_size_div2;

	symbols = move(_gapped_sequence.symbols);
	n_gaps = move(_gapped_sequence.n_gaps);
	dps = move(_gapped_sequence.dps);
	uppercase = move(_gapped_sequence.uppercase);
}

// *******************************************************************
CGappedSequence::~CGappedSequence()
{
}

// *******************************************************************
bool CGappedSequence::operator==(const CGappedSequence& gs)
{
	return id == gs.id &&
		gapped_size == gs.gapped_size &&
		symbols == gs.symbols &&
		n_gaps == gs.n_gaps;
}

// *******************************************************************
void CGappedSequence::RecalculateDPS()
{
	// Last level - from n_gaps
	size_t n = dps_size_div2;

	for (size_t i = 0; i < size / 2; ++i)
		dps[n + i] = n_gaps[2 * i] + n_gaps[2 * i + 1] + 2;

	if (size & 1)
		dps[n + size / 2] = n_gaps[size - 1] + n_gaps[size] + 2;
	else
		dps[n + size / 2] = n_gaps[size] + 1;

	// Last but one level - from dps, but take case of size of dps last level
	n = dps_size_div2 / 2;
	for (size_t i = 0; i < size / 4; ++i)
		dps[n + i] = dps[2 * n + 2 * i] + dps[2 * n + 2 * i + 1];

	if ((size / 2) & 1)
		dps[n + size / 4] = dps[2 * n + 2 * (size / 4)] + dps[2 * n + 2 * (size / 4) + 1];
	else
		dps[n + size / 4] = dps[2 * n + 2 * (size / 4)];

	// Remaining levels
	for (n = dps_size_div2 / 4; n; n /= 2)
		for (size_t i = 0; i < n; ++i)
			dps[n + i] = dps[2 * n + 2 * i] + dps[2 * n + 2 * i + 1];
}

// *******************************************************************
void CGappedSequence::InitialiseDPS()
{
	// Round size of n_gap to the nearest not smaller power of 2
	dps_size = size + 1;
	if ((dps_size & (dps_size - 1)) != 0)
	{
		while ((dps_size & (dps_size - 1)) != 0)
			dps_size &= dps_size - 1;
		dps_size <<= 1;
	}

	dps_size_div2 = dps_size / 2;
//	dps.resize(dps_size_div2 + size + 1, 0);
	dps.resize(dps_size_div2 + size / 2 + 1, 0);

	RecalculateDPS();
}

// *******************************************************************
string CGappedSequence::Decode()
{
	string s;

	s.reserve(gapped_size);

	// Starting gaps
	s.append(n_gaps[0], '-');

	for (int i = 1; i <= size; ++i)
	{
		char symbol = mapping_table[symbols[i]];

		if (!uppercase[i - 1])
			symbol += 32;		// change to lowercase

		s.push_back(symbol);

		s.append(n_gaps[i], '-');
	}

	return s;
}

// *******************************************************************
void CGappedSequence::DecodeRaw(symbol_t* seq)
{
	uint32_t seq_pos = 1;

	seq[0] = GUARD;

	// Starting gaps
	for (int j = 0; j < n_gaps[0]; ++j)
		seq[seq_pos++] = GAP;
	for (int i = 1; i <= size; ++i)
	{
		seq[seq_pos++] = symbols[i];

		for (int j = 0; j < n_gaps[i]; ++j)
			seq[seq_pos++] = GAP;
	}
}

// *******************************************************************
void CGappedSequence::InsertGap(size_t pos)
{
	// Look for the place to insert the gap
	size_t x = 1;

	++dps[1];

	while (x < dps_size_div2)
	{
		x *= 2;
		if (dps[x] < pos)
			pos -= dps[x++];
		++dps[x];

		if (x >= dps_size_div2)
			break;

		x *= 2;
		if (dps[x] < pos)
			pos -= dps[x++];
		++dps[x];

		if (x >= dps_size_div2)
			break;

		x *= 2;
		if (dps[x] < pos)
			pos -= dps[x++];
		++dps[x];

		if (x >= dps_size_div2)
			break;

		x *= 2;
		if (dps[x] < pos)
			pos -= dps[x++];
		++dps[x];
	}

	x = 2 * x - dps_size;

	if (n_gaps[x] + 1 < pos)
		++x;

	++n_gaps[x];
	++gapped_size;
}

// *******************************************************************
void CGappedSequence::InsertGaps(size_t pos, uint32_t n)
{
	// Look for the place to insert the gap
	size_t x = 1;

	dps[1] += n;

	while (x < dps_size_div2)
	{
		x *= 2;
		if (dps[x] < pos)
			pos -= dps[x++];
		dps[x] += n;

		if (x >= dps_size_div2)
			break;

		x *= 2;
		if (dps[x] < pos)
			pos -= dps[x++];
		dps[x] += n;

		if (x >= dps_size_div2)
			break;

		x *= 2;
		if (dps[x] < pos)
			pos -= dps[x++];
		dps[x] += n;

		if (x >= dps_size_div2)
			break;

		x *= 2;
		if (dps[x] < pos)
			pos -= dps[x++];
		dps[x] += n;
	}

	x = 2 * x - dps_size;

	if (n_gaps[x] + 1 < pos)
		++x;

	n_gaps[x] += n;
	gapped_size += n;
}

// *******************************************************************
void CGappedSequence::InsertGapsVector(const vector<pair<uint32_t, uint32_t>>& v_gaps)
{
	uint32_t c_pos = 0;
	uint32_t data_idx = 0;

	for (auto& x : v_gaps)
	{
		while (x.first > c_pos + 1 + n_gaps[data_idx])
		{
			c_pos += n_gaps[data_idx] + 1;
			++data_idx;
		}

		if (data_idx == n_gaps.size())
			--data_idx;

		n_gaps[data_idx] += x.second;

		gapped_size += x.second;
	}

	RecalculateDPS();
}

// *******************************************************************
uint32_t CGappedSequence::NoSymbols()
{
	return (uint32_t) symbols.size();
}

// *******************************************************************
void CGappedSequence::RemoveGap(size_t pos)
{
	// Look for the place to remove the gap
	size_t x = 1;

	while (x < dps_size_div2)
	{
		if (dps[2 * x] >= pos)
			x = 2 * x;
		else
		{
			pos -= dps[2 * x];
			x = 2 * x + 1;
		}
	}

	x = 2 * x - dps_size;
	if (n_gaps[x] + 1 < pos)
		++x;

	// Decrement the no. of gaps
	--n_gaps[x];

	x += dps_size;
	x /= 2;

	// Update DSP
	for (; x; x /= 2)
		--dps[x];

	--gapped_size;
}

// *******************************************************************
void CGappedSequence::RemoveGaps(size_t pos, uint32_t n)
{
	// Look for the place to remove the gap
	size_t x = 1;

	while (x < dps_size_div2)
	{
		if (dps[2 * x] >= pos)
			x = 2 * x;
		else
		{
			pos -= dps[2 * x];
			x = 2 * x + 1;
		}
	}

	x = 2 * x - dps_size;
	if (n_gaps[x] + 1 < pos)
		++x;

	// Decrement the no. of gaps
	n_gaps[x] -= n;

	x += dps_size;
	x /= 2;

	// Update DSP
	for (; x; x /= 2)
		dps[x] -= n;

	gapped_size -= n;
}

// *******************************************************************
symbol_t CGappedSequence::GetSymbol(size_t pos)
{
	// Look for the place to remove the gap
	size_t x = 1;

	while (x < dps_size_div2)
	{
		if (dps[2 * x] >= pos)
			x = 2 * x;
		else
		{
			pos -= dps[2 * x];
			x = 2 * x + 1;
		}
	}

	x = 2 * x - dps_size;
	if (n_gaps[x] + 1 < pos)
	{
		pos -= n_gaps[x] + 1;
		++x;
	}

	if (pos == n_gaps[x] + 1)
		return symbols[x + 1];
	else
		return GAP;
}

// *******************************************************************
void CGappedSequence::Clear()
{
	clear_vector(symbols);
	clear_vector(uppercase);
	clear_vector(n_gaps);
	clear_vector(dps);
}

// *******************************************************************
void CGappedSequence::ClearDPS()
{
	clear_vector(dps);
}
