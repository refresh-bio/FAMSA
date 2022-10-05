/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

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
CSequence::CSequence(const string& _id, const string& seq, int sequence_no, memory_monotonic_safe* mma)
	: 
	length(0), 
	original_no(sequence_no),
	sequence_no(sequence_no),
	id(_id),
	mma(mma)
{
	// omit gaps
	for (const char c : seq) {
		if (c != '-') {
			++length;
		}
	}
	
	uppercase.resize(length);
	data_size = length;

	if (length)
	{
		if (mma)
			data = (symbol_t*)mma->allocate(data_size + 1);
		else
			data = new symbol_t[data_size + 1];
	}
	else
		data = nullptr;

	p_bit_masks = nullptr;
	p_bv_len = 0;

	for(uint32_t i = 0, i_out = 0; i < seq.length(); ++i)
	{
		char c = seq[i];
		
		if (c == '-') {
			continue;
		}
		
		if(c > 'Z') {
			c -= 32;
			uppercase[i_out] = false;
		}
		else {
			uppercase[i_out] = true;
		}

		char *q = find(mapping_table, mapping_table+25, c);
		if (q == mapping_table + 25) {	
			extra_symbols.emplace_back(i_out, c); // save non-standard symbol
			data[i_out] = (symbol_t)UNKNOWN_SYMBOL;
		}
		else {
			data[i_out] = (symbol_t)(q - mapping_table);
		}

		++i_out;
	}
}

// *******************************************************************
CSequence::CSequence(CSequence&& x) noexcept 
	:
	original_no(x.original_no)
{
	sequence_no = move(x.sequence_no);
	length = move(x.length);
	id = move(x.id);

	data = x.data;
	x.data = nullptr;
	data_size = x.data_size;

	mma = x.mma;
	x.mma = nullptr;

	uppercase = move(x.uppercase);
	extra_symbols = move(x.extra_symbols);

	p_bit_masks = move(x.p_bit_masks);
	x.p_bit_masks = nullptr;
	p_bv_len = x.p_bv_len;
}


// *******************************************************************
/*
CSequence& CSequence::operator=(CSequence&& x) noexcept
{
	this->sequence_no = move(x.sequence_no);
	this->length = move(x.length);
	this->id = move(x.id);

	delete_arr_ptr(this->data);
	this->data = x.data;
	x.data = nullptr;
	data_size = x.data_size;
	this->mma = x.mma;
	x.mma = nullptr;

	this->uppercase = move(x.uppercase);

	this->p_bit_masks = x.p_bit_masks;
	x.p_bit_masks = nullptr;
	this->p_bv_len = x.p_bv_len;

	return *this;
}
*/

// *******************************************************************
CSequence::~CSequence()
{
	delete_arr_ptr(p_bit_masks);

	if (mma)
		mma->deallocate(data);
	else
		delete_arr_ptr(data);
}

// *******************************************************************
/*
string CSequence::DecodeSequence()
{
	string s;

	s.reserve(length + 1);
	for(uint32_t i = 0; i < length; ++i)
		if(data[i] == GUARD)
			continue;
		else if(data[i] == GAP)
			s += "-";
		else
				s += mapping_table[data[i]];

	return s;
}
*/

// *******************************************************************
void CSequence::DataResize(uint32_t new_size, symbol_t new_symbol)
{
	assert(new_size >= length);

	symbol_t* ptr;
	
	if (mma)
		ptr = (symbol_t*)mma->allocate(new_size + 1);
	else
		ptr = new symbol_t[new_size + 1];

	copy_n(data, min(data_size, new_size), ptr);

	if(new_size > data_size)
		fill(ptr + data_size, ptr + new_size, new_symbol);
	swap(data, ptr);

	data_size = new_size;

	if (mma)
		mma->deallocate(ptr);
	else
		delete[] ptr;
}

// *******************************************************************
// Compute Bit-mask Vectors for bit-parallel computation of LCS for the sequences
void CSequence::ComputeBitMasks()
{
	p_bv_len = (uint32_t) ((data_size + bv_size - 1) / bv_size);

	delete_arr_ptr(p_bit_masks);
	p_bit_masks = new bit_vec_t[p_bv_len * NO_SYMBOLS];
	fill_n(p_bit_masks, p_bv_len * NO_SYMBOLS, (bit_vec_t)0);

	for(size_t i = 0; i < length; ++i)
		if(data[i] >= 0 && data[i] < NO_VALID_AMINOACIDS)
			p_bit_masks[data[i] * p_bv_len + i / bv_size] |= ((bit_vec_t) 1) << (i % bv_size);
}

// *******************************************************************
void CSequence::ReleaseBitMasks() 
{
	delete_arr_ptr(p_bit_masks);
}

// *******************************************************************
//
// *******************************************************************

// *******************************************************************
CGappedSequence::CGappedSequence(const string& _id, const string& seq, int seq_no, memory_monotonic_safe* mma) 
	: mma(mma), sequence_no(seq_no), id(_id)
{
	gapped_size = seq.size();
	size = 0;
	for (const auto &ch : seq) {
		if (ch!='-') {
			++size;
		}
	}

	symbols_size = size;
	uppercase.resize(symbols_size);
	n_gaps.resize(symbols_size + 1, 0);

	if (size) {
		if (mma) { 
			symbols=(symbol_t*)mma->allocate(symbols_size + 1);}
		else {
			symbols=new symbol_t[symbols_size +1];}
	}
	else {
		symbols=nullptr;
	}
  
	int is = 0;
	for (int i = 0; i < (int)gapped_size; ++i) {
		char c = seq[i];
		
		if (c == '-') {
			++n_gaps[is];
		}
		else {
			
			if (c > 'Z') {
				c -= 32;
				uppercase[is] = false;
			}
			else {
				uppercase[is] = true;
			}
			char* q = find(mapping_table, mapping_table + 25, c);
			if (q == mapping_table + 25) {
				extra_symbols.emplace_back(is, c); // save non-standard symbol
				symbols[is] = (symbol_t)UNKNOWN_SYMBOL;
			}
			else {
				symbols[is] = (symbol_t)(q - mapping_table);
			}

			++is;

		}
	}

	InitialiseDPS();
}
//*******************88******************************************************************************************************************************************************************************************************************************************************

CGappedSequence::CGappedSequence(CSequence&& _sequence) :
	mma(_sequence.mma),
	symbols(std::move(_sequence.data)),
	size(_sequence.data_size),
	original_no(_sequence.original_no),
	sequence_no(_sequence.sequence_no),
	id(std::move(_sequence.id)),
	uppercase(_sequence.uppercase),
	extra_symbols(_sequence.extra_symbols)
{
	
	_sequence.data = nullptr;
	_sequence.mma = nullptr;
	delete_arr_ptr(_sequence.p_bit_masks);
	
	symbols_size = size;
	gapped_size = size;
	n_gaps.resize(size + 1, 0);

	InitialiseDPS();
}

// *******************************************************************
CGappedSequence::CGappedSequence(const CGappedSequence& _gapped_sequence) :
	original_no(_gapped_sequence.original_no)
{
	id = _gapped_sequence.id;
	sequence_no = _gapped_sequence.sequence_no;

	size = _gapped_sequence.size;
	symbols_size = _gapped_sequence.symbols_size;
	gapped_size = _gapped_sequence.gapped_size;
	dps_size = _gapped_sequence.dps_size;
	dps_size_div2 = _gapped_sequence.dps_size_div2;
	mma = _gapped_sequence.mma;


	if (mma)
		symbols = (symbol_t*)mma->allocate(symbols_size + 1);
	else
		symbols = new symbol_t[symbols_size + 1];

	
	copy_n(_gapped_sequence.symbols, symbols_size, symbols);

	n_gaps = _gapped_sequence.n_gaps;
	dps = _gapped_sequence.dps;
	uppercase = _gapped_sequence.uppercase;
	extra_symbols = _gapped_sequence.extra_symbols;
}

// *******************************************************************
CGappedSequence::CGappedSequence(CGappedSequence&& _gapped_sequence) noexcept 
	:
	original_no(_gapped_sequence.original_no)
{
	id = move(_gapped_sequence.id);
	sequence_no = _gapped_sequence.sequence_no;

	size = _gapped_sequence.size;
	symbols_size = _gapped_sequence.symbols_size;
	gapped_size = _gapped_sequence.gapped_size;
	dps_size = _gapped_sequence.dps_size;
	dps_size_div2 = _gapped_sequence.dps_size_div2;

	symbols = _gapped_sequence.symbols;
	mma = _gapped_sequence.mma;

	_gapped_sequence.symbols = nullptr;
	_gapped_sequence.mma = nullptr;

	n_gaps = move(_gapped_sequence.n_gaps);
	dps = move(_gapped_sequence.dps);
	uppercase = move(_gapped_sequence.uppercase);
	extra_symbols = move(_gapped_sequence.extra_symbols);
}

// *******************************************************************
CGappedSequence::~CGappedSequence()
{
	if (mma)
		mma->deallocate(symbols);
	else
		delete_arr_ptr(symbols);
}

// *******************************************************************
bool CGappedSequence::operator==(const CGappedSequence& gs) const
{
	bool r = 
		id == gs.id &&
		gapped_size == gs.gapped_size &&
		size == gs.size &&
		symbols_size == gs.symbols_size &&
		n_gaps == gs.n_gaps;

	if (!r)
		return false;

	return equal(symbols, symbols + symbols_size, gs.symbols);
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
	dps.resize(dps_size_div2 + size / 2 + 1, 0);

	RecalculateDPS();
}

// *******************************************************************
string CGappedSequence::Decode()
{
	string s;

	// decode symbols
	for (int i = 1; i <= (int)size; ++i) {
		symbols[i] = mapping_table[symbols[i]];
	}

	// replace extra characters
	for (auto it = extra_symbols.begin(); it < extra_symbols.end(); ++it) {
		symbols[it->first + 1] = it->second;
	}

	s.reserve(gapped_size);

	// Starting gaps
	s.append(n_gaps[0], '-');

	for (int i = 1; i <= (int) size; ++i)
	{	
		char symbol = symbols[i];
		
		if (!uppercase[i - 1])
			symbol += 32;		// change to lowercase

		s.push_back(symbol);

		s.append(n_gaps[i], '-');
	}

	return s;
}

// *******************************************************************
/*
void CGappedSequence::DecodeRaw(symbol_t* seq)
{
	uint32_t seq_pos = 1;

	seq[0] = GUARD;

	// Starting gaps
	for (uint32_t j = 0; j < n_gaps[0]; ++j)
		seq[seq_pos++] = GAP;
	for (uint32_t i = 1; i <= (uint32_t) size; ++i)
	{
		seq[seq_pos++] = symbols[i];

		for (uint32_t j = 0; j < n_gaps[i]; ++j)
			seq[seq_pos++] = GAP;
	}
}
*/
// *******************************************************************
void CGappedSequence::InsertGap(uint32_t pos)
{
	// Look for the place to insert the gap
	uint32_t x = 1;

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
void CGappedSequence::InsertGaps(uint32_t pos, uint32_t n)
{
	// Look for the place to insert the gap
	uint32_t x = 1;

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
	return (uint32_t) symbols_size;
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
	if (mma)
		mma->deallocate(symbols);
	else
		delete_arr_ptr(symbols);

	clear_vector(uppercase);
	clear_vector(n_gaps);
	clear_vector(dps);
}

// *******************************************************************
void CGappedSequence::ClearDPS()
{
	clear_vector(dps);
}

// *******************************************************************
void CGappedSequence::InsertFront(symbol_t new_symbol)
{
	// It is assumed that there is a room in memory for extension by 1 symbol
	copy_backward(symbols, symbols + symbols_size, symbols + symbols_size + 1);
	++symbols_size;
	*symbols = new_symbol;
}
