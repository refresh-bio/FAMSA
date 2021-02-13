/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

Note: The file contains code for the UPGMA method (of high memory consumption)
These functions was borrowed and adpoted from MUSCLE 3.8.1551 by Robert Edgar
The time and memory consumption is O(k^2)

*/
#include "UPGMA.h"
#include "AbstractTreeGenerator.hpp"

#include "lcsbp.h"

#include <thread>
#include <algorithm>
#include <fstream>

using namespace std;


// general average
template <class T, bool is_modified>
struct Average {
	T operator()(T x, T y) const { return (x + y) * 0.5f; }
};

// modified average (MAFFT)
template <class T>
struct Average<T, true> {
	T operator()(T x, T y) const { return 0.05f * (x + y) + 0.9f * std::min(x, y); }
};


// *******************************************************************
void UPGMA::run(std::vector<CSequence>& sequences, tree_structure& tree) {
	UPGMA_dist_t* distances = TriangleMatrix::allocate<UPGMA_dist_t>(sequences.size());
	computeDistances(sequences, distances);

	if (is_modified) {
		computeTree<true>(distances, sequences.size(), tree);
	}
	else {
		computeTree<false>(distances, sequences.size(), tree);
	}

	delete[] distances;
}

// *******************************************************************
void UPGMA::runPartial(std::vector<CSequence*>& sequences, tree_structure& tree) {
	UPGMA_dist_t* distances = TriangleMatrix::allocate<UPGMA_dist_t>(sequences.size());
	CLCSBP lcsbp(instruction_set);
	calculateSimilarityMatrix<CSequence*, UPGMA_dist_t, Measure::DistanceReciprocal>(sequences.data(), sequences.size(), distances, lcsbp);

	if (is_modified) {
		computeTree<true>(distances, sequences.size(), tree);
	}
	else {
		computeTree<false>(distances, sequences.size(), tree);
	}

	delete[] distances;
}


// *******************************************************************
void UPGMA::computeDistances(std::vector<CSequence>& sequences, UPGMA_dist_t *dist_matrix)
{
	size_t n_seq = sequences.size();

	CUPGMAQueue slq(&sequences, n_seq, dist_matrix);
	vector<thread *> workers(n_threads, nullptr);

	mutex mtx;

	// Calculation of similarities is made in working threads
	for (uint32_t i = 0; i < n_threads; ++i)
		workers[i] = new thread([&] {
		CLCSBP lcsbp(instruction_set);
		int row_id;
		vector<CSequence> *sequences;
		UPGMA_dist_t *dist_row;

		while (slq.GetTask(row_id, sequences, dist_row))
		{
			calculateSimilarityVector<CSequence, UPGMA_dist_t, Measure::DistanceReciprocal>(
				(*sequences)[row_id],
				sequences->data(),
				row_id,
				dist_row,
				lcsbp);
		}
	});

	for (auto p : workers)
	{
		p->join();
		delete p;
	}
	workers.clear();
}

// *******************************************************************
template <bool MODIFIED>
void UPGMA::computeTree(UPGMA_dist_t* distances, size_t n_seq, tree_structure& tree)
{
	Average<UPGMA_dist_t, MODIFIED> average;

	uint64_t g_uLeafCount;
	uint64_t g_uTriangleSize;
	uint64_t g_uInternalNodeCount;
	uint64_t g_uInternalNodeIndex;

	// Triangular distance matrix is g_Dist, which is allocated
	// as a one-dimensional vector of length g_uTriangleSize.
	// TriangleSubscript(i,j) maps row,column=i,j to the subscript
	// into this vector.
	// Row / column coordinates are a bit messy.
	// Initially they are leaf indexes 0..N-1.
	// But each time we create a new node (=new cluster, new subtree),
	// we re-use one of the two rows that become available (the children
	// of the new node). This saves memory.
	// We keep track of this through the g_uNodeIndex vector.
	UPGMA_dist_t *g_Dist;

	// Distance to nearest neighbor in row i of distance matrix.
	// Subscript is distance matrix row.
	UPGMA_dist_t *g_MinDist;

	// Nearest neighbor to row i of distance matrix.
	// Subscript is distance matrix row.
	uint64_t *g_uNearestNeighbor;

	// Node index of row i in distance matrix.
	// Node indexes are 0..N-1 for leaves, N..2N-2 for internal nodes.
	// Subscript is distance matrix row.
	uint64_t *g_uNodeIndex;

	// The following vectors are defined on internal nodes,
	// subscripts are internal node index 0..N-2.
	// For g_uLeft/Right, value is the node index 0 .. 2N-2
	// because a child can be internal or leaf.
	uint64_t *g_uLeft;
	uint64_t *g_uRight;
	UPGMA_dist_t *g_Height;
	UPGMA_dist_t *g_LeftLength;
	UPGMA_dist_t *g_RightLength;

	//	void UPGMA2(const DistCalc &DC, Tree &tree, LINKAGE Linkage)
	g_uLeafCount = n_seq;

	g_uTriangleSize = (g_uLeafCount*(g_uLeafCount - 1)) / 2;
	g_uInternalNodeCount = g_uLeafCount - 1;

	g_Dist = distances;

	g_uNodeIndex = new uint64_t[g_uLeafCount];
	g_uNearestNeighbor = new uint64_t[g_uLeafCount];
	g_MinDist = new UPGMA_dist_t[g_uLeafCount];

	g_uLeft = new uint64_t[g_uInternalNodeCount];
	g_uRight = new uint64_t[g_uInternalNodeCount];
	g_Height = new UPGMA_dist_t[g_uInternalNodeCount];
	g_LeftLength = new UPGMA_dist_t[g_uInternalNodeCount];
	g_RightLength = new UPGMA_dist_t[g_uInternalNodeCount];

	for (uint64_t i = 0; i < g_uLeafCount; ++i)
	{
		g_MinDist[i] = BIG_DIST;
		g_uNodeIndex[i] = i;
		g_uNearestNeighbor[i] = 0x7FFFFFFF;
	}

	for (uint64_t i = 0; i < g_uInternalNodeCount; ++i)
	{
		g_uLeft[i] = 0x7FFFFFFF;
		g_uRight[i] = 0x7FFFFFFF;
		g_LeftLength[i] = BIG_DIST;
		g_RightLength[i] = BIG_DIST;
		g_Height[i] = BIG_DIST;
	}

	// Compute initial NxN triangular distance matrix.
	// Store minimum distance for each full (not triangular) row.
	// Loop from 1, not 0, because "row" is 0, 1 ... i-1,
	// so nothing to do when i=0.
	for (uint64_t i = 1; i < g_uLeafCount; ++i)
	{
		const UPGMA_dist_t *Row = g_Dist + TriangleMatrix::access(i, 0);
		for (uint64_t j = 0; j < i; ++j)
		{
			const UPGMA_dist_t d = Row[j];
			if (d < g_MinDist[i])
			{
				g_MinDist[i] = d;
				g_uNearestNeighbor[i] = j;
			}
			if (d < g_MinDist[j])
			{
				g_MinDist[j] = d;
				g_uNearestNeighbor[j] = i;
			}
		}
	}

	for (g_uInternalNodeIndex = 0; g_uInternalNodeIndex < g_uLeafCount - 1; ++g_uInternalNodeIndex)
	{
		// Find nearest neighbors
		uint64_t Lmin = 0x7FFFFFFF;
		uint64_t Rmin = 0x7FFFFFFF;
		UPGMA_dist_t dtMinDist = BIG_DIST;
		for (uint64_t j = 0; j < g_uLeafCount; ++j)
		{
			if (0x7FFFFFFF == g_uNodeIndex[j])
				continue;

			UPGMA_dist_t d = g_MinDist[j];
			if (d < dtMinDist)
			{
				dtMinDist = d;
				Lmin = j;
				Rmin = g_uNearestNeighbor[j];
			}
		}

		// Compute distances to new node
		// New node overwrites row currently assigned to Lmin
		UPGMA_dist_t dtNewMinDist = BIG_DIST;
		uint64_t uNewNearestNeighbor = 0x7FFFFFFF;
		for (uint64_t j = 0; j < g_uLeafCount; ++j)
		{
			if (j == Lmin || j == Rmin)
				continue;
			if (0x7FFFFFFF == g_uNodeIndex[j])
				continue;

			const uint64_t vL = TriangleMatrix::access(Lmin, j);
			const uint64_t vR = TriangleMatrix::access(Rmin, j);
			const UPGMA_dist_t dL = g_Dist[vL];
			const UPGMA_dist_t dR = g_Dist[vR];
			UPGMA_dist_t dtNewDist;

			dtNewDist = average(dL, dR);

			if (g_uNearestNeighbor[j] == Rmin)
				g_uNearestNeighbor[j] = Lmin;

			g_Dist[vL] = dtNewDist;
			if (dtNewDist < dtNewMinDist)
			{
				dtNewMinDist = dtNewDist;
				uNewNearestNeighbor = j;
			}
		}

		const uint64_t v = TriangleMatrix::access(Lmin, Rmin);
		const UPGMA_dist_t dLR = g_Dist[v];
		const UPGMA_dist_t dHeightNew = dLR / 2;
		const uint64_t uLeft = g_uNodeIndex[Lmin];
		const uint64_t uRight = g_uNodeIndex[Rmin];
		const UPGMA_dist_t HeightLeft = uLeft < g_uLeafCount ? 0 : g_Height[uLeft - g_uLeafCount];
		const UPGMA_dist_t HeightRight = uRight < g_uLeafCount ? 0 : g_Height[uRight - g_uLeafCount];

		g_uLeft[g_uInternalNodeIndex] = uLeft;
		g_uRight[g_uInternalNodeIndex] = uRight;
		g_LeftLength[g_uInternalNodeIndex] = dHeightNew - HeightLeft;
		g_RightLength[g_uInternalNodeIndex] = dHeightNew - HeightRight;
		g_Height[g_uInternalNodeIndex] = dHeightNew;

		// Row for left child overwritten by row for new node
		g_uNodeIndex[Lmin] = g_uLeafCount + g_uInternalNodeIndex;
		g_uNearestNeighbor[Lmin] = uNewNearestNeighbor;
		g_MinDist[Lmin] = dtNewMinDist;

		// Delete row for right child
		g_uNodeIndex[Rmin] = 0x7FFFFFFF;
	}

	for (int i = 0; i < n_seq - 1; ++i)
		tree.push_back(make_pair(g_uLeft[i], g_uRight[i]));

	delete[] g_uNodeIndex;
	delete[] g_uNearestNeighbor;
	delete[] g_MinDist;
	delete[] g_Height;

	delete[] g_uLeft;
	delete[] g_uRight;
	delete[] g_LeftLength;
	delete[] g_RightLength;
}



// *******************************************************************
bool UPGMA::saveDistances(const std::string& file_name, std::vector<CSequence>& sequences) {

	// open output file
	ofstream file(file_name);
	if (!file) {
		return false;
	}

	// calculate distances first
	size_t g_uLeafCount = sequences.size();
	size_t g_uTriangleSize = (g_uLeafCount*(g_uLeafCount - 1)) / 2;
	UPGMA_dist_t *g_Dist = new UPGMA_dist_t[g_uTriangleSize];

	// Temporarily resize the sequences by adding unknown symbols (just to simplify the pairwise comparison)
	size_t max_seq_len =
		std::max_element(sequences.begin(), sequences.end(), [](const CSequence &x, CSequence &y) {return x.length < y.length; })->length;

	for (int i = 0; i < sequences.size(); ++i) {
		sequences[i].data.resize(max_seq_len, UNKNOWN_SYMBOL);
	}

	computeDistances(sequences, g_Dist);

	for (int i = 0; i < sequences.size(); ++i) {
		// Bring the sequences to the valid length
		sequences[i].data.resize(sequences[i].length, UNKNOWN_SYMBOL);
	
		// store distance in a file
		file << sequences[i].id << ", ";
			for (size_t j = 0; j < i; ++j) {
				const size_t id = TriangleMatrix::access(i, j);
				UPGMA_dist_t d = g_Dist[id];
				file << d << ", ";
			}
		file << std::endl;
	}

	delete[] g_Dist;
	return true;
}


// *******************************************************************
// Queue for UPGMA
// *******************************************************************
CUPGMAQueue::CUPGMAQueue(vector<CSequence> *_sequences, uint32_t _n_rows, UPGMA_dist_t *_dist_matrix)
{
	sequences = _sequences;
	n_rows = _n_rows;
	lowest_uncomputed_row = 0;
	eoq_flag = false;
	dist_matrix = _dist_matrix;
}

CUPGMAQueue::~CUPGMAQueue()
{
}

bool CUPGMAQueue::GetTask(int &row_id, vector<CSequence> *&_sequences, UPGMA_dist_t *&dist_row)
{
	unique_lock<mutex> lck(mtx);

	if (eoq_flag)
		return false;	// End of data in the profiles queue

	row_id = lowest_uncomputed_row++;

	if (lowest_uncomputed_row >= n_rows)
		eoq_flag = true;

	_sequences = sequences;

	dist_row = dist_matrix + TriangleMatrix::access(row_id, 0);

	return true;
}

