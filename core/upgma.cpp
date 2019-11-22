/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys


Note: The file contains code for the UPGMA method (of high memory consumption)
These functions was borrowed and adpoted from MUSCLE 3.8.1551 by Robert Edgar
The time and memory consumption is O(k^2)
*/

#include "../core/msa.h"
#include <thread>

#include "../libs/vectorclass.h"

using namespace std;

#define	AVG(x, y)	(((x) + (y))/2)

// *******************************************************************
inline uint64_t CFAMSA::UPGMA_TriangleSubscript(uint64_t uIndex1, uint64_t uIndex2)
{
	if (uIndex1 >= uIndex2)
		return uIndex2 + (uIndex1 * (uIndex1 - 1)) / 2;
	else
		return uIndex1 + (uIndex2 * (uIndex2 - 1)) / 2;
}

// *******************************************************************
void CFAMSA::UPGMA_CalculateDistances(UPGMA_dist_t *dist_matrix)
{
	int next;
	int n_seq = sequences.size();

	double indel_exp = params.indel_exp;

	// Temporarily resize the sequences by adding unknown symbols (just to simplify the pairwise comparison)
	size_t max_seq_len = max_element(sequences.begin(), sequences.end(), [](const CSequence &x, CSequence &y) {return x.length < y.length; })->length;

	for (int i = 1; i < n_seq; ++i)
		sequences[i].data.resize(max_seq_len, UNKNOWN_SYMBOL);

	CUPGMAQueue slq(&sequences, n_seq, dist_matrix);
	vector<thread *> workers(n_threads, nullptr);

	uint32_t computed_prof = 0;
	mutex mtx;

	// Calculation of similarities is made in working threads
	for (uint32_t i = 0; i < n_threads; ++i)
		workers[i] = new thread([&] {
		CLCSBP lcsbp(instruction_set);
		int row_id;
		vector<CSequence> *sequences;
		UPGMA_dist_t *dist_row;
		uint32_t lcs_lens[4];

		while (slq.GetTask(row_id, sequences, dist_row))
		{
			for (int j = 0; j < row_id / 4; ++j)
			{
				lcsbp.GetLCSBP(&(*sequences)[row_id], &(*sequences)[j * 4 + 0], &(*sequences)[j * 4 + 1], &(*sequences)[j * 4 + 2], &(*sequences)[j * 4 + 3],
					lcs_lens[0], lcs_lens[1], lcs_lens[2], lcs_lens[3]);
				for (int k = 0; k < 4; ++k)
				{
					double indel = (*sequences)[row_id].length + (*sequences)[j * 4 + k].length - 2 * lcs_lens[k];
					dist_row[j * 4 + k] = 1.0 / (lcs_lens[k] / pow(indel, indel_exp));
				}
			}

			if (row_id / 4 * 4 < row_id)
			{
				lcsbp.GetLCSBP(&(*sequences)[row_id],
					(row_id / 4 * 4 + 0 < row_id) ? &(*sequences)[row_id / 4 * 4 + 0] : nullptr,
					(row_id / 4 * 4 + 1 < row_id) ? &(*sequences)[row_id / 4 * 4 + 1] : nullptr,
					(row_id / 4 * 4 + 2 < row_id) ? &(*sequences)[row_id / 4 * 4 + 2] : nullptr,
					(row_id / 4 * 4 + 3 < row_id) ? &(*sequences)[row_id / 4 * 4 + 3] : nullptr,
					lcs_lens[0], lcs_lens[1], lcs_lens[2], lcs_lens[3]);
				for (int k = 0; k < 4 && row_id / 4 * 4 + k < row_id; ++k)
				{
					double indel = (*sequences)[row_id].length + (*sequences)[row_id / 4 * 4 + k].length - 2 * lcs_lens[k];
					dist_row[row_id / 4 * 4 + k] = 1.0 / (lcs_lens[k] / pow(indel, indel_exp));
				}
			}
		}
	});

	for (auto p : workers)
	{
		p->join();
		delete p;
	}
	workers.clear();

	// Bring the sequences to the valid length
	for (int i = 1; i < n_seq; ++i)
		sequences[i].data.resize(sequences[i].length, UNKNOWN_SYMBOL);
}

// *******************************************************************
void CFAMSA::UPGMA()
{
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
	g_uLeafCount = sequences.size();

	g_uTriangleSize = (g_uLeafCount*(g_uLeafCount - 1)) / 2;
	g_uInternalNodeCount = g_uLeafCount - 1;

	g_Dist = new UPGMA_dist_t[g_uTriangleSize];

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

	// Calculate complete (triangular) matrix of distances for UPGMA algorithm
	UPGMA_CalculateDistances(g_Dist);

	if (params.distance_matrix_out_file.length() > 0) {
		ExportDistanceMatrix(g_Dist, sequences.size(), params.distance_matrix_out_file);
	}
		
	// Compute initial NxN triangular distance matrix.
	// Store minimum distance for each full (not triangular) row.
	// Loop from 1, not 0, because "row" is 0, 1 ... i-1,
	// so nothing to do when i=0.
	for (uint64_t i = 1; i < g_uLeafCount; ++i)
	{
		UPGMA_dist_t *Row = g_Dist + UPGMA_TriangleSubscript(i, 0);
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

			const uint64_t vL = UPGMA_TriangleSubscript(Lmin, j);
			const uint64_t vR = UPGMA_TriangleSubscript(Rmin, j);
			const UPGMA_dist_t dL = g_Dist[vL];
			const UPGMA_dist_t dR = g_Dist[vR];
			UPGMA_dist_t dtNewDist;

			dtNewDist = AVG(dL, dR);

			if (g_uNearestNeighbor[j] == Rmin)
				g_uNearestNeighbor[j] = Lmin;

			g_Dist[vL] = dtNewDist;
			if (dtNewDist < dtNewMinDist)
			{
				dtNewMinDist = dtNewDist;
				uNewNearestNeighbor = j;
			}
		}

		const uint64_t v = UPGMA_TriangleSubscript(Lmin, Rmin);
		const UPGMA_dist_t dLR = g_Dist[v];
		const UPGMA_dist_t dHeightNew = dLR / 2;
		const uint64_t uLeft = g_uNodeIndex[Lmin];
		const uint64_t uRight = g_uNodeIndex[Rmin];
		const UPGMA_dist_t HeightLeft  = uLeft < g_uLeafCount ? 0 : g_Height[uLeft - g_uLeafCount];
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

	int n_seq = sequences.size();
	for (int i = 0; i < n_seq - 1; ++i)
		guide_tree.push_back(make_pair(g_uLeft[i], g_uRight[i]));

	delete[] g_Dist;

	delete[] g_uNodeIndex;
	delete[] g_uNearestNeighbor;
	delete[] g_MinDist;
	delete[] g_Height;

	delete[] g_uLeft;
	delete[] g_uRight;
	delete[] g_LeftLength;
	delete[] g_RightLength;
}
