/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

Note: The file contains code for the UPGMA method (of high memory consumption)
These functions was borrowed and adpoted from MUSCLE 3.8.1551 by Robert Edgar
The time and memory consumption is O(k^2)

*/
#include "../core/guide_tree.h"
#include "../core/guide_tree.hpp"

#include "../libs/instrset.h"
#include "../core/queues.h"
#include "../core/NewickTree.h"

#include "../core/lcsbp.h"
#include "../core/lcsbp_classic.h"

#include "../core/log.h"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <random>

using namespace std;

#define	AVG(x, y)	(((x) + (y))/2)

// *******************************************************************
GuideTree::GuideTree(double indel_exp, size_t n_threads, int seed, int parttree_size, bool are_sorted)
	: indel_exp(indel_exp), n_threads(n_threads), seed(seed), parttree_size(parttree_size), are_sorted(are_sorted) {
	int x = instrset_detect();

	if (x >= 0 && x <= 8)
		instruction_set = (instruction_set_t)x;
	else if (x < 0)
		instruction_set = instruction_set_t::none;
	else
		instruction_set = instruction_set_t::avx2;
}

// *******************************************************************
void GuideTree::compute(std::vector<CSequence>& sequences, GT_method::Value method) {
	
	guide_tree.clear();
	guide_tree.resize(sequences.size(), std::make_pair<int, int>(-1, -1));

	// Temporarily resize the sequences by adding unknown symbols (just to simplify the pairwise comparison)
	size_t max_seq_len = are_sorted ? sequences[0].length
		: max_element(sequences.begin(), sequences.end(), [](const CSequence &x, CSequence &y) {return x.length < y.length; })->length;
	
	for (int i = 0; i < sequences.size(); ++i) {
		sequences[i].data.resize(max_seq_len, UNKNOWN_SYMBOL);
	}

	// perform the algorithm
	if (method == GT_method::SLINK) {
		computeSingleLinkage(sequences, guide_tree);
	}
	else if (method == GT_method::UPGMA) {
		UPGMA_dist_t* distances = TriangleMatrix::allocate<UPGMA_dist_t>(sequences.size());
		calculateDistances(sequences, distances);
		computeUPGMA(distances, sequences.size(), guide_tree);
		delete[] distances;
	}
	else if (method == GT_method::parttree_UPGMA || method == GT_method::parttree_SLINK) {
		computePartTree(sequences, method, guide_tree);
	}
#ifdef DEVELOPER_MODE
	else if (method == GT_method::chained)
		computeChained(sequences);
#endif

	// Bring the sequences to the valid length
	for (int i = 0; i < sequences.size(); ++i) {
		sequences[i].data.resize(sequences[i].length, UNKNOWN_SYMBOL);
	}
}

// *******************************************************************
bool GuideTree::loadNewick(
	const std::string& file_name,
	std::vector<CSequence>& sequences)
{
	// Load newick description
	ifstream newickFile;
	newickFile.open(file_name);
	if (!newickFile.good()) {
		return false;
	}

	std::stringstream ss;
	ss << newickFile.rdbuf();
	std::string description(ss.str());
	auto newend = std::remove_if(description.begin(), description.end(),
		[](char c)->bool { return c == '\r' || c == '\n';  });
	description.erase(newend, description.end());

	// Load guide tree
	NewickTree nw(false);
	nw.parse(sequences, description, guide_tree);

	return true;
}

// *******************************************************************
bool GuideTree::saveNewick(
	const std::string& file_name,
	const std::vector<CSequence>& sequences) const
{
	// store guide tree
	string description;
	NewickTree nw(false);
	nw.store(sequences, guide_tree, description);

	// Open file
	ofstream newickFile;
	newickFile.open(file_name);
	if (!newickFile.good()) {
		return false;
	}

	newickFile << description;

	return true;
}


// *******************************************************************
bool GuideTree::saveDistances(const std::string& file_name, std::vector<CSequence>& sequences) {
	
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
	size_t max_seq_len = sequences[0].length;
	for (int i = 1; i < sequences.size(); ++i) {
		sequences[i].data.resize(max_seq_len, UNKNOWN_SYMBOL);
	}

	calculateDistances(sequences, g_Dist);

	// Bring the sequences to the valid length
	for (int i = 1; i < sequences.size(); ++i) {
		sequences[i].data.resize(sequences[i].length, UNKNOWN_SYMBOL);
	}

	// store distances in a file
	for (size_t i = 0; i < sequences.size(); ++i) {
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
void GuideTree::computeSingleLinkage(std::vector<CSequence>& sequences, tree_structure& tree) {
	int next;
	int n_seq = sequences.size();

	int prefetch_offset = 64 * 2;
	int prefetch_offset_2nd = 128 * 2;

	vector<int> pi(n_seq + max(prefetch_offset, prefetch_offset_2nd), 0);
	vector<double> lambda(n_seq);
	vector<double> *sim_vector;

	CSingleLinkageQueue slq(&sequences, n_seq, n_threads * 8);
	vector<thread *> workers(n_threads, nullptr);

	mutex mtx;

	// Calculation of similarities is made in working threads
	for (uint32_t i = 0; i < n_threads; ++i)
		workers[i] = new thread([&] {
		CLCSBP lcsbp(instruction_set);
		int row_id;
		vector<CSequence> *sequences;
		vector<double> *sim_vector;
	
		while (slq.GetTask(row_id, sequences, sim_vector))
		{
			calculateSimilarityVector<CSequence, double>(
				&(*sequences)[row_id], 
				sequences->data(), 
				row_id,
				sim_vector->data(),
				lcsbp);
			
			slq.RegisterSolution(row_id);
		}
	});

	// Single linkage algorithm is here
	for (int i = 0; i < n_seq; ++i)
	{
		pi[i] = i;
		lambda[i] = -infty_double;

		if (i % (100) == 0) {
			LOG_DEBUG << "Computing guide tree - " << fixed << setprecision(1) 
				<< 100.0 * ((double)i * (i + 1) / 2) / ((double)n_seq * (n_seq + 1) / 2) << "\%    (" << i << " of " << n_seq << ")  \r";
		}

		slq.GetSolution(i, sim_vector);

		auto p_lambda = lambda.begin();
		auto p_sim_vector = (*sim_vector).begin();
		auto p_pi = pi.begin();

		for (int j = 0; j < i; ++j)
		{
			next = pi[j];

#ifdef _MSC_VER					// Visual C++
			_mm_prefetch((const char*)&(*sim_vector)[*(p_pi + prefetch_offset)], 2);
#endif
#ifdef __GNUC__
			//			__builtin_prefetch((&(*sim_vector)[pi[j + prefetch_offset]]), 1, 2);
			__builtin_prefetch((&(*sim_vector)[*(p_pi + prefetch_offset)]), 1, 2);
#endif

			auto &x = (*sim_vector)[next];

			if (isgreater(*p_lambda, *p_sim_vector))
			{
				x = max(x, *p_sim_vector);
			}
			else
			{
				x = max(x, *p_lambda);
				*p_pi = i;
				*p_lambda = *p_sim_vector;
			}

			++p_pi;
			++p_lambda;
			++p_sim_vector;
		}

		slq.ReleaseSolution(i);

		p_pi = pi.begin();
		p_lambda = lambda.begin();
		for (int j = 0; j < i; ++j)
		{
#ifdef _MSC_VER					// Visual C++
			_mm_prefetch((const char*)&lambda[*(p_pi + prefetch_offset_2nd)], 0);
#endif
#ifdef __GNUC__
			__builtin_prefetch((&lambda[*(p_pi + prefetch_offset_2nd)]), 1, 0);
#endif

			next = *p_pi;
			if (isgreaterequal(lambda[next], *p_lambda))
				*p_pi = i;

			++p_pi;
			++p_lambda;
		}
	}

	for (auto p : workers)
	{
		p->join();
		delete p;
	}
	workers.clear();

	LOG_DEBUG << "Computing guide tree - 100.0\%                                        \r";
	
	vector<int> elements(n_seq - 1);
	for (int i = 0; i < n_seq - 1; ++i)
		elements[i] = i;

#ifdef DEBUG_MODE
	identity /= n_seq * (n_seq - 1) / 2.0;
#endif

	stable_sort(elements.begin(), elements.end(), [&](int x, int y) {
		return lambda[x] > lambda[y];
	});

	vector<int> index(n_seq);
	for (int i = 0; i < n_seq; ++i)
		index[i] = i;

	for (int i = 0; i < n_seq - 1; ++i)
	{
		int j = elements[i];
		next = pi[j];
		tree.push_back(make_pair(index[j], index[next]));
		index[next] = n_seq + i;
	}
}

// *******************************************************************
void GuideTree::computeSingleLinkage_serial(std::vector<CSequence*>sequences, tree_structure& tree) 
{
	int next;
	int n_seq = sequences.size();

	int prefetch_offset = 64 * 2;
	int prefetch_offset_2nd = 128 * 2;

	vector<int> pi(n_seq + max(prefetch_offset, prefetch_offset_2nd), 0);
	vector<double> lambda(n_seq);
	vector<double> sim_vector(n_seq);

	CLCSBP lcsbp(instruction_set);


	// Single linkage algorithm is here
	for (int i = 0; i < n_seq; ++i)
	{
		pi[i] = i;
		lambda[i] = -infty_double;

		if (i % (100) == 0) {
			LOG_DEBUG << "Computing guide tree - " << fixed << setprecision(1)
				<< 100.0 * ((double)i * (i + 1) / 2) / ((double)n_seq * (n_seq + 1) / 2) << "\%    (" << i << " of " << n_seq << ")  \r";
		}

		calculateSimilarityVector<CSequence*, double>(
			sequences[i],
			sequences.data(),
			i,
			sim_vector.data(),
			lcsbp);

		auto p_lambda = lambda.begin();
		auto p_sim_vector = sim_vector.begin();
		auto p_pi = pi.begin();

		for (int j = 0; j < i; ++j)
		{
			next = pi[j];

#ifdef _MSC_VER					// Visual C++
			_mm_prefetch((const char*)&(sim_vector)[*(p_pi + prefetch_offset)], 2);
#endif
#ifdef __GNUC__
			//			__builtin_prefetch((&(*sim_vector)[pi[j + prefetch_offset]]), 1, 2);
			__builtin_prefetch((&(sim_vector)[*(p_pi + prefetch_offset)]), 1, 2);
#endif

			auto &x = (sim_vector)[next];

			if (isgreater(*p_lambda, *p_sim_vector))
			{
				x = max(x, *p_sim_vector);
			}
			else
			{
				x = max(x, *p_lambda);
				*p_pi = i;
				*p_lambda = *p_sim_vector;
			}

			++p_pi;
			++p_lambda;
			++p_sim_vector;
		}

		p_pi = pi.begin();
		p_lambda = lambda.begin();
		for (int j = 0; j < i; ++j)
		{
#ifdef _MSC_VER					// Visual C++
			_mm_prefetch((const char*)&lambda[*(p_pi + prefetch_offset_2nd)], 0);
#endif
#ifdef __GNUC__
			__builtin_prefetch((&lambda[*(p_pi + prefetch_offset_2nd)]), 1, 0);
#endif

			next = *p_pi;
			if (isgreaterequal(lambda[next], *p_lambda))
				*p_pi = i;

			++p_pi;
			++p_lambda;
		}
	}

	LOG_DEBUG << "Computing guide tree - 100.0\%                                        \r";

	vector<int> elements(n_seq - 1);
	for (int i = 0; i < n_seq - 1; ++i)
		elements[i] = i;

#ifdef DEBUG_MODE
	identity /= n_seq * (n_seq - 1) / 2.0;
#endif

	stable_sort(elements.begin(), elements.end(), [&](int x, int y) {
		return lambda[x] > lambda[y];
	});

	vector<int> index(n_seq);
	for (int i = 0; i < n_seq; ++i)
		index[i] = i;

	for (int i = 0; i < n_seq - 1; ++i)
	{
		int j = elements[i];
		next = pi[j];
		tree.push_back(make_pair(index[j], index[next]));
		index[next] = n_seq + i;
	}
}

// *******************************************************************
void GuideTree::calculateDistances(std::vector<CSequence>& sequences, UPGMA_dist_t *dist_matrix)
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
			calculateSimilarityVector<CSequence, UPGMA_dist_t, true>(
				&(*sequences)[row_id],
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
void GuideTree::computeUPGMA(UPGMA_dist_t* distances, size_t n_seq, tree_structure& tree)
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
void GuideTree::computePartTree(std::vector<CSequence>& sequences, GT_method::Value method, tree_structure& tree) {

	// Temporarily resize the sequences by adding unknown symbols (just to simplify the pairwise comparison)
	const auto& longest = sequences[0]; // sequences are sorted w.r.t. length
	for (int i = 1; i < sequences.size(); ++i)
		sequences[i].data.resize(longest.length, UNKNOWN_SYMBOL);

	// create vector of pointers to be passed to the recursion
	std::vector<CSequence*> sequencePtrs(sequences.size());
	std::transform(sequences.begin(), sequences.end(), sequencePtrs.begin(), [](CSequence& s)->CSequence* { return &s; });

	partTreeStep(sequencePtrs, method, tree);

	// Bring the sequences to the valid length
	for (int i = 1; i < sequences.size(); ++i)
		sequences[i].data.resize(sequences[i].length, UNKNOWN_SYMBOL);
}

// *******************************************************************
void GuideTree::partTreeStep(std::vector<CSequence*>& sequences, GT_method::Value method, tree_structure& tree)
{
	size_t n_seqs = sequences.size();
	CLCSBP lcsbp(instruction_set);
	
	if (n_seqs > parttree_size) {
		// calculate distances 0'th (longest) vs all (first row)
		float* similarities = new float[n_seqs * 2]; // second row will be used later
		
		float* similarity_row = similarities;
		calculateSimilarityVector<CSequence*, float, false>(sequences[0], sequences.data(), n_seqs, similarity_row, lcsbp);
		
		// two first seeds: the longest (0'th), the furthest 
		size_t* seed_ids = new size_t[parttree_size];
		seed_ids[0] = 0;
		seed_ids[1] = std::min_element(similarity_row + 1, similarity_row + n_seqs) - similarity_row;;
		
		// select randomly other seeds 
		std::uniform_int_distribution<size_t> dist(1, n_seqs - 1);
		std::mt19937 mt;
		std::generate(seed_ids + 2, seed_ids + parttree_size, [&dist, &mt]()->size_t { return dist(mt); });
		std::stable_sort(seed_ids, seed_ids + parttree_size);
		size_t n_seeds = std::unique(seed_ids, seed_ids + parttree_size) - seed_ids;

		//
		// --- Clustering ---
		//

		// assume all sequences are clustered to 0'th seed at the beginning
		std::vector<CSequence*> seeds(n_seeds);
		int* assignments = new int[n_seqs];
		std::fill_n(assignments, n_seqs, 0);

		// make assignments of all sequences to seeds
		float* current_row = similarities + n_seqs;
		seeds[0] = sequences[seed_ids[0]];
		for (int k = 1; k < seeds.size(); ++k) {
			seeds[k] = sequences[seed_ids[k]];
			calculateSimilarityVector<CSequence*, float, false>(seeds[k], sequences.data(), n_seqs, current_row, lcsbp);
			
			for (size_t j = 0; j < n_seqs; ++j) {
				if (current_row[j] > similarity_row[j]) { 	// use similarity_row for storing maximum similarities 
					similarity_row[j] = current_row[j];
					assignments[j] = k;
				}
			}
		}

		// get histogram
		int* histogram = new int[seeds.size()];
		std::fill_n(histogram, seeds.size(), 0);
		for (size_t j = 0; j < n_seqs; ++j) {
			++histogram[assignments[j]];
		}

		// reserve memory for subgroups
		std::vector<std::vector<CSequence*>> subgroups(seeds.size());
		for (int k = 0; k < seeds.size(); ++k) {
			subgroups[k].reserve(histogram[k]);
			assignments[seed_ids[k]] = -1; // mark seeds as not assigned to anything
		}

		// add sequences to subgroups
		for (size_t j = 0; j < n_seqs; ++j) {
			if (assignments[j] >= 0) { // do not assign seeds themselves
				subgroups[assignments[j]].push_back(sequences[j]);
			}
		}

		// delete everything before proceeding with the recursion
		delete[] histogram;
		delete[] assignments;
		delete[] seed_ids;
		delete[] similarities;

		// process child nodes
		std::vector<size_t> subroots(seeds.size());
		for (int k = 0; k < seeds.size(); ++k) {
			auto &subgroup = subgroups[k];

			if (subgroup.size() > 1) {
				// create a subtree
				partTreeStep(subgroup, method, tree);
				// make an intermediate node to join a subtree with a seed
				tree.push_back(node_t(seeds[k]->sequence_no, tree.size() - 1));
				subroots[k] = tree.size() - 1;
			}
			else if (subgroup.size() == 1) {
				// merge seed with the only element in the group
				tree.push_back(node_t(seeds[k]->sequence_no, subgroup.front()->sequence_no));
				subroots[k] = tree.size() - 1;
			}
		}

		size_t previousTop = tree.size();

		// generate partial tree for the seeds
		if (method == GT_method::parttree_SLINK) {
			computeSingleLinkage_serial(seeds, tree);
		}
		else if (method == GT_method::parttree_UPGMA){
			UPGMA_dist_t* matrix = TriangleMatrix::allocate<UPGMA_dist_t>(seeds.size());
			calculateSimilarityMatrix<CSequence*, UPGMA_dist_t, true>(seeds.data(), seeds.size(), matrix, lcsbp);
			computeUPGMA(matrix, seeds.size(), tree);
			delete[] matrix;
		}
		else {
			throw std::runtime_error("Tree type not supported by PartTree");
		}

		// correct node identifiers in the guide tree
		for (size_t node_id = previousTop; node_id < tree.size(); ++node_id) {
			auto& node  = tree[node_id];
		
			node.first = node.first < seeds.size()	
				? (subgroups[node.first].size() ? subroots[node.first] : seeds[node.first]->sequence_no)	// case: seed id - change to subroot or seq id 
				: node.first + previousTop - seeds.size();	 // case: intermediate node

			node.second = node.second < seeds.size()
				? (subgroups[node.second].size() ? subroots[node.second] : seeds[node.second]->sequence_no)	// case: seed id - change to subroot or seq id						 // case: seed id - change to subroot
				: node.second + previousTop - seeds.size();	 // case: intermediate node
		}
	}
	else {

		size_t previousTop = tree.size();

		// generate tree from all sequences
		if (method == GT_method::parttree_SLINK) {
			computeSingleLinkage_serial(sequences, tree);
		}
		else if (method == GT_method::parttree_UPGMA) {
			UPGMA_dist_t* matrix = TriangleMatrix::allocate<UPGMA_dist_t>(sequences.size());
			calculateSimilarityMatrix<CSequence*, UPGMA_dist_t, true>(sequences.data(), sequences.size(), matrix, lcsbp);
			computeUPGMA(matrix, sequences.size(), tree);
			delete[] matrix;
		}
		else {
			throw std::runtime_error("Tree type not supported by PartTree");
		}

		// correct node identifiers in the guide tree
		if (previousTop > sequences.size()) {
			for (size_t node_id = previousTop; node_id < tree.size(); ++node_id) {
				auto& node = tree[node_id];
				
				node.first = node.first < sequences.size()
					? sequences[node.first]->sequence_no				// case: sequence id 
					: node.first + previousTop - sequences.size();	// case: intermediate node

				node.second = node.second < sequences.size()
					? sequences[node.second]->sequence_no				// case: sequence id 
					: node.second + previousTop - sequences.size();	// case: intermediate node
			}
		}
	}
}

// *******************************************************************
#ifdef DEVELOPER_MODE
void GuideTree::computeChained(std::vector<CSequence>& sequences)
{
	mt19937 rnd;

	if (sequences.size() < 2)
		return;

	vector<int> idx(sequences.size());

	for (int i = 0; i < sequences.size(); ++i)
		idx[i] = i;

	random_device rd;

	// Skip some number of initial values
	for (int i = 0; i < seed; ++i)
		rd();

	mt19937 g(rd());

	shuffle(idx.begin(), idx.end(), g);

	guide_tree.push_back(make_pair(idx[0], idx[1]));

	for (int i = 2; i < sequences.size(); ++i)
		guide_tree.push_back(make_pair(idx[i], guide_tree.size() - 1));
}


// *******************************************************************
size_t GuideTree::subTreeSize(
	const std::vector<CSequence>& sequences, 
	const std::vector<CSequence>& ref_sequences,
	const std::set<int> &seq_ids)
{
	vector<pair<int, set<int>>> node_stats;
	set<int> tmp_set;
	int n_seq = sequences.size();

	// Calculate stats for single sequence nodes
	for (int i = 0; i < n_seq; ++i)
	{
		tmp_set.clear();
		if (seq_ids.count(i) > 0)
			tmp_set.insert(i);

		node_stats.push_back(make_pair(1, tmp_set));
	}

	// Calculate stats for internal nodes
	for (int i = n_seq; i < 2 * n_seq - 1; ++i)
	{
		tmp_set.clear();
		auto &x = node_stats[guide_tree[i].first];
		auto &y = node_stats[guide_tree[i].second];
		tmp_set.insert(x.second.begin(), x.second.end());
		tmp_set.insert(y.second.begin(), y.second.end());

		node_stats.push_back(make_pair(x.first + y.first, tmp_set));

		if (tmp_set.size() == ref_sequences.size())
			return node_stats.back().first;
	}

	return 0;
}


// *******************************************************************
size_t GuideTree::refSequencesSubTreeSize(
	const vector<CSequence>& sequences, 
	const vector<CSequence>& ref_sequences,
	double *monte_carlo_subtree_size)
{
	const int monte_carlo_trials = 1000;
	
	set<int> ref_seq_ids;
	int n_seq = sequences.size();
	int r = 0;

	if (ref_sequences.size() == 1)
		return 1;

	// Find the ids of the referential sequences in the input file
	for (int i = 0; i < n_seq; ++i)
	{
		bool is_ref = false;
		for (auto &y : ref_sequences)
			if (sequences[i].id == y.id)
				is_ref = true;

		if (is_ref)
			ref_seq_ids.insert(i);
	}

	r = subTreeSize(sequences, ref_sequences, ref_seq_ids);

	if (monte_carlo_subtree_size)
	{
		mt19937 mt;
		double mc_r = 0;

		for (int i = 0; i < monte_carlo_trials; ++i)
		{
			set<int> mc_seq_ids;

			while (mc_seq_ids.size() < ref_seq_ids.size())
				mc_seq_ids.insert(mt() % n_seq);
			mc_r += subTreeSize(sequences, ref_sequences, mc_seq_ids);
		}

		*monte_carlo_subtree_size = mc_r / (double)monte_carlo_trials;
	}

	return r;
}

#endif