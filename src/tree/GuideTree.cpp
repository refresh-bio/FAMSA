/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

Note: The file contains code for the UPGMA method (of high memory consumption)
These functions was borrowed and adpoted from MUSCLE 3.8.1551 by Robert Edgar
The time and memory consumption is O(k^2)

*/
#include "GuideTree.h"
#include "NewickParser.h"

#include <sstream>
#include <algorithm>

using namespace std;

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
	NewickParser nw(false);
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
	NewickParser nw(false);
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

//	calculateDistances(sequences, g_Dist);

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


