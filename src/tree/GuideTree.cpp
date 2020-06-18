/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

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
uint64_t GuideTree::calculateSackinIndex() {
	
	uint64_t idx = 0;
	size_t n_sequences = getSequenceCount();
	
	if (n_sequences) {
		std::vector<size_t> depths(this->guide_tree.size());
		for (size_t i = guide_tree.size() - 1; i >= n_sequences; --i)
		{
			depths[guide_tree[i].first] = depths[i] + 1;
			depths[guide_tree[i].second] = depths[i] + 1;
		}


		for (int i = 0; i < n_sequences; ++i) {
			idx += depths[i] + 1;
		}
	}

	return idx;
}




