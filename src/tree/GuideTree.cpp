/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#include "GuideTree.h"
#include "NewickParser.h"

#include <sstream>
#include <algorithm>
#include <numeric>

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
		throw std::runtime_error("Unable to open Newick file: " + file_name);
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
int64_t GuideTree::calculateSackinIndex() {
	
	uint64_t idx = 0;
	int n_sequences = getSequenceCount();
	
	if (n_sequences) {
		std::vector<int64_t> depths(this->guide_tree.size());
		for (int i = (int)guide_tree.size() - 1; i >= n_sequences; --i)
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

// *******************************************************************
void GuideTree::toUnique(const std::vector<int>& original2unique, int n_uniques) {

	int n_total_seqs = original2unique.size();
	auto& vt = guide_tree;
	
	// remove duplicated sequences from the imported tree
	int offset = n_total_seqs - n_uniques;
	
	vt.erase(vt.begin() + n_uniques, vt.begin() + n_total_seqs);

	std::vector<int> out_ids(vt.size());
	std::iota(out_ids.begin(), out_ids.begin() + n_uniques, 0);
	int n_dups = 0;

	auto is_duplicate = [&out_ids, n_uniques](int node_id) ->bool { return out_ids[node_id] < n_uniques; };

	for (int i = n_uniques; i < (int)vt.size(); ++i) {
		auto& node = vt[i];
		// correct indices
		node.first = (node.first < n_total_seqs) ? original2unique[node.first] : node.first - offset;
		node.second = (node.second < n_total_seqs) ? original2unique[node.second] : node.second - offset;

		if (node.first == node.second) {
			// merge a leaf with itself - make a duplicate node
			++n_dups;
			out_ids[i] = node.second;
		}
		else if (is_duplicate(node.first) && node.second == out_ids[node.first]) {
			// merge a duplicate node (first) with itself (second)
			++n_dups;
			out_ids[i] = node.second;
		}
		else if (is_duplicate(node.second) && node.first == out_ids[node.second]) {
			// merge a duplicate node (second) with itself (left)
			++n_dups;
			out_ids[i] = node.first;
		}
		else {
			// merge two different nodes
			node.first = out_ids[node.first];
			node.second = out_ids[node.second];
			out_ids[i] = i - n_dups;
		}
	}

	for (int i = n_uniques; i < (int)vt.size(); ++i) {
		if (!is_duplicate(i)) {
			vt[out_ids[i]] = vt[i];
		}
	}

	vt.erase(vt.end() - n_dups, vt.end());
}

// *******************************************************************
void GuideTree::fromUnique(const std::vector<int>& original2unique) {
	
	int n_total_seqs = (int)original2unique.size();
	int n_uniques = this->getSequenceCount();
	int n_dups = n_total_seqs - n_uniques;
	auto& vt = guide_tree;

	std::vector<std::vector<int>> unique2original(n_uniques, std::vector<int>());
	std::vector<int> out_ids(n_uniques, -1);
	std::iota(out_ids.begin(), out_ids.end(), 0);

	for (int i = 0; i < n_total_seqs; ++i) {
		unique2original[original2unique[i]].push_back(i);
	}

	// add duplicated leafs (n_dups) and nodes joining duplicated leafs (n_dups)
	vt.insert(vt.begin() + n_uniques, 2 * n_dups, node_t(-1, -1));

	// add subtrees joining duplicated leafs [n_unique + dups, n_unique + 2 * dups]
	int node_id = n_uniques + n_dups;
	for (int iu = 0; iu < n_uniques; ++iu) {
		const vector<int>& occs = unique2original[iu];

		// iterate over occurrences of unique
		for (int i = 1; i < (int)occs.size(); ++i, ++node_id) {
			if (i == 1) {
				// first duplication - join original leaf with duplicated leaf
				vt[node_id].first = occs[0];
				vt[node_id].second = occs[1];
			}
			else {
				// following duplicates - join original leaf with previous node
				vt[node_id].first = occs[i];
				vt[node_id].second = node_id - 1;
			}
		}

		if (occs.size() > 1) {
			out_ids[iu] = node_id - 1; // represent sequence by a subtree
		}
		else {
			out_ids[iu] = occs[0];		// represent by the original id
		}
	}

	// in the previously-existing nodes replace duplicated sequences with subtrees 
	for (int i = node_id; i < (int)vt.size(); ++i) {
		auto& node = vt[i];
		if (node.first < n_uniques) {
			node.first = out_ids[node.first];
		}
		else {
			node.first += 2 * n_dups; 
		}

		if (node.second < n_uniques) {
			node.second = out_ids[node.second];
		}
		else {
			node.second += 2 * n_dups;
		}
	}
}