/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#pragma once

#include "TreeDefs.h"
#include "../core/sequence.h"

#include <vector>


class GuideTree {
public:
	
	tree_structure & raw() { return guide_tree; }

	GuideTree() {}

	int getSequenceCount() const { return (int)(guide_tree.size() + 1) / 2; }

	bool loadNewick(
		const std::string& file_name,
		std::vector<CSequence>& sequences);

	bool saveNewick(
		const std::string& file_name,
		const std::vector<CSequence>& sequences) const;

	int64_t calculateSackinIndex();

	void toUnique(const std::vector<int>& original2unique, int n_uniques);
	void fromUnique(const std::vector<int>& original2unique);

protected:

	tree_structure guide_tree;

};


