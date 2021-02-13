/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#pragma once

#include "TreeDefs.h"
#include "sequence.h"

#include <vector>


class GuideTree {
public:

	tree_structure & raw() { return guide_tree; }

	GuideTree() {}

	size_t getSequenceCount() const { return (guide_tree.size() + 1) / 2; }

	bool loadNewick(
		const std::string& file_name,
		std::vector<CSequence>& sequences);

	bool saveNewick(
		const std::string& file_name,
		const std::vector<CSequence>& sequences) const;

	uint64_t calculateSackinIndex();

protected:

	tree_structure guide_tree;

};


