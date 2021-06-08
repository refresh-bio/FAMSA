#pragma once

#include "AbstractTreeGenerator.h"



#ifdef DEVELOPER_MODE
void GuideTree::computeChained(std::vector<CSequence>& sequences);
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

	guide_tree.emplace_back(idx[0], idx[1]);

	for (int i = 2; i < sequences.size(); ++i)
		guide_tree.emplace_back(idx[i], guide_tree.size() - 1);
}

#endif

