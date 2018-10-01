/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#ifndef _NEWICK_TREE_H
#define _NEWICK_TREE_H
#include <string>
#include <vector>
#include <ostream>

class CSequence;

class NewickTree {

protected:
	bool verbose;

public:
	NewickTree(bool verbose) : verbose(verbose) {}

	void parse(
		const std::vector<CSequence>& sequences,
		const std::string& description,
		std::vector<std::pair<int, int>>& guideTree);

	void store(
		const std::vector<CSequence>& sequences,
		const std::vector<std::pair<int, int>>& guideTree,
		std::string& description);

private:
	std::ostream&  storeBranch(const std::vector<CSequence>& sequences, const std::vector<std::pair<int, int>>& guideTree, int index, std::ostream& oss);
};


#endif