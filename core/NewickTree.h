/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

Version: 1.1
Date   : 2016-06-29
*/

#ifndef _NEWICK_TREE_H
#define _NEWICK_TREE_H
#include <string>
#include <vector>

class CSequence;

void parseNewickTree(
	const std::vector<CSequence>& sequences, 
	const std::string& description,
	std::vector<std::pair<int,int>>& guideTree,
	bool verbose);

#endif