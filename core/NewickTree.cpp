/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/

#include <stdexcept>
#include <algorithm>
#include <ostream>
#include <sstream>
#include <string.h>

#include <boost/spirit/home/classic.hpp>

#include "NewickTree.h"
#include "TreeGrammar.h"

#undef min

namespace bs = boost::spirit::classic;
using namespace quickprobs;



void NewickTree::parse(
	const std::vector<CSequence>& sequences,
	const std::string& description,
	std::vector<pair<int, int>>& guideTree)
{
	if (description.length() == 0) {
		throw std::runtime_error("Error while parsing Newick tree: empty description.");
	}

	if (verbose) {
		cout << "Parsing guide tree..." << endl
			<< "Description length: " << description.length() << endl;
	}

	TreeGrammar grammar(sequences);
	
	string newDesc;

	if (description.back() == ';') {
		newDesc = description.substr(0, description.size() - 1);
	}
	else {
		newDesc = description;
	}

	auto info = bs::parse(newDesc.c_str(), grammar);
	guideTree = std::move(grammar.wrapper.guideTree);
	guideTree.resize(grammar.wrapper.currentInternalId);

	if (verbose) {
		string unparsed(info.stop, std::min((::size_t)50, strlen(info.stop)));

		cout << "Unparsed characters: " << strlen(info.stop) << endl
			<< "Unparsed context: " << unparsed << endl
			<< "Parsed anything: " << info.hit << endl
			<< "Parsed everything: " << info.full << endl;
	}

	if (info.full == false) {
		throw std::runtime_error("Error while parsing Newick tree: invalid format.");
	}
}

void NewickTree::store(
	const std::vector<CSequence>& sequences,
	const std::vector<pair<int, int>>& guideTree,
	std::string& description) {

	ostringstream oss;
	oss << "(";
	storeBranch(sequences, guideTree, guideTree.back().first, oss);
	oss << ",";
	storeBranch(sequences, guideTree, guideTree.back().second, oss);
	oss << ");";
	description = std::move(oss.str());
}

std::ostream& NewickTree::storeBranch(const std::vector<CSequence>& sequences, const std::vector<pair<int, int>>& guideTree, int index, std::ostream& oss) {

	if (index < sequences.size()) {

		if (sequences[index].id[0] == '>') {
			oss << sequences[index].id.substr(1) << ":1.0";
		}
		else {
			oss << sequences[index].id << ":1.0";
		}
	}
	else {
		oss << "(";
		storeBranch(sequences, guideTree, guideTree[index].first, oss) << ",";
		storeBranch(sequences, guideTree, guideTree[index].second, oss) << "):1.0";
	}
	
	return oss;
}
