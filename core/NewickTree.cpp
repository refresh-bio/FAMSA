/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

Version: 1.1
Date   : 2016-06-29
*/

#include <stdexcept>
#include <algorithm>
#include <string.h>

#include <boost/spirit/home/classic.hpp>

#include "NewickTree.h"
#include "TreeGrammar.h"

#undef min

namespace bs = boost::spirit::classic;
using namespace quickprobs;



void parseNewickTree(
	const std::vector<CSequence>& sequences,
	const std::string& description,
	std::vector<pair<int, int>>& guideTree,
	bool verbose)
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

