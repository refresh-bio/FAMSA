/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#include <stdexcept>
#include <algorithm>
#include <ostream>
#include <sstream>
#include <string.h>

#include <boost/spirit/home/classic.hpp>

#include "NewickParser.h"
#include "TreeGrammar.h"
#include "../utils/log.h"

#undef min

namespace bs = boost::spirit::classic;
using namespace quickprobs;



void NewickParser::parse(
	const std::vector<CSequence>& sequences,
	const std::string& description,
	std::vector<pair<int, int>>& guideTree)
{
	if (description.length() == 0) {
		throw std::runtime_error("Error while parsing Newick tree: empty description.");
	}

	LOG_VERBOSE << "Parsing guide tree..." << endl
		<< "Description length: " << description.length() << endl;
	
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

		LOG_VERBOSE << "Unparsed characters: " << strlen(info.stop) << endl
			<< "Unparsed context: " << unparsed << endl
			<< "Parsed anything: " << info.hit << endl
			<< "Parsed everything: " << info.full << endl;
	}

	if (info.full == false) {
		throw std::runtime_error("Error while parsing Newick tree: invalid format.");
	}
}

void NewickParser::store(
	const std::vector<CSequence>& sequences,
	const std::vector<pair<int, int>>& guideTree,
	std::string& description) {

	ostringstream oss;
	oss << "(";
	storeBranch(sequences, guideTree, guideTree.back().first, oss);
	oss << ",";
	storeBranch(sequences, guideTree, guideTree.back().second, oss);
	oss << ");";
	description = oss.str();
}

std::ostream& NewickParser::storeBranch(const std::vector<CSequence>& sequences, const std::vector<pair<int, int>>& guideTree, int index, std::ostream& oss) {

	if (index < (int)sequences.size()) {

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
