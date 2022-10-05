/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#include <stdexcept>
#include <algorithm>
#include <ostream>
#include <sstream>
#include <map>
#include <string.h>


#include "NewickParser.h"
#include "../utils/log.h"

#undef min

using namespace std;

void NewickParser::parse(
	const std::vector<CSequence>& sequences,
	const std::string& description,
	std::vector<std::pair<int, int>>& guideTree)
{
	if (description.length() == 0) {
		throw std::runtime_error("Error while parsing Newick tree: empty description.");
	}

	LOG_VERBOSE << endl << "Newick description length: " << description.length() << endl;

	// map sequence names to ids
	std::map<std::string, int> sequencesToIds;
	guideTree.resize(2 * sequences.size(), std::pair<int, int>(-1, -1)); // add extra node at the end - will be removed

	// fill in mappings
	int n_seqs = (int)sequences.size();
	for (int i = 0; i < n_seqs; ++i) {
		const auto& seq = sequences[i];
		if (seq.id[0] == '>') {
			sequencesToIds[seq.id.substr(1, seq.id.size())] = i; // omit >
		}
		else {
			sequencesToIds[seq.id] = i;
		}
	}

	const char* p = description.c_str();
	const char* end = p + description.size();

	int cur_pos = (int)guideTree.size() - 1; 
	int free_pos = cur_pos - 1;
	bool secondBranch = false;
	std::vector<int> prevs(guideTree.size() + 1, -1);

	while (p < end) {
		// subtree begin
		if (*p == '(') {
			auto& out_branch = secondBranch ? guideTree[cur_pos].second : guideTree[cur_pos].first;
			out_branch = free_pos;
			prevs[free_pos] = cur_pos;
			cur_pos = free_pos;
			
			++p;
			--free_pos;
			secondBranch = false;
		}
		else if (*p == ',') {
			++p;
			secondBranch = true;
		}
		else if (*p == ')') {
			++p;
			cur_pos = prevs[cur_pos];
		}
		else if (*p == ':') {
			// branch length
			++p;
			char* len_end;
			strtof(p, &len_end);
			p = len_end;
		}
		else if (isspace(*p)) {
			++p; // ignore whitespaces
		}
		else {
			// find end of the sequence name
			const char* name_end = std::find_if(p, end, [](char c) { return c == ')' || c == ',' || c == ':' || c == '(';  });
			string name(p, name_end);
			int id = sequencesToIds[name];

			auto& out_branch = secondBranch ? guideTree[cur_pos].second : guideTree[cur_pos].first;
			out_branch = id;
			p = name_end;
		}
	}

	guideTree.resize(guideTree.size() - 1); // remove an extra node
}

void NewickParser::store(
	const std::vector<CSequence>& sequences,
	const std::vector<std::pair<int, int>>& guideTree,
	std::string& description) {

	ostringstream oss;
	
	std::vector<int> prevs(guideTree.size() + 1, -1);
	std::vector<int> num_visits(guideTree.size() + 1, 0);
	int last_pos = guideTree.size() - 1;
	int cur_pos = last_pos;
	
	while (true) {
		
		if (cur_pos < (int)sequences.size()) {
			// if sequence was reached
			
			const char* begin = sequences[cur_pos].id.c_str();
			
			// remove trailing '<' if present
			if (*begin == '>') { ++begin;  }

			oss << begin << ":1.0";
			cur_pos = prevs[cur_pos]; 

		}
		else {
			// if internal node

			if (num_visits[cur_pos] == 0) {
				// no visits - left branch
				oss << '(';
				int dest_pos = guideTree[cur_pos].first;
				++num_visits[cur_pos];
				prevs[dest_pos] = cur_pos;
				cur_pos = dest_pos;
				
			}
			else if (num_visits[cur_pos] == 1) {
				// one visit - right branch
				oss << ',';
				int dest_pos = guideTree[cur_pos].second;
				++num_visits[cur_pos];
				prevs[dest_pos] = cur_pos;
				cur_pos = dest_pos;
			}
			else {
				// two visits - node processed
				if (cur_pos == last_pos) {
					// root processed
					oss << ");";
					break;
				}

				oss << "):1.0";
				++num_visits[cur_pos];
				cur_pos = prevs[cur_pos];
			}
		}
	}

	description = oss.str();
}
