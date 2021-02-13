/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#pragma once

#include <string>
#include <stdexcept>
#include <vector>

// tree structure
using node_t = std::pair<int, int>;
using tree_structure = std::vector<node_t>;



// Class representing guide tree method
class GT {
public:
	enum Method { SLINK, UPGMA, UPGMA_modified, NJ, chained, imported };
	enum Heuristic { None, PartTree, ClusterTree };

	static std::string toString(Method v) {
		switch (v) {
		case SLINK:				return "sl (single linkage)";
		case UPGMA:				return "upgma";
		case UPGMA_modified:	return "upgma_modified";
		case NJ:				return "nj";
		case chained:			return "chained";
		case imported:			return "import";
		}

		throw new std::runtime_error("Error: Illegal guide tree method.");
		return "Unknown";
	}

	static std::string toString(Heuristic v) {
		switch (v) {
		case None: return "None";
		case PartTree:	return "PartTree";
		case ClusterTree:	return "MedoidTree";
		}

		// something went wrong
		throw new std::runtime_error("Error: Illegal guide tree heuristic.");
		return "Unknown";
	}

	static Method fromString(const std::string& name) {
		if (name == "sl") { return SLINK; }
		if (name == "upgma") { return UPGMA; }
		if (name == "upgma_modified") { return UPGMA_modified; }
		if (name == "nj") { return NJ; }
		if (name == "import") { return imported; }
#ifdef DEVELOPER_MODE
		if (name == "chained") { return chained; }
#endif
		// something went wrong
		throw new std::runtime_error("Error: Illegal guide tree method.");
		
		return SLINK;
	}
};


class TriangleMatrix {
public:
	template <class T>
	static T* allocate(size_t size) {
		return new T[(size * (size - 1)) / 2];
	}

	static size_t access(int64_t i, int64_t j) {
		if (i >= j)
			return j + (i * (i - 1)) / 2;
		else
			return i + (j * (j - 1)) / 2;
	}
};