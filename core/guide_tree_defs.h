/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is http://sun.aei.polsl.pl/REFRESH/famsa

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#pragma once

#include <string>
#include <stdexcept>

// UPGMA defines and consts
typedef float UPGMA_dist_t;
const UPGMA_dist_t BIG_DIST = (UPGMA_dist_t) 1e29;

// Class representing guide tree method
class GT_method {
public:
	enum Value { SLINK, UPGMA, parttree_SLINK, parttree_UPGMA, chained, imported };

	static std::string toString(Value v) {
		switch (v) {
		case SLINK:				return "sl (single linkage)";
		case UPGMA:				return "upgma";
		case parttree_SLINK:	return "parttree_sl (parttree + single linkage)";
		case parttree_UPGMA:	return "patttree_upgma (parttree + upgma)";
		case chained:			return "chained";
		case imported:			return "import";
		}
	}

	static Value fromString(const std::string& name) {
		if (name == "sl") { return SLINK; }
		if (name == "upgma") { return UPGMA; }
		if (name == "parttree_sl") { return parttree_SLINK; }
		if (name == "parttree_upgma") { return parttree_UPGMA; }
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

	static size_t access(size_t i, size_t j) {
		if (i >= j)
			return j + (i * (i - 1)) / 2;
		else
			return i + (j * (j - 1)) / 2;
	}
};