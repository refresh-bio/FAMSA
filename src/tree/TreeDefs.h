/*
This file is a part of FAMSA software distributed under GNU GPL 3 licence.
The homepage of the FAMSA project is https://github.com/refresh-bio/FAMSA

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys

*/
#pragma once

#include <string>
#include <stdexcept>
#include <vector>

// tree structure
using node_t = std::pair<int, int>;
using tree_structure = std::vector<node_t>;


// class representing distance
enum class Distance {
	indel_div_lcs,
	sqrt_indel_div_lcs,
	neg_lcs_div_indel,
	neg_lcs_div_minlen,
	neg_lcs_div_len_corrected,
	pairwise_identity
};

inline static Distance str2dist(const std::string& s) {
	if (s == "indel_div_lcs") { return Distance::indel_div_lcs; }
	else if (s == "sqrt_indel_div_lcs") { return Distance::sqrt_indel_div_lcs; }
	else if (s == "neg_lcs_div_indel") { return Distance::neg_lcs_div_indel; }
	else if (s == "neg_lcs_div_minlen") { return Distance::neg_lcs_div_minlen; }
	else if (s == "neg_lcs_div_len_corrected") { return Distance::neg_lcs_div_len_corrected; }
	else {
		throw new std::runtime_error("Error: Illegal pairwise distance measure.");
	}

	return Distance::indel_div_lcs;
}

inline static std::string dist2str(Distance d) {
	switch (d) {
	case Distance::indel_div_lcs:				return "indel_div_lcs";
	case Distance::sqrt_indel_div_lcs:			return "sqrt_indel_div_lcs";
	case Distance::neg_lcs_div_indel:			return "neg_lcs_div_indel";
	case Distance::neg_lcs_div_minlen:			return "neg_lcs_div_minlen";
	case Distance::neg_lcs_div_len_corrected:	return "neg_lcs_div_len_corrected";
	default:
		throw new std::runtime_error("Error: Illegal pairwise distance measure.");
	}

	return "Unknown";
}

// Class representing guide tree method
class GT {
public:
	enum Method { SLINK, MST_Prim, UPGMA, UPGMA_modified, NJ, chained, imported };
	enum Heuristic { None, PartTree, ClusterTree };
	
	static std::string toString(Method v) {
		switch (v) {
		case SLINK:				return "single linkage (SLINK)";
		case MST_Prim:			return "single linkage (MST+Prim)";
		case UPGMA:				return "upgma";
		case UPGMA_modified:	return "upgma_modified";
		case NJ:				return "nj";
		case chained:			return "chained";
		case imported:			return "import";
		default:
			throw new std::runtime_error("Error: Illegal guide tree method.");
		}

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
		if (name == "sl") { return MST_Prim; }
		if (name == "slink") { return SLINK; }
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