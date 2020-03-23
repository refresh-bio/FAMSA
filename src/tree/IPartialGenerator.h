#pragma once

#include "TreeDefs.h"

#include <vector>

class CSequence;


class IPartialGenerator {

public:
	virtual void runPartial(std::vector<CSequence*>& sequences, tree_structure& tree) = 0;

};
