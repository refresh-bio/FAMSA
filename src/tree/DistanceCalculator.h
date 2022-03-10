#pragma once
#include "AbstractTreeGenerator.h"

class DistanceCalculator : public AbstractTreeGenerator {
private:
	string out_file;
	bool generate_square_matrix;
	bool calculate_pid;
public:
	DistanceCalculator(double indel_exp, size_t n_threads, const string& out_file, bool generate_square_matrix, bool calculate_pid)
		: AbstractTreeGenerator(indel_exp, n_threads), 
		out_file(out_file), generate_square_matrix(generate_square_matrix), calculate_pid(calculate_pid) {}

protected:
	void run(std::vector<CSequence>& sequences, tree_structure& tree) override;

};