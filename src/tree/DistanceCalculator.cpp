#include "DistanceCalculator.h"
#include "AbstractTreeGenerator.hpp"
#include "SingleLinkage.h"

#include "../utils/conversion.h"
#include "SingleLinkageQueue.h"

#include <thread>

void DistanceCalculator::run(std::vector<CSequence>& sequences, tree_structure& tree) {
	
	std::ofstream ofs(out_file);
	// put header line only when full matrix is needed
	if (generate_square_matrix) {
		for (const CSequence& seq : sequences) {
			ofs << ',' << seq.id.c_str() + 1;
		}
		ofs << endl;
	}

	//
	CSingleLinkageQueue<float> queue(&sequences, sequences.size(), n_threads * 8);
	std::vector<std::thread> workers(n_threads);

	// run workers
	for (int tid = 0; tid < n_threads; ++tid) {
		workers[tid] = thread([&queue, tid, this]() {

			CLCSBP lcsbp(instruction_set);
			int row_id;
			std::vector<CSequence>* sequences;
			vector<float>* sim_vector;
			vector<float> loc_sim_vector;

			if (calculate_pid) {
				Transform<double, Measure::SimilarityDefault> transform;

				while (queue.GetTask(row_id, sequences, sim_vector)) {
					loc_sim_vector.resize(sim_vector->size());
					int to_calculate = generate_square_matrix ? sequences->size() : row_id;

					calculateSimilarityVector<CSequence, float, decltype(transform)>(
						transform,
						(*sequences)[row_id],
						sequences->data(),
						to_calculate,
						loc_sim_vector.data(),
						lcsbp);

					//loc_sim_vector[row_id] = 1.0;
					swap(*sim_vector, loc_sim_vector);

					//cout << "push " << row_id << endl;
					queue.RegisterSolution(row_id);
				}
			}
			else {
				Transform<double, Measure::DistanceReciprocal> transform;

				while (queue.GetTask(row_id, sequences, sim_vector)) {
					loc_sim_vector.resize(sim_vector->size());
					int to_calculate = generate_square_matrix ? sequences->size() : row_id;

					calculateSimilarityVector<CSequence, float, decltype(transform)>(
						transform,
						(*sequences)[row_id],
						sequences->data(),
						to_calculate,
						loc_sim_vector.data(),
						lcsbp);

					swap(*sim_vector, loc_sim_vector);

					//cout << "push " << row_id << endl;
					queue.RegisterSolution(row_id);
				}
			}
		
		});
	}

	char* out_row = new char[10000 + sequences.size() * 100];
	char *ptr = out_row;

	cout << endl;

	// Gather results in one thread
	for (int row_id = 0; row_id < sequences.size(); ++row_id) {

		vector<float>* sim_vector;

		queue.GetSolution(row_id, sim_vector);
		//cout << "pop " << row_id << endl;
		//if ((row_id + 1) % 100 == 0) {
		//	cout << "\r" << row_id + 1 << "...                      " << std::flush;
		//}

		ptr = out_row;
		ptr += sprintf(ptr, "%s,", sequences[row_id].id.c_str() + 1);

		if (generate_square_matrix) {
			ptr += num2str(sim_vector->data(), sim_vector->size(), ',', ptr);
		}
		else {
			ptr += num2str(sim_vector->data(), row_id, ',', ptr);
		}

		queue.ReleaseSolution(row_id);
		//cout << "return " << row_id << endl;
		--ptr; 
		*ptr++ = '\n';
		ofs.write(out_row, ptr - out_row);
	}

	delete[] out_row;

	// make sure all threads have finished
	for (auto &w : workers) {
		w.join();
	}

}