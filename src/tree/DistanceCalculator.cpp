#include "DistanceCalculator.h"
#include "AbstractTreeGenerator.hpp"
#include "SingleLinkage.h"

#include "../utils/conversion.h"
#include "SingleLinkageQueue.h"

#include <thread>

template <Distance _distance>
void DistanceCalculator<_distance>::run(std::vector<CSequence*>& sequences, tree_structure& tree) {
	
	std::ofstream ofs(out_file);
	// put header line only when full matrix is needed
	if (generate_square_matrix) {
		for (const auto seq : sequences) {
			ofs << ',' << seq->id.c_str() + 1;
		}
		ofs << endl;
	}

	//
	int n_seqs = (int)sequences.size();
	CSingleLinkageQueue<float> queue(&sequences, (uint32_t) sequences.size(), n_threads * 8);
	std::vector<std::thread> workers(n_threads);

	// run workers
	for (int tid = 0; tid < n_threads; ++tid) {
		workers[tid] = thread([&queue, this]() {

			CLCSBP lcsbp(instruction_set);
			int row_id;
			std::vector<CSequence*>* sequences;
			vector<float>* dist_vector;
			vector<float> loc_dist_vector;

			if (calculate_pid) {
				Transform<float, Distance::pairwise_identity> transform;

				while (queue.GetTask(row_id, sequences, dist_vector)) {
					loc_dist_vector.resize(dist_vector->size());
					int to_calculate = generate_square_matrix ? (int)sequences->size() : row_id;

					calculateDistanceVector<CSequence*, float, decltype(transform)>(
						transform,
						(*sequences)[row_id],
						sequences->data(),
						to_calculate,
						loc_dist_vector.data(),
						lcsbp);

					//loc_dist_vector[row_id] = 1.0;
					swap(*dist_vector, loc_dist_vector);

					//cout << "push " << row_id << endl;
					queue.RegisterSolution(row_id);
				}
			}
			else {
				Transform<double, _distance> transform;

				while (queue.GetTask(row_id, sequences, dist_vector)) {
					loc_dist_vector.resize(dist_vector->size());
					int to_calculate = generate_square_matrix ? (int)sequences->size() : row_id;

					calculateDistanceVector<CSequence*, float, decltype(transform)>(
						transform,
						(*sequences)[row_id],
						sequences->data(),
						to_calculate,
						loc_dist_vector.data(),
						lcsbp);

					swap(*dist_vector, loc_dist_vector);

					//cout << "push " << row_id << endl;
					queue.RegisterSolution(row_id);
				}
			}
		
		});
	}

	char* out_row = new char[10000 + sequences.size() * 100];
	char *ptr = out_row;

	// Gather results in one thread
	for (int row_id = 0; row_id < n_seqs; ++row_id) {

		vector<float>* dist_vector;

		queue.GetSolution(row_id, dist_vector);
		//cout << "pop " << row_id << endl;
		//if ((row_id + 1) % 100 == 0) {
		//	cout << "\r" << row_id + 1 << "...                      " << std::flush;
		//}

		ptr = out_row;
		ptr += sprintf(ptr, "%s,", sequences[row_id]->id.c_str() + 1);

		if (generate_square_matrix) {
			ptr += num2str(dist_vector->data(), dist_vector->size(), ',', ptr);
		}
		else {
			ptr += num2str(dist_vector->data(), row_id, ',', ptr);
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



// *******************************************************************
// Explicit template specializations for specified distance measures

template class DistanceCalculator<Distance::indel_div_lcs>;
template class DistanceCalculator<Distance::sqrt_indel_div_lcs>;
