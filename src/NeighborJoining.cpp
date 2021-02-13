#include "NeighborJoining.h"
#include "AbstractTreeGenerator.hpp"

#include "lcsbp.h"

#include <limits>

// *******************************************************************
void NeighborJoining::run(std::vector<CSequence>& sequences, tree_structure& tree) {
	float* distances = TriangleMatrix::allocate<float>(sequences.size());
	CLCSBP lcsbp(instruction_set);
	calculateSimilarityMatrix<CSequence, float, Measure::DistanceReciprocal>(sequences.data(), sequences.size(), distances, lcsbp);

	computeTree(distances, sequences.size(), tree);

	delete[] distances;
}

// *******************************************************************
void NeighborJoining::runPartial(std::vector<CSequence*>& sequences, tree_structure& tree) {
	float* distances = TriangleMatrix::allocate<float>(sequences.size());
	CLCSBP lcsbp(instruction_set);
	calculateSimilarityMatrix<CSequence*, float, Measure::DistanceReciprocal>(sequences.data(), sequences.size(), distances, lcsbp);

	computeTree(distances, sequences.size(), tree);

	delete[] distances;
}


// *******************************************************************
void NeighborJoining::computeTree(float* distances, size_t n_seq, tree_structure& tree) {

	struct cluster {
		float sum_of_dists;
		int row_id;
		int node_id;
	};

	float *D = distances;
	std::vector<cluster> clusters(n_seq);

	// initialize clusters
	for (size_t i = 0; i < n_seq; ++i) {
		auto& ci = clusters[i];
		ci.row_id = ci.node_id = i;
		ci.sum_of_dists = 0;

		for (size_t j = 0; j < n_seq; ++j) {
			if (i != j) {
				ci.sum_of_dists += D[TriangleMatrix::access(i, j)];
			}
		}
	}

	// merge clusters as long as there are two left
	for (int iter = 0, n_clusters = n_seq; n_clusters > 2; ++iter, --n_clusters) {

		// find minimum element in Q matrix
		float min_q = std::numeric_limits<float>::max();
		int min_i, min_j;


		for (int i = 0; i < n_clusters; ++i) {
			for (int j = i + 1; j < n_clusters; ++j) {
				const auto & ci = clusters[i];
				const auto & cj = clusters[j];

				float q = (n_clusters - 2) * D[TriangleMatrix::access(ci.row_id, cj.row_id)] - ci.sum_of_dists - cj.sum_of_dists;
				if (q < min_q) {
					min_q = q;
					min_i = i;
					min_j = j;
				}
			}
		}

		// merge two resulting clusters
		auto & ci = clusters[min_i];
		auto & cj = clusters[min_j];
		float Dij = D[TriangleMatrix::access(ci.row_id, cj.row_id)];


		// ci is going to be replaced by a new cluster
		tree.push_back(node_t(ci.node_id, cj.node_id));

		ci.sum_of_dists = 0;
		ci.node_id = n_seq + iter;

		// recalculate distances
		for (int k = 0; k < clusters.size(); ++k) {

			if (k != min_i && k != min_j) {
				auto & ck = clusters[k];

				float Dik = D[TriangleMatrix::access(ci.row_id, ck.row_id)];
				float Djk = D[TriangleMatrix::access(cj.row_id, ck.row_id)];

				// remove contribution of ci and cj from ck sums
				ck.sum_of_dists -= Dik + Djk;

				// ci is replaced
				Dik = (Dik + Djk - Dij) / 2;	// updated disance
				ck.sum_of_dists += Dik;
				ci.sum_of_dists += Dik;

				D[TriangleMatrix::access(ci.row_id, ck.row_id)] = Dik;
			}
		}

		// remove cj
		clusters.erase(clusters.begin() + min_j);
	}

	// join two remanining clusters
	tree.push_back(node_t(clusters[0].node_id, clusters[1].node_id));
}
