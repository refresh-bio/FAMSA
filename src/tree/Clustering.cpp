#include "Clustering.h"
#include "TreeDefs.h"

#include "../utils/deterministic_random.h"
#include "../utils/log.h"

#include <algorithm>
#include <numeric>
#include <cassert>
#include <iterator>

struct solution_t {
	float cost;
	int* candidate;
};

void CLARANS::operator()(const float* distanceMatrix, size_t n_elems, size_t n_medoids, int* medoids) {

	// CLARANS paper formula
	int n_swaps = (n_elems - n_medoids) * n_medoids;
	int maxNeighbor = n_swaps < minMaxNeighbor
		? n_swaps
		: std::max((int)(exploreFraction * n_swaps), minMaxNeighbor);
	
	// FastCLARANS paper correction 
	// (FastPAM always investigates all possible medoids for swap, so we need to decrease a number of steps)
	int correctedMaxNeighbour = maxNeighbor / n_medoids;

	const float* D = distanceMatrix;

	// array of indices:
	// medoids - first n_centers elements,
	// non-medoids - rest of elements
	int * candidate = new int[n_elems];

	std::iota(candidate, candidate + n_elems, 0);

	// structures for storing solutions
	float INF = std::numeric_limits<float>::max();
	solution_t best{ INF, medoids };
	solution_t current{ INF, new int[n_elems] };
	
	float *raw_distances = new float[n_elems * 3];
	float *dists_nearest	= raw_distances + 0 * n_elems;
	float *dists_second		= raw_distances + 1 * n_elems;
	float *deltas = raw_distances + 2 * n_elems;
	
	int * assignments = new int[n_elems];

	// select centers randomly
	std::mt19937 gen_nodes;
	std::mt19937 gen_positions;
	
	det_uniform_int_distribution<int> non_medoid_sampling(n_medoids, n_elems - 1);

	// test NUM_LOCAL starting nodes
	for (int iter = 0; iter < numLocal; ++iter) {
		LOG_DEBUG << std::endl << "=====================================================" << std::endl;
		// select starting node randomly
		std::shuffle(candidate, candidate + n_elems, gen_nodes);
		
		std::copy_n(candidate, n_elems, current.candidate);
		current.cost = std::numeric_limits<float>::max();

		// initialize info for medoids
		for (size_t mm = 0; mm < n_medoids; ++mm) {
			int m = candidate[mm];
			dists_nearest[m] = 0;	// medoid is at 0 distance from itself
			dists_second[m] = -1;	// dummy value
			assignments[m] = -1;	// dummy value
		}

		// initialize info for non-medoids
		current.cost = 0;
		for (size_t xx = n_medoids; xx < n_elems; ++xx) {
			int x = candidate[xx];

			float dn = INF;
			float ds = INF;
			int assign = -1;
			
			// find two closest medoids
			for (size_t mm = 0; mm < n_medoids; ++mm) {
				int m = candidate[mm];
				float d = D[TriangleMatrix::access(m, x)];
				if (d < dn) {
					ds = dn;
					dn = d;
					assign = mm;
				}
				else if (d < ds) {
					ds = d;
				}
			}
			current.cost += dn;

			dists_nearest[x] = dn;
			dists_second[x] = ds;
			assignments[x] = assign;
		}

		// search partition graph
		for (int step = 0; step < correctedMaxNeighbour; ++step) {
			
			if (Log::getInstance(Log::LEVEL_DEBUG).isEnabled()) {
				std::copy_n(candidate, n_medoids, std::ostream_iterator<int>(std::cerr, ","));
				float realCost = calculateCost(D, candidate, n_elems, n_medoids);
				LOG_DEBUG << " (" << current.cost  << ", " << realCost << ", delta = " << realCost - current.cost << ")" << std::endl;
			}

			// select non-medoid randomly as a candidate for swap
			int xx = non_medoid_sampling(gen_positions);
			int x = candidate[xx];
			float dx = dists_nearest[x];
	//		int assign_x = assignments[x];

	//		std::fill_n(deltas, n_medoids, -dx); // gain for making x a new medoid (same for all current medoids)
			// bug in the psuedocode?
			std::fill_n(deltas, n_medoids, 0); // gain for making x a new medoid (same for all current medoids)
			
			// iterate over other non-medoids y to establish which medoid should be swapped with x
			for (int yy = n_medoids; yy < n_elems; ++yy) {
				if (yy == xx) {
					continue;
				}
				int y = candidate[yy];
				float dxy = D[TriangleMatrix::access(x, y)];

				// cache y info
				int nn = assignments[y];
				float dn = dists_nearest[y];
				float ds = dists_second[y];
 				int n = candidate[nn];
				
				deltas[nn] += std::min(dxy, ds) - dn; // cost/gain of reassigning y to x or to its second best
				float change = dxy - dn;
				// if y will be reassigned to x (dxy < dn < ds)
				if (change < 0) { 
					// iterate over other medoids
					for (int kk = 0; kk < nn; ++kk) {
						deltas[kk] += change; // gain for swapping x with y
					}
					for (int kk = nn + 1; kk < n_medoids; ++kk) {
						deltas[kk] += change; // gain for swapping x with y
					}
				}
			}
			// get id of best medoid to swap (with smallest delta)
			auto mm_new = std::min_element(deltas, deltas + n_medoids) - deltas;
			float delta = deltas[mm_new];

			// if there is an improvement over current medoid
			if (delta < 0) {
				// swap non-medoid with medoid
				std::swap(candidate[mm_new], candidate[xx]);
				int m_new = candidate[mm_new];
				int x_new = candidate[xx];

				// update stats for new medoid
				dists_nearest[m_new] = 0;
				dists_second[m_new] = -1;
				assignments[m_new] = -1;

				// update stats for former medoid
				// (find nearest and second nearest)
				float dn_new = INF;
				float ds_new = INF;
				int assign_new = -1;

				for (size_t kk = 0; kk < n_medoids; ++kk) {
					int k = candidate[kk];
					float d = D[TriangleMatrix::access(k, x_new)];
					if (d < dn_new) {
						ds_new = dn_new;
						dn_new = d;
						assign_new = kk;
					}
					else if (d < ds_new) {
						ds_new = d;
					}
				}

				dists_nearest[x_new] = dn_new;
				dists_second[x_new] = ds_new;
				assignments[x_new] = assign_new;
				
				// update distance tables for other non-medoids
				// (check if new medoid is closer than previously assigned)
				for (size_t yy = n_medoids; yy < n_elems; ++yy) {
					int y = candidate[yy];
					float d = D[TriangleMatrix::access(m_new, y)];
					if (d < dists_nearest[y]) {
						dists_second[y] = dists_nearest[y];
						dists_nearest[y] = d;
						assignments[y] = mm_new;
					}
					else if (d < dists_second[y]) {
						dists_second[y] = d;
					}
				}

				// save current solution
				std::swap(current.candidate[mm_new], current.candidate[xx]);
				current.cost += delta;
				step = 0;

				//assert(current.cost > 0);
				LOG_DEBUG << "Accept: " << candidate[xx] << " [" << xx << "] -> " << candidate[mm_new] << " [" << mm_new << "]    (" << delta << ")" << std::endl;
			}
			else {
				LOG_DEBUG << step << ". Reject: ? -> " << candidate[xx] << " [" << xx << "]     (" << delta << ")" << std::endl;
			}
		}

		// fixme: recalculate the cost of current(it should be ok)
		current.cost = calculateCost(D, current.candidate, n_elems, n_medoids);

		// if current solution is better than the best solution
		if (current.cost < best.cost) {
			best.cost = current.cost;
			std::copy_n(current.candidate, n_medoids, best.candidate);
		}
	}

	delete[] candidate;
	delete[] current.candidate;
	delete[] raw_distances;
	delete[] assignments;

}



float CLARANS::calculateCost(const float* distanceMatrix, int *candidate, size_t n_elems, size_t n_medoids) {
	
	float cost = 0.0;
	
	// iterate over non-medoids
	for (size_t xx = n_medoids; xx < n_elems; ++xx) {
		int x = candidate[xx];

		// find closest medoid
		float dmin = std::numeric_limits<float>::max();
		for (size_t mm = 0; mm < n_medoids; ++mm) {
			int m = candidate[mm];
			float d = distanceMatrix[TriangleMatrix::access(m, x)];
			if (d < dmin) {
				dmin = d;
			}
		}

		cost += dmin;
	}

	return cost;
}