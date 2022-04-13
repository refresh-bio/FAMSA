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

void CLARANS::operator()(const float* distanceMatrix, int n_elems, int n_medoids, int n_fixed_medoids, int* medoids) {

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
	constexpr float INF = std::numeric_limits<float>::max();
	solution_t best{ INF, medoids };
	solution_t current{ INF, new int[n_elems] };
	
	float *raw_distances = new float[(size_t)n_elems * 3];
	float *dists_nearest	= raw_distances + 0 * n_elems;
	float *dists_second		= raw_distances + 1 * n_elems;
	float *deltas = raw_distances + 2 * (size_t)n_elems;
	
	int * raw_assign = new int[(size_t)n_elems * 2];
	int * assign_nearest = raw_assign + 0 * n_elems;
	int * assign_second = raw_assign + 1 * n_elems;

	// select centers randomly
	std::mt19937 gen_nodes;
	std::mt19937 gen_positions;
	
	det_uniform_int_distribution<int> non_medoid_sampling(n_medoids, n_elems - 1);

	// test NUM_LOCAL starting nodes
	for (int iter = 0; iter < numLocal; ++iter) {
	//	LOG_DEBUG << std::endl << "=====================================================" << std::endl;
		// select starting node randomly
		partial_shuffle(candidate + n_fixed_medoids, candidate + n_elems, candidate + n_elems, gen_nodes);
		
		std::copy_n(candidate, n_elems, current.candidate);
		current.cost = std::numeric_limits<float>::max();

		// initialize info for medoids
		for (int mm = 0; mm < n_medoids; ++mm) {
			int m = candidate[mm];
			dists_nearest[m] = 0;	// medoid is at 0 distance from itself
			dists_second[m] = -1;	// dummy value
			assign_nearest[m] = -1;	// dummy value
			assign_second[m] = -1;	// dummy value
		}

		// initialize info for non-medoids
		current.cost = 0;
		for (int xx = n_medoids; xx < n_elems; ++xx) {
			
			int x = candidate[xx];
			updateAssignment(x, candidate, n_medoids, D, dists_nearest[x], dists_second[x], assign_nearest[x], assign_second[x]);
			current.cost += dists_nearest[x];
		}
		LOG_DEBUG << "Attempt " << iter << "[" << candidate[1] << "]: " << current.cost;

		// search partition graph
		for (int step = 0; step < correctedMaxNeighbour; ++step) {
			
	/*		if (Log::getInstance(Log::LEVEL_DEBUG).isEnabled()) {
				std::copy_n(candidate, n_medoids, std::ostream_iterator<int>(std::cerr, ","));
				float realCost = calculateCost(D, candidate, n_elems, n_medoids);
				LOG_DEBUG << " (cuur = " << current.cost  << ", real = " << realCost << ", delta = " << realCost - current.cost << ")" << std::endl;
				getchar();
			}
*/
			// select non-medoid randomly as a candidate for swap
			int xx = non_medoid_sampling(gen_positions);
			int x = candidate[xx];
	//		float dx = dists_nearest[x];
	//		int assign_x = assignments[x];

		//	std::fill_n(deltas, n_medoids, -dx); // gain for making x a new medoid (same for all current medoids)
			// bug in the psuedocode?
			std::fill_n(deltas, n_medoids, 0); // gain for making x a new medoid (same for all current medoids)
			
			// establish which medoid should be swapped with x
			// iterate over other non-medoids y
			for (int yy = n_medoids; yy < n_elems; ++yy) {
				if (yy == xx) {
					continue;
				}
				int y = candidate[yy];
				float dxy = D[TriangleMatrix::access(x, y)];

				// cache y info
				int nn = assign_nearest[y];
				float dn = dists_nearest[y];
				float ds = dists_second[y];
 				//int n = candidate[nn];
				
				// nn is removed so y will be reassigned
				deltas[nn] += std::min(dxy, ds) - dn;
				
				// if y is reassigned to x
				float change = dxy - dn;
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
			int mm_new = (int)(std::min_element(deltas + n_fixed_medoids, deltas + n_medoids) - deltas);
			float delta = deltas[mm_new];

			// if there is an improvement over current medoid
			if (delta < 0) {
				// swap non-medoid with medoid
				std::swap(candidate[mm_new], candidate[xx]);
				int m_new = candidate[mm_new];
				//int x_new = candidate[xx];

				// new medoid no longer contributes to cost
				current.cost -= dists_nearest[m_new];

				// update stats for new medoid
				dists_nearest[m_new] = 0;
				dists_second[m_new] = -1;
				assign_nearest[m_new] = -1;
				assign_second[m_new] = -1;

				// update distance tables for non-medoids
				for (int yy = n_medoids; yy < n_elems; ++yy) {
					int y = candidate[yy];
					float d_new = D[TriangleMatrix::access(m_new, y)];
					float dn = dists_nearest[y];
					int an = assign_nearest[y];

					// former medoid - update everything
					if (yy == xx) {
						updateAssignment(
							y, candidate, n_medoids, D,
							dists_nearest[y], dists_second[y], 
							assign_nearest[y], assign_second[y]);

						current.cost += dists_nearest[y];
						continue;
					}
					
					if (an == mm_new) {
						// previous nearest medoid was removed
						float ds = dists_second[y];

						if (d_new < ds) {
							// new one is closer than second best - just replace the nearest for new one
							dists_nearest[y] = d_new;
							assign_nearest[y] = mm_new;

							current.cost += d_new - dn;
						}
						else {
							// second nearest is the best assignment
							// recompute
							updateAssignment(y, candidate, n_medoids, D, 
								dists_nearest[y], dists_second[y], 
								assign_nearest[y], assign_second[y]);

							current.cost += ds - dn;
						}
					} else if (d_new < dn) {
						// previous nearest wasn't removed
						// and new medoid is closer than previous nearest
						dists_second[y] = dn;
						assign_second[y] = an;
						dists_nearest[y] = d_new;
						assign_nearest[y] = mm_new;

						current.cost += d_new - dn;
					}
					else { 
						// previous nearest wasn't removed
						// new medoid further than previous nearest
						
						// possible change at second best
						float ds = dists_second[y];
						float as = assign_second[y];
						
						if (as != mm_new && d_new < ds) {
							// second best wasn't removed
							// new medoid closer than second best
							dists_second[y] = d_new;
							assign_second[y] = mm_new;
						}
						else {
							// second nearest is the best assignment
							// recompute
							updateAssignment(y, candidate, n_medoids, D,
								dists_nearest[y], dists_second[y],
								assign_nearest[y], assign_second[y]);

						}
					}
				}

				// save current solution
				std::swap(current.candidate[mm_new], current.candidate[xx]);
				//current.cost += delta;
				step = 0;
				
				LOG_DEBUG << " -> " << current.cost;
		//		LOG_DEBUG << "Accept: " << candidate[xx] << " [" << xx << "] -> " << candidate[mm_new] << " [" << mm_new << "]    (" << delta << ")" << std::endl;
			}
			else {
		//		LOG_DEBUG << step << ". Reject: ? -> " << candidate[xx] << " [" << xx << "]     (" << delta << ")" << std::endl;
			}
		}

		// fixme: recalculate the cost of current(it should be ok)
		// current.cost = calculateCost(D, current.candidate, n_elems, n_medoids);
		
		// if current solution is better than the best solution
		if (current.cost < best.cost) {
			best.cost = current.cost;
			std::copy_n(current.candidate, n_medoids, best.candidate);
		}

		LOG_DEBUG << std::endl;
	}

	LOG_DEBUG << "\t\tBEST: " << best.cost << std::endl;

	delete[] candidate;
	delete[] current.candidate;
	delete[] raw_distances;
	delete[] raw_assign;

}


void CLARANS::updateAssignment(
	int x,
	int *candidate,
	int n_medoids,
	const float* D,
	float& dist_nearest,
	float& dist_second,
	int& assign_nearest,
	int& assign_second) 
{
	float dn = std::numeric_limits<float>::max();
	float ds = std::numeric_limits<float>::max();
	int an = -1;
	int as = -1;

	// find two closest medoids
	for (int mm = 0; mm < n_medoids; ++mm) {
		int m = candidate[mm];
		float d = D[TriangleMatrix::access(m, x)];
		if (d < dn) {
			ds = dn;
			as = an;
			dn = d;
			an = mm;
		}
		else if (d < ds) {
			ds = d;
			as = mm;
		}
	}

	dist_nearest = dn;
	dist_second = ds;
	assign_nearest = an;
	assign_second = as;
}


float CLARANS::calculateCost(const float* distanceMatrix, int *candidate, int n_elems, int n_medoids) {
	
	float cost = 0.0;
	
	// iterate over non-medoids
	for (int xx = n_medoids; xx < n_elems; ++xx) {
		int x = candidate[xx];

		// find closest medoid
		float dmin = std::numeric_limits<float>::max();
		for (int mm = 0; mm < n_medoids; ++mm) {
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