#pragma once

#include <cstddef>

class IClustering {
public:
	virtual void operator()(const float* distanceMatrix, size_t n_elems, size_t n_medoids, size_t n_fixed_medoids, int* centers) = 0;
};

class CLARANS : public IClustering {
	
	const int minMaxNeighbor = 250;
	const float exploreFraction;
	const int numLocal;
		
public:

	CLARANS(float exploreFraction, int numLocal) : exploreFraction(exploreFraction), numLocal(numLocal) {}

	void operator()(const float* distanceMatrix, size_t n_elems, size_t n_medoids, size_t n_fixed_medoids, int* medoids) override;

protected:

	float calculateCost(const float* distanceMatrix, int *candidate, size_t n_elems, size_t n_medoids);

	void updateAssignment(
		int x,
		int *candidate,
		size_t n_medoids,
		const float* D,
		float& dist_nearest,
		float& dist_second,
		int& assign_nearest,
		int& assign_second);

};

