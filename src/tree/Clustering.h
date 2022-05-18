#pragma once

#include <cstddef>

class IClustering {
public:
	virtual void operator()(
		const float* distanceMatrix, 
		int n_elems, 
		int n_medoids, 
		int n_fixed_medoids, 
		int* centers) = 0;

	virtual ~IClustering() {}
};

class CLARANS : public IClustering {
	
	const int minMaxNeighbor = 250;
	const float exploreFraction;
	const int numLocal;
		
public:

	CLARANS(float exploreFraction, int numLocal) : exploreFraction(exploreFraction), numLocal(numLocal) {}

	void operator()(const float* distanceMatrix, int n_elems, int n_medoids, int n_fixed_medoids, int* medoids) override;

protected:

	float calculateCost(const float* distanceMatrix, int *candidate, int n_elems, int n_medoids);

	void updateAssignment(
		int x,
		int *candidate,
		int n_medoids,
		const float* D,
		float& dist_nearest,
		float& dist_second,
		int& assign_nearest,
		int& assign_second);

};

