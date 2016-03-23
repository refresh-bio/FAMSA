#pragma once
#include <algorithm>

template <class InputIt, class OutputIt, class BinaryOperation>
void rank_range(InputIt first1, InputIt last1, OutputIt first2, BinaryOperation op) {
	
	std::vector<InputIt> ptrs(last1 - first1);
	int j = 0;
	for (auto it = first1; it != last1; ++it, ++j) {
		ptrs[j] = it;
	}

	std::stable_sort(ptrs.begin(), ptrs.end(), 
		[](InputIt x, InputIt y) -> bool {
			return *x < *y;
	});

	for (int i = 0; i < ptrs.size(); ++i) {
		auto d = ptrs[i] - first1;
		first2[d] = i;
	}
}