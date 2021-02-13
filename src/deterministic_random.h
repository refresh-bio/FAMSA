#pragma once
#include <random>
#include <tuple>
#include <iostream>

#undef min
#undef max

template<class IntType = int>
class det_uniform_int_distribution {
public:
	// types
	typedef IntType result_type;
	typedef std::pair<int, int> param_type;

	// constructors and reset functions
	explicit det_uniform_int_distribution(IntType a = 0, IntType b = std::numeric_limits<IntType>::max());
	explicit det_uniform_int_distribution(const param_type& parm);
	void reset();

	// generating functions
	template<class URNG>
	result_type operator()(URNG& g);
	template<class URNG>
	result_type operator()(URNG& g, const param_type& parm);

	// property functions
	result_type a() const;
	result_type b() const;
	param_type param() const;
	void param(const param_type& parm);
	result_type min() const;
	result_type max() const;

private:
	typedef typename std::make_unsigned<IntType>::type diff_type;

	IntType lower;
	IntType upper;
};

template<class IntType>
det_uniform_int_distribution<IntType>::det_uniform_int_distribution(IntType a, IntType b) {
	lower = a;
	upper = b;
}

template<class IntType>
det_uniform_int_distribution<IntType>::det_uniform_int_distribution(const param_type& parm) {
	param(parm);
}

template<class IntType>
void det_uniform_int_distribution<IntType>::reset() {}

template<class IntType>
template<class URNG>
auto det_uniform_int_distribution<IntType>::operator()(URNG& g) -> result_type {
	return operator()(g, param());
}

template<class IntType>
template<class URNG>
auto det_uniform_int_distribution<IntType>::operator()(URNG& g, const param_type& parm) -> result_type {
	diff_type diff = (diff_type)parm.second - (diff_type)parm.first + 1;
	if (diff == 0) // If the +1 overflows we are using the full range, just return g()
		return g();

	diff_type badDistLimit = std::numeric_limits<diff_type>::max() / diff;
	do {
		diff_type generatedRand = g();

		if (generatedRand / diff < badDistLimit)
			return (IntType)((generatedRand % diff) + (diff_type)parm.first);
	} while (true);
}

template<class IntType>
auto det_uniform_int_distribution<IntType>::a() const -> result_type {
	return lower;
}

template<class IntType>
auto det_uniform_int_distribution<IntType>::b() const -> result_type {
	return upper;
}

template<class IntType>
auto det_uniform_int_distribution<IntType>::param() const -> param_type {
	return param_type(lower, upper);
}

template<class IntType>
void det_uniform_int_distribution<IntType>::param(const param_type& parm) {
	std::tie(lower, upper) = parm;
	if (upper < lower)
		throw std::exception();
}

template<class IntType>
auto det_uniform_int_distribution<IntType>::min() const -> result_type {
	return lower;
}

template<class IntType>
auto det_uniform_int_distribution<IntType>::max() const -> result_type {
	return upper;
};




template<class RandomIt, class URBG>
void partial_shuffle(RandomIt first, RandomIt middle, RandomIt last, URBG&& g)
{
	typedef typename std::iterator_traits<RandomIt>::difference_type diff_t;
	typedef det_uniform_int_distribution<diff_t> distr_t;
	typedef typename distr_t::param_type param_t;

	distr_t D;
	diff_t n = middle - first;
	diff_t N = last - first - 1;
	for (diff_t i =0; i < n; ++i) {
		using std::swap;
		swap(first[i], first[D(g, param_t(i, N))]);
	}
}


