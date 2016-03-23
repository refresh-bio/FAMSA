#pragma once
#include <xmmintrin.h>
#include <boost/type_traits.hpp>
#include "mathex.h"

namespace mathsse
{
	// alignment type 
	typedef boost::align::a16 sse_t;
	
	inline void add(float* y, const float* x, const ::size_t n)
	{
		::size_t i = 0;
		::size_t end = mathex::floorround(n, (::size_t)8);
		
		for(i = 0; i < end; i += 8) {
			// load 8 data elements at a time
			__m128 x1 = _mm_loadu_ps(x + i + 0);
			__m128 x2 = _mm_loadu_ps(x + i + 4);
			__m128 y1 = _mm_loadu_ps(y + i + 0);
			__m128 y2 = _mm_loadu_ps(y + i + 4);
			// do the computations
			__m128 z1 = _mm_add_ps(y1, x1);
			__m128 z2 = _mm_add_ps(y2, x2);
			// store the results
			_mm_storeu_ps(y + i + 0, z1);
			_mm_storeu_ps(y + i + 4, z2);
		}

		for(; i < n; i++) {
			y[n] += x[n];
		}
	}
};
