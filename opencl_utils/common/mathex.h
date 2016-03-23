#pragma once
#include <cstdint>
#include <cmath>
#include <limits>

#undef min
#undef max

namespace mathex
{

template <typename T>
inline const T log2(T x) 
{
	return std::log(x)/std::log(2);
}

template <typename T>
inline const T max(const T& a, const T& b)
{
	return b > a ? b : a;
}

template <typename T>
inline const T min(const T& a, const T& b)
{
	return a > b ? b : a;
}

template <typename T>
inline const T max(const T& a, const T& b, const T& c)
{
	return max(max(a,b),c);
}

template <typename T>
inline const T min(const T& a, const T& b, const T& c)
{
	return min(min(a,b),c);
}


template <typename T>
inline const T max(const T& a, const T& b, const T& c, const T& d)
{
	return max(max(a,b),max(c,d));
}

template <typename T>
inline const T max(const T& a, const T& b, const T& c, const T& d, const T& e)
{
	return max(max(a,b),max(c,d,e));
}

template <typename T>
inline const T min(const T& a, const T& b, const T& c, const T& d)
{
	return min(min(a,b),min(c,d));
}

template <typename T>
inline const T ceildiv(const T& a, const T& b)
{
	return (a + b - 1) / b;
}

template <typename T>
inline const T ceilround(const T& a, const T& b)
{
	return ((a + b - 1) / b) * b;
}

template <typename T>
inline const T floorround(const T& a, const T& b)
{
	return ((a - b + 1) / b) * b;
}

template <typename T, typename U> 
T max3idx(T a, T b, T c, U &index)
{
	if (a >= b)
		if (a >= c)	{index = 0; return a;}		
		else		{index = 2; return c;}			
	else 								
		if (b >= c)	{index = 1; return b;}		
		else		{index = 2; return c;}	
}

template <typename T, typename U> 
T max4idx(T a, T b, T c, T d, U &index)
{
	if (a >= b)
		if (a >= c)									 
			if (a >= d)	{index = 0; return a;}	
			else		{index = 3; return d;}			
		else /* c > a */							
			if (c >= d)	{index = 2; return c;}	
			else		{index = 3; return d;}			
	else /* b > a */								
		if (b >= c)									 
			if (b >= d)	{index = 1; return b;}	
			else		{index = 3; return d;}			
		else /* c > b */							
			if (c >= d)	{index = 2; return c;}	
			else		{index = 3; return d;}
}


void half2float(float* out, const uint16_t in);
void float2half(uint16_t* out, const float in);

bool almost_eq(float A, float B, int maxUlps);
int diff_ulps(float A, float B);

}