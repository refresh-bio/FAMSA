#include <cmath>
#include <cassert>
#include <cstring>
#include <memory>
#include <iostream>

#include "mathex.h"

int mathex::diff_ulps(float A, float B)
{
	int aInt;
	int bInt;
	std::memcpy(&aInt, &A, sizeof(int));
	std::memcpy(&bInt, &B, sizeof(int));

	// Make aInt lexicographically ordered as a twos-complement int
	if (aInt < 0)
		aInt = 0x80000000 - aInt;

	// Make bInt lexicographically ordered as a twos-complement int
	if (bInt < 0)
		bInt = 0x80000000 - bInt;

	int diff = std::abs(aInt - bInt);

	return diff;
}

bool mathex::almost_eq(float A, float B, int maxUlps)
{
	// Make sure maxUlps is non-negative and small enough that the
	// default NAN won't compare as equal to anything.
	assert(maxUlps >= 0 && maxUlps < 4 * 1024 * 1024);

	if (diff_ulps(A, B) <= maxUlps)
		return true;

	return false;
}




void mathex::half2float(float* out, const uint16_t in) {
	uint32_t t1;
	uint32_t t2;
	uint32_t t3;

	t1 = in & 0x7fff;                       // Non-sign bits
	t2 = in & 0x8000;                       // Sign bit
	t3 = in & 0x7c00;                       // Exponent

	t1 <<= 13;                              // Align mantissa on MSB
	t2 <<= 16;                              // Shift sign bit into position

	t1 += 0x38000000;                       // Adjust bias

	t1 = (t3 == 0 ? 0 : t1);                // Denormals-as-zero

	t1 |= t2;                               // Re-insert sign bit

	*((uint32_t*)out) = t1;
};

void mathex::float2half(uint16_t* out, const float in) {
	uint32_t inu = *((uint32_t*)&in);
	uint32_t t1;
	uint32_t t2;
	uint32_t t3;

	t1 = inu & 0x7fffffff;                 // Non-sign bits
	t2 = inu & 0x80000000;                 // Sign bit
	t3 = inu & 0x7f800000;                 // Exponent

	t1 >>= 13;                             // Align mantissa on MSB
	t2 >>= 16;                             // Shift sign bit into position

	t1 -= 0x1c000;                         // Adjust bias

	t1 = (t3 > 0x38800000) ? 0 : t1;       // Flush-to-zero
	t1 = (t3 < 0x8e000000) ? 0x7bff : t1;  // Clamp-to-max
	t1 = (t3 == 0 ? 0 : t1);               // Denormals-as-zero

	t1 |= t2;                              // Re-insert sign bit

	*((uint16_t*)out) = t1;
};
