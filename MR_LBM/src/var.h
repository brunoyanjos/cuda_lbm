#ifndef __VAR_H
#define __VAR_H

#include <builtin_types.h> // for devices variables
#include <stdint.h>		   // for uint32_t
#include <map>

#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include <limits>
#include <cstring>

/* ----------------------------- PROBLEM DEFINE ---------------------------- */
typedef float dfloat;

#define GPU_INDEX 0
/* -------------------------------  TEMPLATE ------------------------------- */

template <typename T>
__host__ __device__ inline T dsqrt(T x);

template <>
__host__ __device__ inline float dsqrt<float>(float x)
{
	return sqrtf(x);
}

template <>
__host__ __device__ inline double dsqrt<double>(double x)
{
	return sqrt(x);
}

template <typename T>
__host__ __device__ inline T dabs(T x);

template <>
__host__ __device__ inline float dabs<float>(float x)
{
	return fabsf(x);
}

template <>
__host__ __device__ inline double dabs<double>(double x)
{
	return fabs(x);
}

/* --------------------------  SIMULATION DEFINES -------------------------- */

#define STR_IMPL(A) #A
#define STR(A) STR_IMPL(A)

// Helper function to get the closest power of 2 under or equal to `n`
constexpr size_t closestPowerOfTwo(size_t n)
{
	size_t power = 1;
	while (power * 2 <= n)
	{
		power *= 2;
	}
	return power;
}

struct BlockDim
{
	size_t x, y;
};

// Compute optimal block dimensions
constexpr BlockDim findOptimalBlockDimensions(size_t maxElements)
{
	size_t bestX = 1, bestY = 1;
	size_t closestVolume = 0;
	float bestForm = 0.0;

	// Iterate over powers of 2 up to `maxElements` to find optimal dimensions
	for (size_t x = closestPowerOfTwo(maxElements); x >= 1; x /= 2)
	{
		for (size_t y = closestPowerOfTwo(maxElements / x); y >= 1; y /= 2)
		{
			if (x * y <= maxElements)
			{
				const size_t volume = x * y;
				float form = 1.0 / (1.0 / x + 1.0 / y);
				if (volume > closestVolume)
				{
					bestX = x;
					bestY = y;
					closestVolume = volume;
					bestForm = form;
				}
				else if (volume == closestVolume && form > bestForm)
				{
					bestX = x;
					bestY = y;
					bestForm = form;
				}
			}
		}
	}
	return {bestX, bestY};
}

#include "cases/outputs_and_model.h"
#include "definitions.h"
#endif //__VAR_H