#ifndef BOUNDARIES_CUH
#define BOUNDARIES_CUH

// CUDA INCLUDE
#include <cuda.h>
#include <cuda_runtime.h>
#include "device_launch_parameters.h"

#include "../../var.h"
#include "../../globalFunctions.h"
#include "../../nodeTypeMap.h"

__host__ inline size_t coarse_boundary_definition(size_t x, size_t y)
{
	if (x == 0 && y == 0)
	{
		return SOUTH_WEST;
	}
	else if (x == 0 && y == (NY_COARSE - 1))
	{
		return NORTH_WEST;
	}
	else if (y == 0)
	{
		return SOUTH;
	}
	else if (y == (NY_COARSE - 1))
	{
		return NORTH;
	}
	else if (x == 0)
	{
		return WEST;
	}
	else
	{
		return BULK;
	}
}

__host__ inline size_t fine_boundary_definition(size_t x, size_t y)
{
	if (x == (NX_FINE - 1) && y == 0)
	{
		return SOUTH_EAST;
	}
	else if (x == (NX_FINE - 1) && y == (NY_FINE - 1))
	{
		return NORTH_EAST;
	}
	else if (y == 0)
	{
		return SOUTH;
	}
	else if (y == (NY_FINE - 1))
	{
		return NORTH;
	}
	else if (x == (NX_FINE - 1))
	{
		return EAST;
	}
	else
	{
		return BULK;
	}
}

#endif // BOUNDARY_FUNCTIONS_CUH
