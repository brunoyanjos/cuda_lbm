#ifndef BOUNDARIES_CUH
#define BOUNDARIES_CUH

// CUDA INCLUDE
#include <cuda.h>
#include <cuda_runtime.h>
#include "device_launch_parameters.h"

#include "../../var.h"
#include "../../globalFunctions.h"
#include "../../nodeTypeMap.h"
#include "numerical_solutions.cuh"
#include "aux_functions.cuh"
#include "interpolation_utilities.cuh"

__host__ __device__ inline void boundary_definition(unsigned int *nodeType, unsigned int x, unsigned int y)
{
	if (x == 0)
	{
		*nodeType = WEST;
	}
	else if (x == (NX - 1))
	{
		*nodeType = EAST;
	}
	else
	{
		*nodeType = BULK;
	}
}

__device__ inline void boundary_calculation(unsigned int nodeType, dfloat *rhoVar, dfloat *ux, dfloat *uy, dfloat *mxx, dfloat *myy, dfloat *mxy, dfloat *pop, dfloat *fMom, int x, int y, dfloat OMEGA)
{
	const dfloat pop_0 = pop[0] + W0;

	const dfloat pop_1 = pop[1] + W1;
	const dfloat pop_2 = pop[2] + W1;
	const dfloat pop_3 = pop[3] + W1;
	const dfloat pop_4 = pop[4] + W1;

	const dfloat pop_5 = pop[5] + W2;
	const dfloat pop_6 = pop[6] + W2;
	const dfloat pop_7 = pop[7] + W2;
	const dfloat pop_8 = pop[8] + W2;

	switch (nodeType)
	{
	case WEST:
	{
		const dfloat rhoIn = pop_0 + pop_2 + pop_3 + pop_4 + pop_6 + pop_7;
		const dfloat inv_rhoIn = static_cast<float>(1) / rhoIn;

		const dfloat mxxIn = (pop_3 + pop_6 + pop_7) * inv_rhoIn - cs2;
		const dfloat mxyIn = (pop_7 - pop_6) * inv_rhoIn;
		const dfloat myyIn = (pop_2 + pop_4 + pop_6 + pop_7) * inv_rhoIn - cs2;

		*uy = static_cast<float>(0);
		*ux = U_MAX;

		const dfloat rho = (static_cast<float>(4) * rhoIn + static_cast<float>(3) * rhoIn * mxxIn) /
						   (static_cast<float>(3) - static_cast<float>(3) * (*ux));
		*mxx = (rho + static_cast<float>(9) * rhoIn * mxxIn + static_cast<float>(3) * rho * (*ux)) / (static_cast<float>(6) * rho);
		*mxy = static_cast<float>(2) * rhoIn * mxyIn / rho;
		*myy = static_cast<float>(6) * rhoIn * myyIn / (static_cast<float>(5) * rho);

		*rhoVar = rho;

		break;
	}
	case EAST:
	{
		const dfloat rhoIn = pop_0 + pop_1 + pop_2 + pop_4 + pop_5 + pop_8;
		const dfloat inv_rhoIn = static_cast<float>(1) / rhoIn;

		const dfloat mxxIn = (pop_1 + pop_5 + pop_8) * inv_rhoIn - cs2;
		const dfloat mxyIn = (pop_5 - pop_8) * inv_rhoIn;
		const dfloat myyIn = (pop_2 + pop_4 + pop_5 + pop_8) * inv_rhoIn - cs2;

		const dfloat rho = RHO_0 + fMom[idxMom(threadIdx.x - 1, threadIdx.y, M_RHO_INDEX, blockIdx.x, blockIdx.y)];
		*ux = fMom[idxMom(threadIdx.x - 1, threadIdx.y, M_UX_INDEX, blockIdx.x, blockIdx.y)] / F_M_I_SCALE;
		*uy = fMom[idxMom(threadIdx.x - 1, threadIdx.y, M_UY_INDEX, blockIdx.x, blockIdx.y)] / F_M_I_SCALE;

		*rhoVar = rho;

		*mxx = (rho + static_cast<float>(9) * rhoIn * mxxIn - static_cast<float>(3) * rho * (*ux)) / (static_cast<float>(6) * rho);
		*mxy = (static_cast<float>(6) * rhoIn * mxyIn - rho * (*uy)) / (static_cast<float>(3) * rho);
		*myy = static_cast<float>(6) * rhoIn * myyIn / (static_cast<float>(5) * rho);

		break;
	}
	default:
		break;
	}
}

#endif // BOUNDARY_FUNCTIONS_CUH
