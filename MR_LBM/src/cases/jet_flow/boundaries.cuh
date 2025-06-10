#ifndef BOUNDARIES_CUH
#define BOUNDARIES_CUH

// CUDA INCLUDE
#include <cuda.h>
#include <cuda_runtime.h>
#include "device_launch_parameters.h"

#include "../../var.h"
#include "../../globalFunctions.h"
#include "../../nodeTypeMap.h"

__host__ __device__
inline void boundary_definition(unsigned int* nodeType, unsigned int x, unsigned int y) {
	if (x == 0 && y == 0)
	{
#ifdef BC_Y_PERIODIC
		* nodeType = BULK;
#else
		* nodeType = SOUTH;
#endif
	}
	else if (x == 0 && y == (NY - 1))
	{
#ifdef BC_Y_PERIODIC
		* nodeType = BULK;
#else
		* nodeType = NORTH;
#endif
	}
	else if (x == (NX - 1) && y == 0)
	{
#ifdef BC_Y_PERIODIC
		* nodeType = BULK;
#else
		* nodeType = SOUTH;
#endif
	}
	else if (x == (NX - 1) && y == (NY - 1))
	{
#ifdef BC_Y_PERIODIC
		* nodeType = BULK;
#else
		* nodeType = NORTH;

#endif
	}
	else if (y == 0)
	{
#ifdef BC_Y_PERIODIC
		* nodeType = BULK;
#else
		* nodeType = SOUTH;
#endif
	}
	else if (y == (NY - 1))
	{
#ifdef BC_Y_PERIODIC
		* nodeType = BULK;
#else
		* nodeType = NORTH;
#endif
	}
	else if (x == 0)
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

__device__
inline void boundary_calculation(unsigned int nodeType, dfloat* rhoVar, dfloat* ux, dfloat* uy, dfloat* mxx, dfloat* myy, dfloat* mxy, dfloat* pop, dfloat* fMom, int x, int y) {
	switch (nodeType)
	{
	case NORTH: {
		const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[3] + pop[5] + pop[6];
		const dfloat inv_rhoIn = 1.0 / rhoIn;

		const dfloat mxxIn = (pop[1] + pop[3] + pop[5] + pop[6]) * inv_rhoIn - cs2;
		const dfloat mxyIn = (pop[5] - pop[6]) * inv_rhoIn;
		const dfloat myyIn = (pop[2] + pop[5] + pop[6]) * inv_rhoIn - cs2;

		*ux = 0.0;
		*uy = 0.0;

		if (IRBC) {
			*rhoVar = 6.0 * rhoIn / 5.0;

			*mxx = 0.0;
			*mxy = 5.0 * mxyIn / 3.0;
			*myy = 0.0;
		}
		else {
			*rhoVar = -3.0 * (-4.0 * rhoIn - 3.0 * myyIn * rhoIn + 3.0 * myyIn * OMEGA * rhoIn) / (9.0 + OMEGA);

			*mxx = -2.0 * mxxIn * (9.0 + OMEGA) / (5.0 * (-4.0 - 3.0 * myyIn + 3.0 * myyIn * OMEGA));
			*mxy = (-18.0 * mxyIn - 2.0 * mxyIn * OMEGA + 4.0 * U_MAX + 3.0 * myyIn * U_MAX - 3.0 * myyIn * OMEGA * U_MAX) / (3.0 * (-4.0 - 3.0 * myyIn + 3.0 * myyIn * OMEGA));
			*myy = (-2.0 - 15.0 * myyIn) / (3.0 * (-4.0 - 3.0 * myyIn + 3.0 * myyIn * OMEGA));
		}

		break;
	}
	case SOUTH: {
		const dfloat rhoIn = pop[0] + pop[1] + pop[3] + pop[4] + pop[7] + pop[8];
		const dfloat inv_rhoIn = 1.0 / rhoIn;

		const dfloat mxxIn = (pop[1] + pop[3] + pop[7] + pop[8]) * inv_rhoIn - cs2;
		const dfloat mxyIn = (pop[7] - pop[8]) * inv_rhoIn;
		const dfloat myyIn = (pop[3] + pop[7] + pop[8]) * inv_rhoIn - cs2;

		*ux = 0.0;
		*uy = 0.0;

		if (IRBC) {
			*rhoVar = 6.0 * rhoIn / 5.0;

			*mxx = 0.0;
			*mxy = 5.0 * mxyIn / 3.0;
			*myy = 0.0;
		}
		else {
			*rhoVar = -3.0 * (-4.0 * rhoIn - 3.0 * myyIn * rhoIn + 3.0 * myyIn * OMEGA * rhoIn) / (9.0 + OMEGA);

			*mxx = -2.0 * mxxIn * (9.0 + OMEGA) / (5.0 * (-4.0 - 3.0 * myyIn + 3.0 * myyIn * OMEGA));
			*mxy = -2.0 * mxyIn * (9.0 + OMEGA) / (3.0 * (-4.0 - 3.0 * myyIn + 3.0 * myyIn * OMEGA));
			*myy = (-2.0 - 15.0 * myyIn) / (3.0 * (-4.0 - 3.0 * myyIn + 3.0 * myyIn * OMEGA));
		}

		break;
	}
	case WEST:
	{
		const dfloat rhoIn = pop[0] + pop[2] + pop[3] + pop[4] + pop[6] + pop[7];
		const dfloat inv_rhoIn = 1.0 / rhoIn;

		const dfloat mxxIn = (pop[3] + pop[6] + pop[7]) * inv_rhoIn - cs2;
		const dfloat mxyIn = (pop[7] - pop[6]) * inv_rhoIn;
		const dfloat myyIn = (pop[2] + pop[4] + pop[6] + pop[7]) * inv_rhoIn - cs2;

		*ux = 0.0;
		*uy = 0.0;

		if (y < NY / 2 + D / 2 && y > NY / 2 - D / 2) {
			*ux = U_MAX;

			const dfloat rho = (4.0 * rhoIn + 3.0 * rhoIn * mxxIn) / (3.0 - 3.0 * U_MAX);

			*mxx = (rho + 9.0 * rhoIn * mxxIn + 3.0 * rho * U_MAX) / (6.0 * rho);
			*mxy = 2.0 * rhoIn * mxyIn / rho;
			*myy = 6.0 * rhoIn * myyIn / (5.0 * rho);

			*rhoVar = rho;
		}
		else {
			*rhoVar = 6.0 * rhoIn / 5.0;

			*mxx = 0.0;
			*mxy = 5.0 * mxyIn / 3.0;
			*myy = 0.0;
		}

		break;
	}
	case EAST:
	{
		const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[4] + pop[5] + pop[8];
		const dfloat inv_rhoIn = 1.0 / rhoIn;

		const dfloat mxxIn = (pop[1] + pop[5] + pop[8]) * inv_rhoIn - cs2;
		const dfloat mxyIn = (pop[5] - pop[8]) * inv_rhoIn;
		const dfloat myyIn = (pop[2] + pop[4] + pop[5] + pop[8]) * inv_rhoIn - cs2;

		*ux = 0.0;
		*uy = 0.0;

		const dfloat rho = RHO_0 + fMom[idxMom(threadIdx.x - 1, threadIdx.y, M_RHO_INDEX, blockIdx.x, blockIdx.y)];
		*ux = fMom[idxMom(threadIdx.x - 1, threadIdx.y, M_UX_INDEX, blockIdx.x, blockIdx.y)] / F_M_I_SCALE;
		*uy = fMom[idxMom(threadIdx.x - 1, threadIdx.y, M_UY_INDEX, blockIdx.x, blockIdx.y)] / F_M_I_SCALE;

		*rhoVar = rho;

		*mxx = (rho + 9.0 * rhoIn * mxxIn - 3.0 * rho * *ux) / (6.0 * rho);
		*mxy = (6.0 * rhoIn * mxyIn - rho * *uy) / (3.0 * rho);
		*myy = 6.0 * rhoIn * myyIn / (5.0 * rho);

		break;
	}
	default:
		break;
	}
}

#endif // BOUNDARY_FUNCTIONS_CUH
