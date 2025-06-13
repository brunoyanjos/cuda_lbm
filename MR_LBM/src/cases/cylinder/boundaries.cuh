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
		* nodeType = WEST;
#else
		* nodeType = SOUTH_WEST;
#endif
	}
	else if (x == 0 && y == (NY - 1))
	{
#ifdef BC_Y_PERIODIC
		* nodeType = WEST;
#else
		* nodeType = NORTH_WEST;
#endif
	}
	else if (x == (NX - 1) && y == 0)
	{
#ifdef BC_Y_PERIODIC
		* nodeType = EAST;
#else
		* nodeType = SOUTH_EAST;
#endif
	}
	else if (x == (NX - 1) && y == (NY - 1))
	{
#ifdef BC_Y_PERIODIC
		* nodeType = EAST;
#else
		* nodeType = NORTH_EAST;

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
inline void boundary_calculation(unsigned int nodeType, dfloat* rhoVar, dfloat* ux, dfloat* uy, dfloat* mxx, dfloat* myy, dfloat* mxy, dfloat* pop, dfloat* fMom, int x, int y, dfloat OMEGA) {
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

		*rhoVar = 3.0 * rhoIn * (4.0 + 3.0 * (1.0 - OMEGA) * myyIn) / (9.0 + OMEGA);

		*mxx = 6.0 * rhoIn * mxxIn/(5.0 * (*rhoVar));
		*mxy = 2.0 * rhoIn * mxyIn/ (*rhoVar);
		*myy = (*rhoVar + 9.0 * rhoIn * myyIn)/(6.0 * (*rhoVar));

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

		*rhoVar = 3.0 * rhoIn * (4.0 + 3.0 * (1.0 - OMEGA) * myyIn) / (9.0 + OMEGA);

		*mxx = 6.0 * rhoIn * mxxIn/(5.0 * (*rhoVar));
		*mxy = 2.0 * rhoIn * mxyIn/ (*rhoVar);
		*myy = (*rhoVar + 9.0 * rhoIn * myyIn)/(6.0 * (*rhoVar));

		break;
	}
	case WEST:
	{
		const dfloat rhoIn = pop[0] + pop[2] + pop[3] + pop[4] + pop[6] + pop[7];
		const dfloat inv_rhoIn = 1.0 / rhoIn;

		const dfloat mxxIn = (pop[3] + pop[6] + pop[7]) * inv_rhoIn - cs2;
		const dfloat mxyIn = (pop[7] - pop[6]) * inv_rhoIn;
		const dfloat myyIn = (pop[2] + pop[4] + pop[6] + pop[7]) * inv_rhoIn - cs2;
		
		const dfloat H= NY-1;
		*uy = 0.0;
		*ux = 4.0 * U_MAX *((y/H) - (y/H)*(y/H));
		// *ux = U_MAX;
		
		const dfloat rho = (4.0 * rhoIn + 3.0 * rhoIn * mxxIn) / (3.0 - 3.0 * (*ux));
		*mxx = (rho + 9.0 * rhoIn * mxxIn + 3.0 * rho * (*ux)) / (6.0 * rho);
		*mxy = 2.0 * rhoIn * mxyIn / rho;
		*myy = 6.0 * rhoIn * myyIn / (5.0 * rho);

		*rhoVar = rho;

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

		*rhoVar = RHO_0;

		*mxx = (rho + 9.0 * rhoIn * mxxIn - 3.0 * rho * *ux) / (6.0 * rho);
		*mxy = (6.0 * rhoIn * mxyIn - rho * *uy) / (3.0 * rho);
		*myy = 6.0 * rhoIn * myyIn / (5.0 * rho);

		break;
	}
	
	case SOUTH_WEST: {
		const dfloat rhoIn = pop[0] + pop[3] + pop[4] + pop[7];
		const dfloat inv_rhoIn = 1.0 / rhoIn;

		const dfloat mxxIn = (pop[3] + pop[7]) * inv_rhoIn - cs2;
		const dfloat mxyIn = pop[7] * inv_rhoIn;
		const dfloat myyIn = (pop[4] + pop[7]) * inv_rhoIn - cs2;

		*ux = 0.0;
		*uy = 0.0;

		if (IRBC) {
			*rhoVar = 36.0 * (rhoIn - mxyIn * rhoIn + mxyIn * OMEGA * rhoIn) / (24.0 + OMEGA);

			*mxx = 0.0;
			*mxy = (36.0 * mxyIn * rhoIn - (*rhoVar)) / (9.0 * (*rhoVar));
			*myy = 0.0;
		}
		else {
			*rhoVar = -12.0 * rhoIn * (-3.0 - 3.0 * mxxIn + 7.0 * mxyIn - 3.0 * myyIn + 3.0 * mxxIn * OMEGA - 7.0 * mxyIn * OMEGA + 3 * myyIn * OMEGA) / (16.0 + 9.0 * OMEGA);

			*mxx = 2.0 * (9.0 * rhoIn * mxxIn - 6.0 * rhoIn * mxyIn + (*rhoVar)) / (9.0 * (*rhoVar));
			*mxy = -(18.0 * rhoIn * mxxIn - 132.0 * rhoIn * mxyIn + 18.0 * rhoIn * myyIn + 7.0 * (*rhoVar)) / (27.0 * (*rhoVar));
			*myy = -2.0 * (6.0 * rhoIn * mxyIn - 9.0 * rhoIn * myyIn - (*rhoVar)) / (9.0 * (*rhoVar));
		}

		break;
	}
	case SOUTH_EAST: {
		const dfloat rhoIn = pop[0] + pop[1] + pop[4] + pop[8];
		const dfloat inv_rhoIn = 1.0 / rhoIn;

		const dfloat mxxIn = (pop[1] + pop[8]) * inv_rhoIn - cs2;
		const dfloat mxyIn = -pop[8] * inv_rhoIn;
		const dfloat myyIn = (pop[4] + pop[8]) * inv_rhoIn - cs2;


		*ux = 0.0;
		*uy = 0.0;

		if (IRBC) {
			*rhoVar = -36.0 * (mxyIn * OMEGA * rhoIn - rhoIn - mxyIn * rhoIn) / (24 + OMEGA);

			*mxx = 0.0;
			*mxy = (36.0 * mxyIn * rhoIn + (*rhoVar)) / (9.0 * (*rhoVar));
			*myy = 0.0;
		}
		else {
			*rhoVar = -12.0 * rhoIn * (-3.0 - 3.0 * mxxIn - 7.0 * mxyIn - 3.0 * myyIn + 3.0 * mxxIn * OMEGA + 7.0 * mxyIn * OMEGA + 3 * myyIn * OMEGA) / (16.0 + 9.0 * OMEGA);

			*mxx = 2.0 * (9.0 * rhoIn * mxxIn + 6.0 * rhoIn * mxyIn + (*rhoVar)) / (9.0 * (*rhoVar));
			*mxy = -(-18.0 * rhoIn * mxxIn - 132.0 * rhoIn * mxyIn - 18.0 * rhoIn * myyIn - 7.0 * (*rhoVar)) / (27.0 * (*rhoVar));
			*myy = 2.0 * (6.0 * rhoIn * mxyIn + 9.0 * rhoIn * myyIn + (*rhoVar)) / (9.0 * (*rhoVar));
		}

		break;
	}
	case NORTH_WEST: {
		const dfloat rhoIn = pop[0] + pop[2] + pop[3] + pop[6];
		const dfloat inv_rhoIn = 1.0 / rhoIn;

		const dfloat mxxIn = (pop[3] + pop[6]) * inv_rhoIn - cs2;
		const dfloat mxyIn = -pop[6] * inv_rhoIn;
		const dfloat myyIn = (pop[2] + pop[6]) * inv_rhoIn - cs2;

		*ux = 0.0;
		*uy = 0.0;

		if (IRBC) {
			*rhoVar = -36.0 * (mxyIn * OMEGA * rhoIn - rhoIn - mxyIn * rhoIn) / (24 + OMEGA);

			*mxx = 0.0;
			*mxy = (36.0 * mxyIn * rhoIn + (*rhoVar)) / (9.0 * (*rhoVar));
			*myy = 0.0;
		}
		else {
			*rhoVar = -12.0 * rhoIn * (-3.0 - 3.0 * mxxIn - 7.0 * mxyIn - 3.0 * myyIn + 3.0 * mxxIn * OMEGA + 7.0 * mxyIn * OMEGA + 3 * myyIn * OMEGA) / (16.0 + 9.0 * OMEGA);

			*mxx = 2.0 * (9.0 * rhoIn * mxxIn + 6.0 * rhoIn * mxyIn + (*rhoVar)) / (9.0 * (*rhoVar));
			*mxy = -(-18.0 * rhoIn * mxxIn - 132.0 * rhoIn * mxyIn - 18.0 * rhoIn * myyIn - 7.0 * (*rhoVar)) / (27.0 * (*rhoVar));
			*myy = 2.0 * (6.0 * rhoIn * mxyIn + 9.0 * rhoIn * myyIn + (*rhoVar)) / (9.0 * (*rhoVar));
		}

		break;
	}
	case NORTH_EAST: {
		const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[5];
		const dfloat inv_rhoIn = 1.0 / rhoIn;

		const dfloat mxxIn = (pop[1] + pop[5]) * inv_rhoIn - cs2;
		const dfloat mxyIn = pop[5] * inv_rhoIn;
		const dfloat myyIn = (pop[2] + pop[5]) * inv_rhoIn - cs2;

		*ux = 0.0;
		*uy = 0.0;

		if (IRBC) {
			*rhoVar = 36.0 * (rhoIn - mxyIn * rhoIn + mxyIn * OMEGA * rhoIn) / (24.0 + OMEGA);

			*mxx = 0.0;
			*mxy = (36.0 * mxyIn * rhoIn - (*rhoVar)) / (9.0 * (*rhoVar));
			*myy = 0.0;
		}
		else {
			*rhoVar = -12.0 * rhoIn * (-3.0 - 3.0 * mxxIn + 7.0 * mxyIn - 3.0 * myyIn + 3.0 * mxxIn * OMEGA - 7.0 * mxyIn * OMEGA + 3 * myyIn * OMEGA) / (16.0 + 9.0 * OMEGA);

			*mxx = 2.0 * (9.0 * rhoIn * mxxIn - 6.0 * rhoIn * mxyIn + (*rhoVar)) / (9.0 * (*rhoVar));
			*mxy = -(18.0 * rhoIn * mxxIn - 132.0 * rhoIn * mxyIn + 18.0 * rhoIn * myyIn + 7.0 * (*rhoVar)) / (27.0 * (*rhoVar));
			*myy = -2.0 * (6.0 * rhoIn * mxyIn - 9.0 * rhoIn * myyIn - (*rhoVar)) / (9.0 * (*rhoVar));
		}

		break;
	}
	default:
		break;
	}
}

#endif // BOUNDARY_FUNCTIONS_CUH
