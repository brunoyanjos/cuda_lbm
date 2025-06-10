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
		*nodeType = SOUTH_WEST;

	}
	else if (x == 0 && y == (NY - 1))
	{
		*nodeType = NORTH_WEST;
	}
	else if (x == (NX - 1) && y == 0)
	{
		*nodeType = SOUTH_EAST;
	}
	else if (x == (NX - 1) && y == (NY - 1))
	{
		*nodeType = NORTH_EAST;
	}
	else if (y == 0)
	{
		*nodeType = SOUTH;
	}
	else if (y == (NY - 1))
	{
		*nodeType = NORTH;
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

		*ux = U_MAX;
		*uy = 0.0;

		if (IRBC) {
			*rhoVar = 6.0 * rhoIn / 5.0;

			*mxx = U_MAX * U_MAX;
			*mxy = 5.0 * mxyIn / 3.0 - U_MAX / 3.0;
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
	case WEST: {
		const dfloat rhoIn = pop[0] + pop[2] + pop[3] + pop[4] + pop[6] + pop[7];
		const dfloat inv_rhoIn = 1.0 / rhoIn;

		const dfloat mxxIn = (pop[3] + pop[6] + pop[7]) * inv_rhoIn - cs2;
		const dfloat mxyIn = (pop[7] - pop[6]) * inv_rhoIn;
		const dfloat myyIn = (pop[2] + pop[4] + pop[6] + pop[7]) * inv_rhoIn - cs2;

		*ux = 0.0;
		*uy = 0.0;

		if (IRBC) {
			*rhoVar = 6.0 * rhoIn / 5.0;

			*mxx = 0.0;
			*mxy = 5.0 * mxyIn / 3.0;
			*myy = 0.0;
		}
		else {
			*rhoVar = -3.0 * (-4.0 * rhoIn - 3.0 * mxxIn * rhoIn + 3.0 * mxxIn * OMEGA * rhoIn) / (9.0 + OMEGA);

			*mxx = (-2.0 - 15.0 * mxxIn) / (3.0 * (-4.0 - 3.0 * mxxIn + 3.0 * mxxIn * OMEGA));
			*mxy = -2.0 * mxyIn * (9.0 + OMEGA) / (3.0 * (-4.0 - 3.0 * mxxIn + 3.0 * mxxIn * OMEGA));
			*myy = -2.0 * myyIn * (9.0 + OMEGA) / (5.0 * (-4.0 - 3.0 * mxxIn + 3.0 * mxxIn * OMEGA));
		}

		break;
	}
	case EAST: {
		const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[4] + pop[5] + pop[8];
		const dfloat inv_rhoIn = 1.0 / rhoIn;

		const dfloat mxxIn = (pop[1] + pop[5] + pop[8]) * inv_rhoIn - cs2;
		const dfloat mxyIn = (pop[5] - pop[8]) * inv_rhoIn;
		const dfloat myyIn = (pop[2] + pop[4] + pop[5] + pop[8]) * inv_rhoIn - cs2;

		*ux = 0.0;
		*uy = 0.0;

		if (IRBC) {
			*rhoVar = 6.0 * rhoIn / 5.0;

			*mxx = 0.0;
			*mxy = 5.0 * mxyIn / 3.0;
			*myy = 0.0;
		}
		else {
			*rhoVar = -3.0 * (-4.0 * rhoIn - 3.0 * mxxIn * rhoIn + 3.0 * mxxIn * OMEGA * rhoIn) / (9.0 + OMEGA);

			*mxx = (-2.0 - 15.0 * mxxIn) / (3.0 * (-4.0 - 3.0 * mxxIn + 3.0 * mxxIn * OMEGA));
			*mxy = -2.0 * mxyIn * (9.0 + OMEGA) / (3.0 * (-4.0 - 3.0 * mxxIn + 3.0 * mxxIn * OMEGA));
			*myy = -2.0 * myyIn * (9.0 + OMEGA) / (5.0 * (-4.0 - 3.0 * mxxIn + 3.0 * mxxIn * OMEGA));
		}


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

		*ux = U_MAX;
		*uy = 0.0;

		if (IRBC) {
			*rhoVar = -36.0 * (mxyIn * OMEGA * rhoIn - rhoIn - mxyIn * rhoIn) / (24.0 + OMEGA + 18.0 * U_MAX - 3.0 * OMEGA * U_MAX - 18.0 * U_MAX * U_MAX + 3.0 * OMEGA * U_MAX * U_MAX);

			*mxx = U_MAX * U_MAX;
			*mxy = (36.0 * mxyIn * rhoIn + (*rhoVar) - 3.0 * U_MAX * (*rhoVar) + 3.0 * U_MAX * U_MAX * (*rhoVar)) / (9.0 * (*rhoVar));
			*myy = 0.0;
		}
		else {
			*rhoVar = 12.0 * rhoIn * (-3.0 - 3.0 * mxxIn - 7.0 * mxyIn - 3.0 * myyIn + 3.0 * mxxIn * OMEGA + 7.0 * mxyIn * OMEGA + 3 * myyIn * OMEGA) / (-16.0 - 9.0 * OMEGA - 14.0 * U_MAX - OMEGA * U_MAX + 15.0 * OMEGA * U_MAX * U_MAX);

			*mxx = 2.0 * (9.0 * rhoIn * mxxIn + 6.0 * rhoIn * mxyIn + (*rhoVar) + 2.0 * U_MAX * (*rhoVar)) / (9.0 * (*rhoVar));
			*mxy = -(-18.0 * rhoIn * mxxIn - 132.0 * rhoIn * mxyIn - 18.0 * rhoIn * myyIn - 7.0 * (*rhoVar) + 7.0 * U_MAX * (*rhoVar)) / (27.0 * (*rhoVar));
			*myy = -2.0 * (-6.0 * rhoIn * mxyIn - 9.0 * rhoIn * myyIn - (*rhoVar) + U_MAX * (*rhoVar)) / (9.0 * (*rhoVar));
		}

		break;
	}
	case NORTH_EAST: {
		const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[5];
		const dfloat inv_rhoIn = 1.0 / rhoIn;

		const dfloat mxxIn = (pop[1] + pop[5]) * inv_rhoIn - cs2;
		const dfloat mxyIn = pop[5] * inv_rhoIn;
		const dfloat myyIn = (pop[2] + pop[5]) * inv_rhoIn - cs2;

		*ux = U_MAX;
		*uy = 0.0;

		if (IRBC) {
			*rhoVar = 36.0 * (mxyIn * OMEGA * rhoIn + rhoIn - mxyIn * rhoIn) / (24.0 + OMEGA - 18.0 * U_MAX + 3.0 * OMEGA * U_MAX - 18.0 * U_MAX * U_MAX + 3.0 * OMEGA * U_MAX * U_MAX);

			*mxx = U_MAX * U_MAX;
			*mxy = (36.0 * mxyIn * rhoIn - (*rhoVar) - 3.0 * U_MAX * (*rhoVar) - 3.0 * U_MAX * U_MAX * (*rhoVar)) / (9.0 * (*rhoVar));
			*myy = 0.0;
		}
		else {
			*rhoVar = 12.0 * rhoIn * (-3.0 - 3.0 * mxxIn + 7.0 * mxyIn - 3.0 * myyIn + 3.0 * mxxIn * OMEGA - 7.0 * mxyIn * OMEGA + 3 * myyIn * OMEGA) / (-16.0 - 9.0 * OMEGA + 14.0 * U_MAX + OMEGA * U_MAX + 15.0 * OMEGA * U_MAX * U_MAX);

			*mxx = -2.0 * (-9.0 * rhoIn * mxxIn + 6.0 * rhoIn * mxyIn - (*rhoVar) + 2.0 * U_MAX * (*rhoVar)) / (9.0 * (*rhoVar));
			*mxy = -(18.0 * rhoIn * mxxIn - 132.0 * rhoIn * mxyIn + 18.0 * rhoIn * myyIn + 7.0 * (*rhoVar) + 7.0 * U_MAX * (*rhoVar)) / (27.0 * (*rhoVar));
			*myy = 2.0 * (-6.0 * rhoIn * mxyIn + 9.0 * rhoIn * myyIn + (*rhoVar) + U_MAX * (*rhoVar)) / (9.0 * (*rhoVar));
		}

		break;
	}
	default:
		break;
	}
}

#endif // BOUNDARY_FUNCTIONS_CUH
