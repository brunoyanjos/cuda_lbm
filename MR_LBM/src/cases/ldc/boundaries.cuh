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
inline void boundary_calculation(unsigned int nodeType, dfloat* rhoVar, dfloat* ux, dfloat* uy, dfloat* mxx, dfloat* myy, dfloat* mxy, dfloat* pop, latticeNode* nodes, int x, int y) {
	switch (nodeType)
	{
	case NORTH: {
		const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[3] + pop[5] + pop[6];
		const dfloat inv_rhoIn = 1.0f / rhoIn;

		const dfloat mxxIn = (pop[1] + pop[3] + pop[5] + pop[6]) * inv_rhoIn - cs2;
		const dfloat mxyIn = (pop[5] - pop[6]) * inv_rhoIn;
		const dfloat myyIn = (pop[2] + pop[5] + pop[6]) * inv_rhoIn - cs2;

		*ux = U_MAX;
		*uy = 0.0f;

		if (IRBC) {
			*rhoVar = 6.0f * rhoIn / 5.0f;

			*mxx = U_MAX * U_MAX;
			*mxy = 5.0f * mxyIn / 3.0f - U_MAX / 3.0f;
			*myy = 0.0f;
		}
		else {
			*rhoVar = -3.0f * (-4.0f * rhoIn - 3.0f * myyIn * rhoIn + 3.0f * myyIn * OMEGA * rhoIn) / (9.0f + OMEGA);

			*mxx = -2.0f * mxxIn * (9.0f + OMEGA) / (5.0f * (-4.0f - 3.0f * myyIn + 3.0f * myyIn * OMEGA));
			*mxy = (-18.0f * mxyIn - 2.0f * mxyIn * OMEGA + 4.0f * U_MAX + 3.0f * myyIn * U_MAX - 3.0f * myyIn * OMEGA * U_MAX) / (3.0f * (-4.0f - 3.0f * myyIn + 3.0f * myyIn * OMEGA));
			*myy = (-2.0f - 15.0f * myyIn) / (3.0f * (-4.0f - 3.0f * myyIn + 3.0f * myyIn * OMEGA));
		}

		break;
	}
	case SOUTH: {
		const dfloat rhoIn = pop[0] + pop[1] + pop[3] + pop[4] + pop[7] + pop[8];
		const dfloat inv_rhoIn = 1.0f / rhoIn;

		const dfloat mxxIn = (pop[1] + pop[3] + pop[7] + pop[8]) * inv_rhoIn - cs2;
		const dfloat mxyIn = (pop[7] - pop[8]) * inv_rhoIn;
		const dfloat myyIn = (pop[3] + pop[7] + pop[8]) * inv_rhoIn - cs2;

		*ux = 0.0f;
		*uy = 0.0f;

		if (IRBC) {
			*rhoVar = 6.0f * rhoIn / 5.0f;

			*mxx = 0.0f;
			*mxy = 5.0f * mxyIn / 3.0f;
			*myy = 0.0f;
		}
		else {
			*rhoVar = -3.0f * (-4.0f * rhoIn - 3.0f * myyIn * rhoIn + 3.0f * myyIn * OMEGA * rhoIn) / (9.0f + OMEGA);

			*mxx = -2.0f * mxxIn * (9.0f + OMEGA) / (5.0f * (-4.0f - 3.0f * myyIn + 3.0f * myyIn * OMEGA));
			*mxy = -2.0f * mxyIn * (9.0f + OMEGA) / (3.0f * (-4.0f - 3.0f * myyIn + 3.0f * myyIn * OMEGA));
			*myy = (-2.0f - 15.0f * myyIn) / (3.0f * (-4.0f - 3.0f * myyIn + 3.0f * myyIn * OMEGA));
		}

		break;
	}
	case WEST: {
		const dfloat rhoIn = pop[0] + pop[2] + pop[3] + pop[4] + pop[6] + pop[7];
		const dfloat inv_rhoIn = 1.0f / rhoIn;

		const dfloat mxxIn = (pop[3] + pop[6] + pop[7]) * inv_rhoIn - cs2;
		const dfloat mxyIn = (pop[7] - pop[6]) * inv_rhoIn;
		const dfloat myyIn = (pop[2] + pop[4] + pop[6] + pop[7]) * inv_rhoIn - cs2;

		*ux = 0.0f;
		*uy = 0.0f;

		if (IRBC) {
			*rhoVar = 6.0f * rhoIn / 5.0f;

			*mxx = 0.0f;
			*mxy = 5.0f * mxyIn / 3.0f;
			*myy = 0.0f;
		}
		else {
			*rhoVar = -3.0f * (-4.0f * rhoIn - 3.0f * mxxIn * rhoIn + 3.0f * mxxIn * OMEGA * rhoIn) / (9.0f + OMEGA);

			*mxx = (-2.0f - 15.0f * mxxIn) / (3.0f * (-4.0f - 3.0f * mxxIn + 3.0f * mxxIn * OMEGA));
			*mxy = -2.0f * mxyIn * (9.0f + OMEGA) / (3.0f * (-4.0f - 3.0f * mxxIn + 3.0f * mxxIn * OMEGA));
			*myy = -2.0f * myyIn * (9.0f + OMEGA) / (5.0f * (-4.0f - 3.0f * mxxIn + 3.0f * mxxIn * OMEGA));
		}

		break;
	}
	case EAST: {
		const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[4] + pop[5] + pop[8];
		const dfloat inv_rhoIn = 1.0f / rhoIn;

		const dfloat mxxIn = (pop[1] + pop[5] + pop[8]) * inv_rhoIn - cs2;
		const dfloat mxyIn = (pop[5] - pop[8]) * inv_rhoIn;
		const dfloat myyIn = (pop[2] + pop[4] + pop[5] + pop[8]) * inv_rhoIn - cs2;

		*ux = 0.0f;
		*uy = 0.0f;

		if (IRBC) {
			*rhoVar = 6.0f * rhoIn / 5.0f;

			*mxx = 0.0f;
			*mxy = 5.0f * mxyIn / 3.0f;
			*myy = 0.0f;
		}
		else {
			*rhoVar = -3.0f * (-4.0f * rhoIn - 3.0f * mxxIn * rhoIn + 3.0f * mxxIn * OMEGA * rhoIn) / (9.0f + OMEGA);

			*mxx = (-2.0f - 15.0f * mxxIn) / (3.0f * (-4.0f - 3.0f * mxxIn + 3.0f * mxxIn * OMEGA));
			*mxy = -2.0f * mxyIn * (9.0f + OMEGA) / (3.0f * (-4.0f - 3.0f * mxxIn + 3.0f * mxxIn * OMEGA));
			*myy = -2.0f * myyIn * (9.0f + OMEGA) / (5.0f * (-4.0f - 3.0f * mxxIn + 3.0f * mxxIn * OMEGA));
		}


		break;
	}
	case SOUTH_WEST: {
		const dfloat rhoIn = pop[0] + pop[3] + pop[4] + pop[7];
		const dfloat inv_rhoIn = 1.0f / rhoIn;

		const dfloat mxxIn = (pop[3] + pop[7]) * inv_rhoIn - cs2;
		const dfloat mxyIn = pop[7] * inv_rhoIn;
		const dfloat myyIn = (pop[4] + pop[7]) * inv_rhoIn - cs2;

		*ux = 0.0f;
		*uy = 0.0f;

		if (IRBC) {
			*rhoVar = 36.0f * (rhoIn - mxyIn * rhoIn + mxyIn * OMEGA * rhoIn) / (24.0f + OMEGA);

			*mxx = 0.0f;
			*mxy = (36.0f * mxyIn * rhoIn - (*rhoVar)) / (9.0f * (*rhoVar));
			*myy = 0.0f;
		}
		else {
			*rhoVar = -12.0f * rhoIn * (-3.0f - 3.0f * mxxIn + 7.0f * mxyIn - 3.0f * myyIn + 3.0f * mxxIn * OMEGA - 7.0f * mxyIn * OMEGA + 3 * myyIn * OMEGA) / (16.0f + 9.0f * OMEGA);

			*mxx = 2.0f * (9.0f * rhoIn * mxxIn - 6.0f * rhoIn * mxyIn + (*rhoVar)) / (9.0f * (*rhoVar));
			*mxy = -(18.0f * rhoIn * mxxIn - 132.0f * rhoIn * mxyIn + 18.0f * rhoIn * myyIn + 7.0f * (*rhoVar)) / (27.0f * (*rhoVar));
			*myy = -2.0f * (6.0f * rhoIn * mxyIn - 9.0f * rhoIn * myyIn - (*rhoVar)) / (9.0f * (*rhoVar));
		}

		break;
	}
	case SOUTH_EAST: {
		const dfloat rhoIn = pop[0] + pop[1] + pop[4] + pop[8];
		const dfloat inv_rhoIn = 1.0f / rhoIn;

		const dfloat mxxIn = (pop[1] + pop[8]) * inv_rhoIn - cs2;
		const dfloat mxyIn = -pop[8] * inv_rhoIn;
		const dfloat myyIn = (pop[4] + pop[8]) * inv_rhoIn - cs2;


		*ux = 0.0f;
		*uy = 0.0f;

		if (IRBC) {
			*rhoVar = -36.0f * (mxyIn * OMEGA * rhoIn - rhoIn - mxyIn * rhoIn) / (24 + OMEGA);

			*mxx = 0.0f;
			*mxy = (36.0f * mxyIn * rhoIn + (*rhoVar)) / (9.0f * (*rhoVar));
			*myy = 0.0f;
		}
		else {
			*rhoVar = -12.0f * rhoIn * (-3.0f - 3.0f * mxxIn - 7.0f * mxyIn - 3.0f * myyIn + 3.0f * mxxIn * OMEGA + 7.0f * mxyIn * OMEGA + 3 * myyIn * OMEGA) / (16.0f + 9.0f * OMEGA);

			*mxx = 2.0f * (9.0f * rhoIn * mxxIn + 6.0f * rhoIn * mxyIn + (*rhoVar)) / (9.0f * (*rhoVar));
			*mxy = -(-18.0f * rhoIn * mxxIn - 132.0f * rhoIn * mxyIn - 18.0f * rhoIn * myyIn - 7.0f * (*rhoVar)) / (27.0f * (*rhoVar));
			*myy = 2.0f * (6.0f * rhoIn * mxyIn + 9.0f * rhoIn * myyIn + (*rhoVar)) / (9.0f * (*rhoVar));
		}

		break;
	}
	case NORTH_WEST: {
		const dfloat rhoIn = pop[0] + pop[2] + pop[3] + pop[6];
		const dfloat inv_rhoIn = 1.0f / rhoIn;

		const dfloat mxxIn = (pop[3] + pop[6]) * inv_rhoIn - cs2;
		const dfloat mxyIn = -pop[6] * inv_rhoIn;
		const dfloat myyIn = (pop[2] + pop[6]) * inv_rhoIn - cs2;

		*ux = U_MAX;
		*uy = 0.0f;

		if (IRBC) {
			*rhoVar = -36.0f * (mxyIn * OMEGA * rhoIn - rhoIn - mxyIn * rhoIn) / (24.0f + OMEGA + 18.0f * U_MAX - 3.0f * OMEGA * U_MAX - 18.0f * U_MAX * U_MAX + 3.0f * OMEGA * U_MAX * U_MAX);

			*mxx = U_MAX * U_MAX;
			*mxy = (36.0f * mxyIn * rhoIn + (*rhoVar) - 3.0f * U_MAX * (*rhoVar) + 3.0f * U_MAX * U_MAX * (*rhoVar)) / (9.0f * (*rhoVar));
			*myy = 0.0f;
		}
		else {
			*rhoVar = 12.0f * rhoIn * (-3.0f - 3.0f * mxxIn - 7.0f * mxyIn - 3.0f * myyIn + 3.0f * mxxIn * OMEGA + 7.0f * mxyIn * OMEGA + 3 * myyIn * OMEGA) / (-16.0f - 9.0f * OMEGA - 14.0f * U_MAX - OMEGA * U_MAX + 15.0f * OMEGA * U_MAX * U_MAX);

			*mxx = 2.0f * (9.0f * rhoIn * mxxIn + 6.0f * rhoIn * mxyIn + (*rhoVar) + 2.0f * U_MAX * (*rhoVar)) / (9.0f * (*rhoVar));
			*mxy = -(-18.0f * rhoIn * mxxIn - 132.0f * rhoIn * mxyIn - 18.0f * rhoIn * myyIn - 7.0f * (*rhoVar) + 7.0f * U_MAX * (*rhoVar)) / (27.0f * (*rhoVar));
			*myy = -2.0f * (-6.0f * rhoIn * mxyIn - 9.0f * rhoIn * myyIn - (*rhoVar) + U_MAX * (*rhoVar)) / (9.0f * (*rhoVar));
		}

		break;
	}
	case NORTH_EAST: {
		const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[5];
		const dfloat inv_rhoIn = 1.0f / rhoIn;

		const dfloat mxxIn = (pop[1] + pop[5]) * inv_rhoIn - cs2;
		const dfloat mxyIn = pop[5] * inv_rhoIn;
		const dfloat myyIn = (pop[2] + pop[5]) * inv_rhoIn - cs2;

		*ux = U_MAX;
		*uy = 0.0f;

		if (IRBC) {
			*rhoVar = 36.0f * (mxyIn * OMEGA * rhoIn + rhoIn - mxyIn * rhoIn) / (24.0f + OMEGA - 18.0f * U_MAX + 3.0f * OMEGA * U_MAX - 18.0f * U_MAX * U_MAX + 3.0f * OMEGA * U_MAX * U_MAX);

			*mxx = U_MAX * U_MAX;
			*mxy = (36.0f * mxyIn * rhoIn - (*rhoVar) - 3.0f * U_MAX * (*rhoVar) - 3.0f * U_MAX * U_MAX * (*rhoVar)) / (9.0f * (*rhoVar));
			*myy = 0.0f;
		}
		else {
			*rhoVar = 12.0f * rhoIn * (-3.0f - 3.0f * mxxIn + 7.0f * mxyIn - 3.0f * myyIn + 3.0f * mxxIn * OMEGA - 7.0f * mxyIn * OMEGA + 3 * myyIn * OMEGA) / (-16.0f - 9.0f * OMEGA + 14.0f * U_MAX + OMEGA * U_MAX + 15.0f * OMEGA * U_MAX * U_MAX);

			*mxx = -2.0f * (-9.0f * rhoIn * mxxIn + 6.0f * rhoIn * mxyIn - (*rhoVar) + 2.0f * U_MAX * (*rhoVar)) / (9.0f * (*rhoVar));
			*mxy = -(18.0f * rhoIn * mxxIn - 132.0f * rhoIn * mxyIn + 18.0f * rhoIn * myyIn + 7.0f * (*rhoVar) + 7.0f * U_MAX * (*rhoVar)) / (27.0f * (*rhoVar));
			*myy = 2.0f * (-6.0f * rhoIn * mxyIn + 9.0f * rhoIn * myyIn + (*rhoVar) + U_MAX * (*rhoVar)) / (9.0f * (*rhoVar));
		}

		break;
	}
	default:
		break;
	}
}

#endif // BOUNDARY_FUNCTIONS_CUH
