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

__host__ __device__ inline void boundary_definition(unsigned int *nodeType, unsigned int x, unsigned int y)
{
	if (x == 0 && y == 0)
	{
#ifdef BC_Y_PERIODIC
		*nodeType = WEST;
#else
		*nodeType = SOUTH_WEST;
#endif
	}
	else if (x == 0 && y == (NY - 1))
	{
#ifdef BC_Y_PERIODIC
		*nodeType = WEST;
#else
		*nodeType = NORTH_WEST;
#endif
	}
	else if (x == (NX - 1) && y == 0)
	{
#ifdef BC_Y_PERIODIC
		*nodeType = EAST;
#else
		*nodeType = SOUTH_EAST;
#endif
	}
	else if (x == (NX - 1) && y == (NY - 1))
	{
#ifdef BC_Y_PERIODIC
		*nodeType = EAST;
#else
		*nodeType = NORTH_EAST;

#endif
	}
	else if (y == 0)
	{
#ifdef BC_Y_PERIODIC
		*nodeType = BULK;
#else
		*nodeType = SOUTH;
#endif
	}
	else if (y == (NY - 1))
	{
#ifdef BC_Y_PERIODIC
		*nodeType = BULK;
#else
		*nodeType = NORTH;
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

__device__ inline void boundary_calculation(unsigned int nodeType, dfloat *rhoVar, dfloat *ux, dfloat *uy, dfloat *mxx, dfloat *myy, dfloat *mxy, dfloat *pop, dfloat *fMom, int x, int y, dfloat OMEGA)
{
	switch (nodeType)
	{
	case NORTH:
	{
		const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[3] + pop[5] + pop[6];
		const dfloat inv_rhoIn = 1.0 / rhoIn;

		const dfloat mxxIn = (pop[1] + pop[3] + pop[5] + pop[6]) * inv_rhoIn - cs2;
		const dfloat mxyIn = (pop[5] - pop[6]) * inv_rhoIn;
		const dfloat myyIn = (pop[2] + pop[5] + pop[6]) * inv_rhoIn - cs2;

		*ux = 0.0;
		*uy = 0.0;

		*rhoVar = 3.0 * rhoIn * (4.0 + 3.0 * (1.0 - OMEGA) * myyIn) / (9.0 + OMEGA);

		*mxx = 6.0 * rhoIn * mxxIn / (5.0 * (*rhoVar));
		*mxy = 2.0 * rhoIn * mxyIn / (*rhoVar);
		*myy = (*rhoVar + 9.0 * rhoIn * myyIn) / (6.0 * (*rhoVar));

		break;
	}
	case SOUTH:
	{
		const dfloat rhoIn = pop[0] + pop[1] + pop[3] + pop[4] + pop[7] + pop[8];
		const dfloat inv_rhoIn = 1.0 / rhoIn;

		const dfloat mxxIn = (pop[1] + pop[3] + pop[7] + pop[8]) * inv_rhoIn - cs2;
		const dfloat mxyIn = (pop[7] - pop[8]) * inv_rhoIn;
		const dfloat myyIn = (pop[3] + pop[7] + pop[8]) * inv_rhoIn - cs2;

		*ux = 0.0;
		*uy = 0.0;

		*rhoVar = 3.0 * rhoIn * (4.0 + 3.0 * (1.0 - OMEGA) * myyIn) / (9.0 + OMEGA);

		*mxx = 6.0 * rhoIn * mxxIn / (5.0 * (*rhoVar));
		*mxy = 2.0 * rhoIn * mxyIn / (*rhoVar);
		*myy = (*rhoVar + 9.0 * rhoIn * myyIn) / (6.0 * (*rhoVar));

		break;
	}
	case WEST:
	{
		const dfloat rhoIn = pop[0] + pop[2] + pop[3] + pop[4] + pop[6] + pop[7];
		const dfloat inv_rhoIn = 1.0 / rhoIn;

		const dfloat mxxIn = (pop[3] + pop[6] + pop[7]) * inv_rhoIn - cs2;
		const dfloat mxyIn = (pop[7] - pop[6]) * inv_rhoIn;
		const dfloat myyIn = (pop[2] + pop[4] + pop[6] + pop[7]) * inv_rhoIn - cs2;

		*uy = 0.0;
		*ux = U_MAX;

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

		const dfloat rho = RHO_0 + fMom[idxMom(threadIdx.x - 1, threadIdx.y, M_RHO_INDEX, blockIdx.x, blockIdx.y)];
		*ux = fMom[idxMom(threadIdx.x - 1, threadIdx.y, M_UX_INDEX, blockIdx.x, blockIdx.y)] / F_M_I_SCALE;
		*uy = fMom[idxMom(threadIdx.x - 1, threadIdx.y, M_UY_INDEX, blockIdx.x, blockIdx.y)] / F_M_I_SCALE;

		*rhoVar = rho;

		*mxx = (rho + 9.0 * rhoIn * mxxIn - 3.0 * rho * (*ux)) / (6.0 * rho);
		*mxy = (6.0 * rhoIn * mxyIn - rho * (*uy)) / (3.0 * rho);
		*myy = 6.0 * rhoIn * myyIn / (5.0 * rho);

		break;
	}

	case SOUTH_WEST:
	{
		const dfloat rhoIn = pop[0] + pop[3] + pop[4] + pop[7];
		const dfloat inv_rhoIn = 1.0 / rhoIn;

		const dfloat mxxIn = (pop[3] + pop[7]) * inv_rhoIn - cs2;
		const dfloat mxyIn = pop[7] * inv_rhoIn;
		const dfloat myyIn = (pop[4] + pop[7]) * inv_rhoIn - cs2;

		*ux = 0.0;
		*uy = 0.0;

		if (IRBC)
		{
			*rhoVar = 36.0 * (rhoIn - mxyIn * rhoIn + mxyIn * OMEGA * rhoIn) / (24.0 + OMEGA);

			*mxx = 0.0;
			*mxy = (36.0 * mxyIn * rhoIn - (*rhoVar)) / (9.0 * (*rhoVar));
			*myy = 0.0;
		}
		else
		{
			*rhoVar = -12.0 * rhoIn * (-3.0 - 3.0 * mxxIn + 7.0 * mxyIn - 3.0 * myyIn + 3.0 * mxxIn * OMEGA - 7.0 * mxyIn * OMEGA + 3 * myyIn * OMEGA) / (16.0 + 9.0 * OMEGA);

			*mxx = 2.0 * (9.0 * rhoIn * mxxIn - 6.0 * rhoIn * mxyIn + (*rhoVar)) / (9.0 * (*rhoVar));
			*mxy = -(18.0 * rhoIn * mxxIn - 132.0 * rhoIn * mxyIn + 18.0 * rhoIn * myyIn + 7.0 * (*rhoVar)) / (27.0 * (*rhoVar));
			*myy = -2.0 * (6.0 * rhoIn * mxyIn - 9.0 * rhoIn * myyIn - (*rhoVar)) / (9.0 * (*rhoVar));
		}

		break;
	}
	case SOUTH_EAST:
	{
		const dfloat rhoIn = pop[0] + pop[1] + pop[4] + pop[8];
		const dfloat inv_rhoIn = 1.0 / rhoIn;

		const dfloat mxxIn = (pop[1] + pop[8]) * inv_rhoIn - cs2;
		const dfloat mxyIn = -pop[8] * inv_rhoIn;
		const dfloat myyIn = (pop[4] + pop[8]) * inv_rhoIn - cs2;

		*ux = 0.0;
		*uy = 0.0;

		if (IRBC)
		{
			*rhoVar = -36.0 * (mxyIn * OMEGA * rhoIn - rhoIn - mxyIn * rhoIn) / (24 + OMEGA);

			*mxx = 0.0;
			*mxy = (36.0 * mxyIn * rhoIn + (*rhoVar)) / (9.0 * (*rhoVar));
			*myy = 0.0;
		}
		else
		{
			*rhoVar = -12.0 * rhoIn * (-3.0 - 3.0 * mxxIn - 7.0 * mxyIn - 3.0 * myyIn + 3.0 * mxxIn * OMEGA + 7.0 * mxyIn * OMEGA + 3 * myyIn * OMEGA) / (16.0 + 9.0 * OMEGA);

			*mxx = 2.0 * (9.0 * rhoIn * mxxIn + 6.0 * rhoIn * mxyIn + (*rhoVar)) / (9.0 * (*rhoVar));
			*mxy = -(-18.0 * rhoIn * mxxIn - 132.0 * rhoIn * mxyIn - 18.0 * rhoIn * myyIn - 7.0 * (*rhoVar)) / (27.0 * (*rhoVar));
			*myy = 2.0 * (6.0 * rhoIn * mxyIn + 9.0 * rhoIn * myyIn + (*rhoVar)) / (9.0 * (*rhoVar));
		}

		break;
	}
	case NORTH_WEST:
	{
		const dfloat rhoIn = pop[0] + pop[2] + pop[3] + pop[6];
		const dfloat inv_rhoIn = 1.0 / rhoIn;

		const dfloat mxxIn = (pop[3] + pop[6]) * inv_rhoIn - cs2;
		const dfloat mxyIn = -pop[6] * inv_rhoIn;
		const dfloat myyIn = (pop[2] + pop[6]) * inv_rhoIn - cs2;

		*ux = 0.0;
		*uy = 0.0;

		if (IRBC)
		{
			*rhoVar = -36.0 * (mxyIn * OMEGA * rhoIn - rhoIn - mxyIn * rhoIn) / (24 + OMEGA);

			*mxx = 0.0;
			*mxy = (36.0 * mxyIn * rhoIn + (*rhoVar)) / (9.0 * (*rhoVar));
			*myy = 0.0;
		}
		else
		{
			*rhoVar = -12.0 * rhoIn * (-3.0 - 3.0 * mxxIn - 7.0 * mxyIn - 3.0 * myyIn + 3.0 * mxxIn * OMEGA + 7.0 * mxyIn * OMEGA + 3 * myyIn * OMEGA) / (16.0 + 9.0 * OMEGA);

			*mxx = 2.0 * (9.0 * rhoIn * mxxIn + 6.0 * rhoIn * mxyIn + (*rhoVar)) / (9.0 * (*rhoVar));
			*mxy = -(-18.0 * rhoIn * mxxIn - 132.0 * rhoIn * mxyIn - 18.0 * rhoIn * myyIn - 7.0 * (*rhoVar)) / (27.0 * (*rhoVar));
			*myy = 2.0 * (6.0 * rhoIn * mxyIn + 9.0 * rhoIn * myyIn + (*rhoVar)) / (9.0 * (*rhoVar));
		}

		break;
	}
	case NORTH_EAST:
	{
		const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[5];
		const dfloat inv_rhoIn = 1.0 / rhoIn;

		const dfloat mxxIn = (pop[1] + pop[5]) * inv_rhoIn - cs2;
		const dfloat mxyIn = pop[5] * inv_rhoIn;
		const dfloat myyIn = (pop[2] + pop[5]) * inv_rhoIn - cs2;

		*ux = 0.0;
		*uy = 0.0;

		if (IRBC)
		{
			*rhoVar = 36.0 * (rhoIn - mxyIn * rhoIn + mxyIn * OMEGA * rhoIn) / (24.0 + OMEGA);

			*mxx = 0.0;
			*mxy = (36.0 * mxyIn * rhoIn - (*rhoVar)) / (9.0 * (*rhoVar));
			*myy = 0.0;
		}
		else
		{
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

__device__ inline void immersed_boundary_treatment(
	unsigned int nodeType,
	dfloat *rhoVar, dfloat *ux, dfloat *uy, dfloat *mxx, dfloat *mxy, dfloat *myy,
	unsigned int tx, unsigned int ty, unsigned int bx, unsigned int by, unsigned int x, unsigned int y,
	cylinderProperties *cylinder_properties, size_t cylinder_counter,
	dfloat *fMom_old, dfloat OMEGA, unsigned int step)
{
	cylinderProperties property = findCylindeProperty(cylinder_properties, cylinder_counter, x, y);

	// for first point

	// dfloat ux1;
	// dfloat uy1;

	// bilinear_velocity_interpolation(property.x1, property.y1, int(property.x1), int(property.y1), int(property.x1) + 1, int(property.y1) + 1, fMom_old, M_UX_INDEX, &ux1);
	// bilinear_velocity_interpolation(property.x1, property.y1, int(property.x1), int(property.y1), int(property.x1) + 1, int(property.y1) + 1, fMom_old, M_UY_INDEX, &uy1);

	// for second point

	// dfloat ux2;
	// dfloat uy2;

	// bilinear_velocity_interpolation(property.x2, property.y2, int(property.x2), int(property.y2), int(property.x2) + 1, int(property.y2) + 1, fMom_old, M_UX_INDEX, &ux2);
	// bilinear_velocity_interpolation(property.x2, property.y2, int(property.x2), int(property.y2), int(property.x2) + 1, int(property.y2) + 1, fMom_old, M_UY_INDEX, &uy2);

	// moment interpolation to first point
	// dfloat mxx1 = 0.0;
	// dfloat myy1 = 0.0;
	// dfloat mxx2 = 0.0;
	// dfloat myy2 = 0.0;

	// if (ROTATIONAL_COORDINATES) {
	// 	/*bilinear_moment_interpolation(property.x1, property.y1, property.lx1, property.ly1, property.lx1p1, property.ly1p1, fMom_old, &mxx1, &myy1);
	// 	bilinear_moment_interpolation(property.x2, property.y2, property.lx2, property.ly2, property.lx2p1, property.ly2p1, fMom_old, &mxx2, &myy2);*/

	// 	bilinear_moment_interpolation(property.x1, property.y1, int(property.x1), int(property.y1), int(property.x1) + 1, int(property.y1) + 1, fMom_old, &mxx1, &myy1);
	// 	bilinear_moment_interpolation(property.x2, property.y2, int(property.x2), int(property.y2), int(property.x2) + 1, int(property.y2) + 1, fMom_old, &mxx2, &myy2);
	// }

	// if (CALCULATE_PRESSURE && step >= N_STEPS - PRESSURE_TIME) {
	// 	dfloat rho1;
	// 	dfloat rho2;
	// 	dfloat rho3;

	// 	/*bilinear_velocity_interpolation(property.x1, property.y1, property.lx1, property.ly1, property.lx1p1, property.ly1p1, fMom_old, M_RHO_INDEX, &rho1);
	// 	bilinear_velocity_interpolation(property.x2, property.y2, property.lx2, property.ly2, property.lx2p1, property.ly2p1, fMom_old, M_RHO_INDEX, &rho2);
	// 	bilinear_velocity_interpolation(property.x3, property.y3, property.lx3, property.ly3, property.lx3p1, property.ly3p1, fMom_old, M_RHO_INDEX, &rho3);*/

	// 	bilinear_velocity_interpolation(property.x1, property.y1, int(property.x1), int(property.y1), int(property.x1) + 1, int(property.y1) + 1, fMom_old, M_RHO_INDEX, &rho1);
	// 	bilinear_velocity_interpolation(property.x2, property.y2, int(property.x2), int(property.y2), int(property.x2) + 1, int(property.y2) + 1, fMom_old, M_RHO_INDEX, &rho2);
	// 	bilinear_velocity_interpolation(property.x3, property.y3, int(property.x3), int(property.y3), int(property.x3) + 1, int(property.y3) + 1, fMom_old, M_RHO_INDEX, &rho3);

	// 	pressure_extrapolation(property.xw, property.yw, property.x1, property.y1, property.x2, property.y2, property.x3, property.y3, rho1, rho2, rho3, &(cylinder_properties[index].ps));
	// }

	// dfloat delta = property.dr;

	// *ux = extrapolation(delta, ux1, ux2);
	// *uy = extrapolation(delta, uy1, uy2);

	// dfloat m_xx_int = extrapolation(delta, mxx1, mxx2);
	// dfloat m_yy_int = extrapolation(delta, myy1, myy2);

	*ux = 0.0;
	*uy = 0.0;

	if (ROTATIONAL_COORDINATES)
	{
		// numericalSolution_rotation(rhoVar, *ux, *uy, mxx, mxy, myy, m_xx_int, m_yy_int, incomings, outgoings, OMEGA, x, y);
	}
	else
	{
		numericalSolution(rhoVar, *ux, *uy, mxx, mxy, myy, property.is, property.os, OMEGA, x, y);
	}
}

#endif // BOUNDARY_FUNCTIONS_CUH
