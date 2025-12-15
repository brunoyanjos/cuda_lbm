#ifndef BOUNDARIES_CUH
#define BOUNDARIES_CUH

// CUDA INCLUDE
#include <cuda.h>
#include <cuda_runtime.h>
#include "device_launch_parameters.h"

#include "../var.h"
#include "../globalFunctions.h"
#include "../nodeTypeMap.h"

[[nodiscard]] __host__ __device__ inline int boundary_definition(unsigned int x, unsigned int y)
{
	if (x == 0 && y == 0)
	{
		return SOUTH_WEST; // SOUTH_WEST;
	}
	else if (x == 0 && y == (NY - 1))
	{
		return NORTH_WEST; // NORTH_WEST;
	}
	else if (x == (NX - 1) && y == 0)
	{
		return SOUTH_EAST; // SOUTH_EAST;
	}
	else if (x == (NX - 1) && y == (NY - 1))
	{
		return NORTH_EAST; // NORTH_EAST;
	}
	else if (y == 0)
	{
		return SOUTH;
	}
	else if (y == (NY - 1))
	{
		return NORTH;
	}
	else if (x == 0)
	{
		return WEST;
	}
	else if (x == (NX - 1))
	{
		return EAST;
	}

	return BULK;
}

__device__ inline void boundary_calculation_old(unsigned int nodeType, dfloat &rhoVar,
												dfloat &ux, dfloat &uy,
												dfloat &mxx, dfloat &myy, dfloat &mxy,
												dfloat *pop)
{
	switch (nodeType)
	{
	case NORTH:
	{
		const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[3] + pop[5] + pop[6];
		const dfloat inv_rhoIn = 1.0f / rhoIn;

		const dfloat mxyIn = (pop[5] - pop[6]) * inv_rhoIn;

		ux = U_MAX;
		uy = 0.0f;

		rhoVar = 6.0f * rhoIn / 5.0f;

		mxx = U_MAX * U_MAX;
		mxy = 5.0f * mxyIn / 3.0f - U_MAX / 3.0f;
		myy = 0.0f;

		break;
	}
	case SOUTH:
	{
		const dfloat rhoIn = pop[0] + pop[1] + pop[3] + pop[4] + pop[7] + pop[8];
		const dfloat inv_rhoIn = 1.0f / rhoIn;

		const dfloat mxyIn = (pop[7] - pop[8]) * inv_rhoIn;

		ux = 0.0f;
		uy = 0.0f;

		rhoVar = 6.0f * rhoIn / 5.0f;

		mxx = 0.0f;
		mxy = 5.0f * mxyIn / 3.0f;
		myy = 0.0f;

		break;
	}
	case WEST:
	{
		const dfloat rhoIn = pop[0] + pop[2] + pop[3] + pop[4] + pop[6] + pop[7];
		const dfloat inv_rhoIn = 1.0f / rhoIn;

		const dfloat mxyIn = (pop[7] - pop[6]) * inv_rhoIn;

		ux = 0.0f;
		uy = 0.0f;

		rhoVar = 6.0f * rhoIn / 5.0f;

		mxx = 0.0f;
		mxy = 5.0f * mxyIn / 3.0f;
		myy = 0.0f;

		break;
	}
	case EAST:
	{
		const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[4] + pop[5] + pop[8];
		const dfloat inv_rhoIn = 1.0f / rhoIn;

		const dfloat mxyIn = (pop[5] - pop[8]) * inv_rhoIn;

		ux = 0.0f;
		uy = 0.0f;

		rhoVar = 6.0f * rhoIn / 5.0f;

		mxx = 0.0f;
		mxy = 5.0f * mxyIn / 3.0f;
		myy = 0.0f;

		break;
	}
	case SOUTH_WEST:
	{
		const dfloat rhoIn = pop[0] + pop[3] + pop[4] + pop[7];
		const dfloat inv_rhoIn = 1.0f / rhoIn;

		const dfloat mxyIn = pop[7] * inv_rhoIn;

		ux = 0.0f;
		uy = 0.0f;

		rhoVar = 36.0f * (rhoIn - mxyIn * rhoIn + mxyIn * OMEGA * rhoIn) / (24.0f + OMEGA);

		mxx = 0.0f;
		mxy = (36.0f * mxyIn * rhoIn - rhoVar) / (9.0f * rhoVar);
		myy = 0.0f;

		break;
	}
	case SOUTH_EAST:
	{
		const dfloat rhoIn = pop[0] + pop[1] + pop[4] + pop[8];
		const dfloat inv_rhoIn = 1.0f / rhoIn;

		const dfloat mxyIn = -pop[8] * inv_rhoIn;

		ux = 0.0f;
		uy = 0.0f;

		rhoVar = -36.0f * (mxyIn * OMEGA * rhoIn - rhoIn - mxyIn * rhoIn) / (24 + OMEGA);

		mxx = 0.0f;
		mxy = (36.0f * mxyIn * rhoIn + rhoVar) / (9.0f * rhoVar);
		myy = 0.0f;

		break;
	}
	case NORTH_WEST:
	{
		const dfloat rhoIn = pop[0] + pop[2] + pop[3] + pop[6];
		const dfloat inv_rhoIn = 1.0f / rhoIn;

		const dfloat mxyIn = -pop[6] * inv_rhoIn;

		ux = U_MAX;
		uy = 0.0f;

		rhoVar = -36.0f * (mxyIn * OMEGA * rhoIn - rhoIn - mxyIn * rhoIn) /
				 (24.0f + OMEGA + 18.0f * U_MAX - 3.0f * OMEGA * U_MAX - 18.0f * U_MAX * U_MAX + 3.0f * OMEGA * U_MAX * U_MAX);

		mxx = U_MAX * U_MAX;
		mxy = (36.0f * mxyIn * rhoIn + rhoVar - 3.0f * U_MAX * rhoVar + 3.0f * U_MAX * U_MAX * rhoVar) /
			  (9.0f * rhoVar);
		myy = 0.0f;

		break;
	}
	case NORTH_EAST:
	{
		const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[5];
		const dfloat inv_rhoIn = 1.0f / rhoIn;

		const dfloat mxyIn = pop[5] * inv_rhoIn;

		ux = U_MAX;
		uy = 0.0f;

		rhoVar = 36.0f * (mxyIn * OMEGA * rhoIn + rhoIn - mxyIn * rhoIn) /
				 (24.0f + OMEGA - 18.0f * U_MAX + 3.0f * OMEGA * U_MAX - 18.0f * U_MAX * U_MAX + 3.0f * OMEGA * U_MAX * U_MAX);

		mxx = U_MAX * U_MAX;
		mxy = (36.0f * mxyIn * rhoIn - rhoVar - 3.0f * U_MAX * rhoVar - 3.0f * U_MAX * U_MAX * rhoVar) /
			  (9.0f * rhoVar);
		myy = 0.0f;

		break;
	}
	default:
		break;
	}
}

__device__ inline void evaluate_dir(uint8_t node_type,
									uint8_t &incoming_mask,
									uint8_t &outgoing_mask)
{
	const unsigned int out_dirs[4][3] = {
		{3, 4, 7}, // bit_1_dir
		{1, 4, 8}, // bit_2_dir
		{2, 3, 6}, // bit_4_dir
		{1, 2, 5}, // bit_8_dir
	};

	const unsigned int in_dirs[4][3] = {
		{1, 2, 5}, // bit_1_dir
		{2, 3, 6}, // bit_2_dir
		{1, 4, 8}, // bit_4_dir
		{3, 4, 7}, // bit_8_dir
	};

#pragma unroll 4
	for (unsigned int j = 1; j <= 8; j <<= 1)
	{
		if (node_type & j)
		{
			unsigned int idx = __ffs(j) - 1;

#pragma unroll 3
			for (unsigned int k = 0; k < 3; ++k)
			{
				unsigned int od = out_dirs[idx][k] - 1;
				unsigned int id = in_dirs[idx][k] - 1;

				// liga o bit correspondente à direção
				outgoing_mask |= (1u << od);
				incoming_mask |= (1u << id);
			}
		}
	}
}

__device__ inline void boundary_calculation(unsigned int nodeType, dfloat &rho,
											dfloat &ux, dfloat &uy,
											dfloat &mxx, dfloat &myy, dfloat &mxy,
											dfloat *pop)
{
	uint8_t incoming_mask = 0;
	uint8_t outgoing_mask = 0;

	ux = 0.0;
	uy = 0.0;

	if (nodeType == NORTH || nodeType == NORTH_WEST || nodeType == NORTH_EAST)
		ux = U_MAX;

	evaluate_dir(nodeType, incoming_mask, outgoing_mask);

	const dfloat omega_var = static_cast<dfloat>(1) - OMEGA;

	const dfloat Hxx_0 = cx[0] * cx[0] - cs2;
	const dfloat Hxy_0 = cx[0] * cy[0];
	const dfloat Hyy_0 = cy[0] * cy[0] - cs2;

	const dfloat A_0 = w[0] * (1 + as2 * ux * cx[0] + as2 * uy * cy[0]);
	const dfloat Bxx_0 = w[0] * as4 * static_cast<dfloat>(0.5) * Hxx_0;
	const dfloat Bxy_0 = w[0] * as4 * static_cast<dfloat>(0.5) * Hxy_0;
	const dfloat Byy_0 = w[0] * as4 * static_cast<dfloat>(0.5) * Hyy_0;

	dfloat rho_I = pop[0];

	dfloat mxx_I = pop[0] * Hxx_0;
	dfloat mxy_I = pop[0] * Hxy_0;
	dfloat myy_I = pop[0] * Hyy_0;

	dfloat A = A_0;

	dfloat A_Hxx = A_0 * Hxx_0;
	dfloat A_Hxy = A_0 * Hxy_0;
	dfloat A_Hyy = A_0 * Hyy_0;

	dfloat Bxx = Bxx_0;
	dfloat Bxy = Bxy_0;
	dfloat Byy = Byy_0;

	dfloat Bxx_Hxx = Bxx_0 * Hxx_0;
	dfloat Bxy_Hxx = Bxy_0 * Hxx_0;
	dfloat Byy_Hxx = Byy_0 * Hxx_0;

	dfloat Bxx_Hxy = Bxx_0 * Hxy_0;
	dfloat Bxy_Hxy = Bxy_0 * Hxy_0;
	dfloat Byy_Hxy = Byy_0 * Hxy_0;

	dfloat Bxx_Hyy = Bxx_0 * Hyy_0;
	dfloat Bxy_Hyy = Bxy_0 * Hyy_0;
	dfloat Byy_Hyy = Byy_0 * Hyy_0;

#pragma unroll 8
	for (int i = 1; i < Q; ++i)
	{
		const dfloat Hxx = cx[i] * cx[i] - cs2;
		const dfloat Hxy = cx[i] * cy[i];
		const dfloat Hyy = cy[i] * cy[i] - cs2;

		const dfloat A_i = w[i] * (1 + as2 * ux * cx[i] + as2 * uy * cy[i]);
		const dfloat Bxx_i = w[i] * as4 * static_cast<dfloat>(0.5) * Hxx;
		const dfloat Bxy_i = w[i] * as4 * static_cast<dfloat>(0.5) * Hxy;
		const dfloat Byy_i = w[i] * as4 * static_cast<dfloat>(0.5) * Hyy;

		if (outgoing_mask & (1u << (i - 1)))
		{
			A += A_i;

			Bxx += Bxx_i;
			Bxy += Bxy_i;
			Byy += Byy_i;
		}

		if (incoming_mask & (1u << (i - 1)))
		{
			rho_I += pop[i];

			mxx_I += pop[i] * Hxx;
			mxy_I += pop[i] * Hxy;
			myy_I += pop[i] * Hyy;

			Bxx_Hxx += Bxx_i * Hxx;
			Bxy_Hxx += Bxy_i * Hxx;
			Byy_Hxx += Byy_i * Hxx;

			Bxx_Hxy += Bxx_i * Hxy;
			Bxy_Hxy += Bxy_i * Hxy;
			Byy_Hxy += Byy_i * Hxy;

			Bxx_Hyy += Bxx_i * Hyy;
			Bxy_Hyy += Bxy_i * Hyy;
			Byy_Hyy += Byy_i * Hyy;

			A_Hxx += A_i * Hxx;
			A_Hxy += A_i * Hxy;
			A_Hyy += A_i * Hyy;
		}
	}

	const dfloat inv_rhoI = static_cast<dfloat>(1) / rho_I;

	mxx_I *= inv_rhoI;
	mxy_I *= inv_rhoI;
	myy_I *= inv_rhoI;

	const dfloat u_sum = ux * ux * Bxx +
						 static_cast<dfloat>(2) * ux * uy * Bxy +
						 uy * uy * Byy;

	// mxx equations

	const dfloat a11 = omega_var * Bxx * mxx_I - Bxx_Hxx;
	const dfloat a12 = static_cast<dfloat>(2) * (omega_var * Bxy * mxx_I - Bxy_Hxx);
	const dfloat a13 = omega_var * Byy * mxx_I - Byy_Hxx;

	const dfloat b1 = A_Hxx - (A + OMEGA * u_sum) * mxx_I;

	// mxy equations

	const dfloat a21 = omega_var * Bxx * mxy_I - Bxx_Hxy;
	const dfloat a22 = static_cast<dfloat>(2) * (omega_var * Bxy * mxy_I - Bxy_Hxy);
	const dfloat a23 = omega_var * Byy * mxy_I - Byy_Hxy;

	const dfloat b2 = A_Hxy - (A + OMEGA * u_sum) * mxy_I;

	// myy equations

	const dfloat a31 = omega_var * Bxx * myy_I - Bxx_Hyy;
	const dfloat a32 = static_cast<dfloat>(2) * (omega_var * Bxy * myy_I - Bxy_Hyy);
	const dfloat a33 = omega_var * Byy * myy_I - Byy_Hyy;

	const dfloat b3 = A_Hyy - (A + OMEGA * u_sum) * myy_I;

	// solving system

	const dfloat denominator = a13 * a22 * a31 - a12 * a23 * a31 - a13 * a21 * a32 + a11 * a23 * a32 + a12 * a21 * a33 - a11 * a22 * a33;
	const dfloat inv_denominator = static_cast<dfloat>(1) / denominator;

	mxx = (a23 * a32 * b1 - a22 * a33 * b1 - a13 * a32 * b2 + a12 * a33 * b2 + a13 * a22 * b3 - a12 * a23 * b3) * inv_denominator;
	mxy = -(a23 * a31 * b1 - a21 * a33 * b1 - a13 * a31 * b2 + a11 * a33 * b2 + a13 * a21 * b3 - a11 * a23 * b3) * inv_denominator;
	myy = (a22 * a31 * b1 - a21 * a32 * b1 - a12 * a31 * b2 + a11 * a32 * b2 + a12 * a21 * b3 - a11 * a22 * b3) * inv_denominator;

	const dfloat mom_sum = mxx * Bxx + static_cast<dfloat>(2) * mxy * Bxy + myy * Byy;

	const dfloat rho_denominator = A + omega_var * mom_sum + OMEGA * u_sum;
	const dfloat inv_rho = static_cast<dfloat>(1) / rho_denominator;

	rho = rho_I * inv_rho;
}

__device__ inline void boundary_calculation_irbc(unsigned int nodeType, dfloat &rho,
												 dfloat &ux, dfloat &uy,
												 dfloat &mxx, dfloat &myy, dfloat &mxy,
												 dfloat *pop)
{
	uint8_t incoming_mask = 0;
	uint8_t outgoing_mask = 0;

	ux = 0.0;
	uy = 0.0;

	if (nodeType == NORTH || nodeType == NORTH_WEST || nodeType == NORTH_EAST)
		ux = U_MAX;

	evaluate_dir(nodeType, incoming_mask, outgoing_mask);

	const dfloat omega_var = static_cast<dfloat>(1) - OMEGA;

	dfloat rho_I = pop[0];

	dfloat mxy_I = static_cast<dfloat>(0);

	dfloat A = w[0];

	dfloat A_Hxx = -cs2 * w[0];
	dfloat A_Hxy = static_cast<dfloat>(0);
	dfloat A_Hyy = -cs2 * w[0];

	dfloat Bxx = -static_cast<dfloat>(0.5) * as2 * w[0];
	dfloat Bxy = static_cast<dfloat>(0);
	dfloat Byy = -static_cast<dfloat>(0.5) * as2 * w[0];

	dfloat Bxx_Hxx = static_cast<dfloat>(0.5) * w[0];
	dfloat Bxy_Hxx = static_cast<dfloat>(0);
	dfloat Byy_Hxx = static_cast<dfloat>(0.5) * w[0];

	dfloat Bxx_Hxy = static_cast<dfloat>(0);
	dfloat Bxy_Hxy = static_cast<dfloat>(0);
	dfloat Byy_Hxy = static_cast<dfloat>(0);

	dfloat Bxx_Hyy = static_cast<dfloat>(0.5) * w[0];
	dfloat Bxy_Hyy = static_cast<dfloat>(0);
	dfloat Byy_Hyy = static_cast<dfloat>(0.5) * w[0];

#pragma unroll 8
	for (int i = 1; i < Q; ++i)
	{
		const dfloat Hxx = cx[i] * cx[i] - cs2;
		const dfloat Hxy = cx[i] * cy[i];
		const dfloat Hyy = cy[i] * cy[i] - cs2;

		const dfloat A_i = w[i] * (1 + as2 * ux * cx[i] + as2 * uy * cy[i]);
		const dfloat Bxx_i = w[i] * as4 * static_cast<dfloat>(0.5) * Hxx;
		const dfloat Bxy_i = w[i] * as4 * static_cast<dfloat>(0.5) * Hxy;
		const dfloat Byy_i = w[i] * as4 * static_cast<dfloat>(0.5) * Hyy;

		if (outgoing_mask & (1u << (i - 1)))
		{
			A += A_i;

			Bxx += Bxx_i;
			Bxy += Bxy_i;
			Byy += Byy_i;
		}

		if (incoming_mask & (1u << (i - 1)))
		{
			rho_I += pop[i];

			mxy_I += pop[i] * Hxy;

			Bxx_Hxx += Bxx_i * Hxx;
			Bxy_Hxx += Bxy_i * Hxx;
			Byy_Hxx += Byy_i * Hxx;

			Bxx_Hxy += Bxx_i * Hxy;
			Bxy_Hxy += Bxy_i * Hxy;
			Byy_Hxy += Byy_i * Hxy;

			Bxx_Hyy += Bxx_i * Hyy;
			Bxy_Hyy += Bxy_i * Hyy;
			Byy_Hyy += Byy_i * Hyy;

			A_Hxx += A_i * Hxx;
			A_Hxy += A_i * Hxy;
			A_Hyy += A_i * Hyy;
		}
	}

	const dfloat inv_rho_I = static_cast<dfloat>(1) / rho_I;

	mxy_I *= inv_rho_I;

	const dfloat u_sum = ux * ux * Bxx +
						 static_cast<dfloat>(2) * ux * uy * Bxy +
						 uy * uy * Byy;

	mxx = ux * ux;
	myy = uy * uy;

	const dfloat mxy_denominator = static_cast<dfloat>(2) * (omega_var * Bxy * mxy_I - Bxy_Hxy);

	const dfloat mxx_xy = omega_var * Bxx * mxy_I - Bxx_Hxy;
	const dfloat myy_xy = omega_var * Byy * mxy_I - Byy_Hxy;

	const dfloat xy_trace = mxx_xy * mxx + myy_xy * myy;

	const dfloat mxy_nominator = A_Hxy - (A + OMEGA * u_sum) * mxy_I - xy_trace;

	mxy = mxy_nominator / mxy_denominator;

	const dfloat mxx_factor = mxx * Bxx;
	const dfloat mxy_factor = static_cast<dfloat>(2) * mxy * Bxy;
	const dfloat myy_factor = myy * Byy;

	const dfloat mom_sum = mxx_factor + mxy_factor + myy_factor;

	const dfloat rho_denominator = A + omega_var * mom_sum + OMEGA * u_sum;

	rho = rho_I / rho_denominator;
}

#endif // BOUNDARY_FUNCTIONS_CUH
