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

[[nodiscard]] __host__ __device__ inline int annular_boundary_definition(unsigned int x, unsigned int y)
{
	const dfloat xc_local = xc;
	const dfloat yc_local = yc;

	// ---------------------------
	// Distance from the cell to the center
	// ---------------------------
	const dfloat dx = dfloat(x) - xc_local;
	const dfloat dy = dfloat(y) - yc_local;
	const dfloat dist = dsqrt(dx * dx + dy * dy);

	// ---------------------------
	// Radii definitions
	// ---------------------------
	const dfloat inner_radius = dfloat(D) / 2.0;
	const dfloat outer_radius = dfloat(NX - 1) / 2.0;

	// Smoothing thickness (LBM-friendly transition band)
	const dfloat smooth = 0.5;

	// ---------------------------
	// INNER BOUNDARY (smooth band)
	// ---------------------------
	const dfloat d_inner = dabs(dist - inner_radius);

	// Inside the inner solid region OR within the smoothing band → solid
	if (dist <= inner_radius || d_inner <= smooth)
	{
		return SOLID_NODE;
	}

	// ---------------------------
	// OUTER BOUNDARY (smooth band)
	// ---------------------------
	const dfloat d_outer = dabs(dist - outer_radius);

	// Outside the outer radius → solid region
	// Within the smoothing band around the outer boundary → solid
	if (d_outer <= smooth || dist >= outer_radius)
	{
		return SOLID_NODE;
	}

	// ---------------------------
	// FLUID REGION (between inner and outer radii)
	// ---------------------------
	return BULK;
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

__device__ inline void evaluate_incomings(uint8_t node_type, dfloat *pop,
										  dfloat &rho, dfloat &ux, dfloat &uy,
										  dfloat &mxx, dfloat &mxy, dfloat &myy)
{
	uint8_t incoming_mask = 0;
	uint8_t outgoing_mask = 0;

	evaluate_dir(node_type, incoming_mask, outgoing_mask);

	dfloat rho_I = pop[0];

	dfloat mxx_I = -cs2 * pop[0];
	dfloat mxy_I = static_cast<dfloat>(0);
	dfloat myy_I = -cs2 * pop[0];

#pragma unroll 8
	for (int i = 1; i < 9; ++i)
	{
		const dfloat Hxx = cx[i] * cx[i] - cs2;
		const dfloat Hxy = cx[i] * cy[i];
		const dfloat Hyy = cy[i] * cy[i] - cs2;

		if (incoming_mask & (1u << (i - 1)))
		{
			rho_I += pop[i];

			mxx_I += pop[i] * Hxx;
			mxy_I += pop[i] * Hxy;
			myy_I += pop[i] * Hyy;
		}
	}

	const dfloat inv_rho_I = static_cast<dfloat>(1) / rho_I;

	rho = rho_I;

	mxx = mxx_I * inv_rho_I;
	mxy = mxy_I * inv_rho_I;
	myy = myy_I * inv_rho_I;
}

__device__ inline void boundary_calculation(unsigned int nodeType, dfloat &rho,
											dfloat &ux, dfloat &uy,
											dfloat &mxx, dfloat &mxy, dfloat &myy,
											dfloat OMEGA)
{
	uint8_t incoming_mask = 0;
	uint8_t outgoing_mask = 0;

	evaluate_dir(nodeType, incoming_mask, outgoing_mask);

	const dfloat omega_var = static_cast<dfloat>(1) - OMEGA;

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

	const dfloat u_sum = ux * ux * Bxx +
						 static_cast<dfloat>(2) * ux * uy * Bxy +
						 uy * uy * Byy;

	// mxx equations

	const dfloat a11 = omega_var * Bxx * mxx - Bxx_Hxx;
	const dfloat a12 = static_cast<dfloat>(2) * (omega_var * Bxy * mxx - Bxy_Hxx);
	const dfloat a13 = omega_var * Byy * mxx - Byy_Hxx;

	const dfloat b1 = A_Hxx - (A + OMEGA * u_sum) * mxx;

	// mxy equations

	const dfloat a21 = omega_var * Bxx * mxy - Bxx_Hxy;
	const dfloat a22 = static_cast<dfloat>(2) * (omega_var * Bxy * mxy - Bxy_Hxy);
	const dfloat a23 = omega_var * Byy * mxy - Byy_Hxy;

	const dfloat b2 = A_Hxy - (A + OMEGA * u_sum) * mxy;

	// myy equations

	const dfloat a31 = omega_var * Bxx * myy - Bxx_Hyy;
	const dfloat a32 = static_cast<dfloat>(2) * (omega_var * Bxy * myy - Bxy_Hyy);
	const dfloat a33 = omega_var * Byy * myy - Byy_Hyy;

	const dfloat b3 = A_Hyy - (A + OMEGA * u_sum) * myy;

	// solving system

	const dfloat denominator = a13 * a22 * a31 - a12 * a23 * a31 - a13 * a21 * a32 + a11 * a23 * a32 + a12 * a21 * a33 - a11 * a22 * a33;
	const dfloat inv_denominator = static_cast<dfloat>(1) / denominator;

	mxx = (a23 * a32 * b1 - a22 * a33 * b1 - a13 * a32 * b2 + a12 * a33 * b2 + a13 * a22 * b3 - a12 * a23 * b3) * inv_denominator;
	mxy = -(a23 * a31 * b1 - a21 * a33 * b1 - a13 * a31 * b2 + a11 * a33 * b2 + a13 * a21 * b3 - a11 * a23 * b3) * inv_denominator;
	myy = (a22 * a31 * b1 - a21 * a32 * b1 - a12 * a31 * b2 + a11 * a32 * b2 + a12 * a21 * b3 - a11 * a22 * b3) * inv_denominator;

	const dfloat mom_sum = mxx * Bxx + static_cast<dfloat>(2) * mxy * Bxy + myy * Byy;

	const dfloat rho_denominator = A + omega_var * mom_sum + OMEGA * u_sum;
	const dfloat inv_rho = static_cast<dfloat>(1) / rho_denominator;

	rho = rho * inv_rho;
}

__device__ inline void boundary_calculation_irbc(unsigned int nodeType, dfloat &rho,
												 dfloat &ux, dfloat &uy,
												 dfloat &mxx, dfloat &myy, dfloat &mxy,
												 dfloat *pop, dfloat OMEGA)
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
