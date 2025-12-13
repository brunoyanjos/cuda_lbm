#include "mlbm.cuh"
#include "globalStructs.h"
#include "globalFunctions.h"
#include "init/state.cuh"
#include "boundaries/node_type.h"
#include "boundaries/boundary_formulation.cuh"
#include "colrec/collision_and_reconstruction.cuh"
#include "interface/interface_handling.cuh"
#include "interpolation_utilities.cuh"

__global__ void streamingAndMom(LBMState state, dfloat OMEGA, ghostInterfaceData ghostInterface)
{
	const int x = threadIdx.x + blockDim.x * blockIdx.x;
	const int y = threadIdx.y + blockDim.y * blockIdx.y;

	if (x >= NX || y >= NY)
		return;

	dfloat pop[Q];

	__shared__ dfloat s_pop[BLOCK_LBM_SIZE * (Q - 1)];

	unsigned int nodeType = state.d_node_type[idxBlock()];

	if (nodeType == SOLID_NODE)
		return;

	dfloat rho = state.d_rho[idxBlock()] + RHO_0;
	dfloat ux = state.d_ux[idxBlock()];
	dfloat uy = state.d_uy[idxBlock()];
	dfloat mxx = state.d_mxx[idxBlock()];
	dfloat mxy = state.d_mxy[idxBlock()];
	dfloat myy = state.d_myy[idxBlock()];

	second::pop_reconstruction(pop, rho, ux, uy, mxx, mxy, myy);

	const unsigned short int xp1 = (threadIdx.x + 1 + BLOCK_NX) % BLOCK_NX;
	const unsigned short int xm1 = (threadIdx.x - 1 + BLOCK_NX) % BLOCK_NX;

	const unsigned short int yp1 = (threadIdx.y + 1 + BLOCK_NY) % BLOCK_NY;
	const unsigned short int ym1 = (threadIdx.y - 1 + BLOCK_NY) % BLOCK_NY;

	// save populations in shared memory
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 0)] = pop[1];
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 1)] = pop[2];
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 2)] = pop[3];
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 3)] = pop[4];
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 4)] = pop[5];
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 5)] = pop[6];
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 6)] = pop[7];
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 7)] = pop[8];

	// sync threads of the block so all populations are saved
	__syncthreads();

	pop[1] = s_pop[idxPopBlock(xm1, threadIdx.y, 0)];
	pop[2] = s_pop[idxPopBlock(threadIdx.x, ym1, 1)];
	pop[3] = s_pop[idxPopBlock(xp1, threadIdx.y, 2)];
	pop[4] = s_pop[idxPopBlock(threadIdx.x, yp1, 3)];
	pop[5] = s_pop[idxPopBlock(xm1, ym1, 4)];
	pop[6] = s_pop[idxPopBlock(xp1, ym1, 5)];
	pop[7] = s_pop[idxPopBlock(xp1, yp1, 6)];
	pop[8] = s_pop[idxPopBlock(xm1, yp1, 7)];

	/* load pop from global in cover nodes */

	pop_load(ghostInterface, threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y, pop);

	dfloat invRho;

	if (nodeType != BULK)
	{
		boundary::eval_incoming_properties(nodeType, pop, rho, ux, uy, mxx, mxy, myy);
	}
	else
	{
		rho = pop[0] + pop[1] + pop[2] + pop[3] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8];
		invRho = static_cast<dfloat>(1) / rho;

		ux = ((pop[1] + pop[5] + pop[8]) - (pop[3] + pop[6] + pop[7])) * invRho;
		uy = ((pop[2] + pop[5] + pop[6]) - (pop[4] + pop[7] + pop[8])) * invRho;

		mxx = (pop[1] + pop[3] + pop[5] + pop[6] + pop[7] + pop[8]) * invRho - cs2;
		mxy = ((pop[5] + pop[7]) - (pop[6] + pop[8])) * invRho;
		myy = (pop[2] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8]) * invRho - cs2;
	}

	state.d_rho[idxBlock()] = rho - RHO_0;

	state.d_ux[idxBlock()] = ux;
	state.d_uy[idxBlock()] = uy;

	state.d_mxx[idxBlock()] = mxx;
	state.d_mxy[idxBlock()] = mxy;
	state.d_myy[idxBlock()] = myy;
}

__global__ void updateBoundaries(LBMState state, dfloat OMEGA, dfloat D_in, dfloat D_out)
{
	const int x = threadIdx.x + blockDim.x * blockIdx.x;
	const int y = threadIdx.y + blockDim.y * blockIdx.y;

	if (x >= NX || y >= NY)
		return;

	unsigned int node_type = state.d_node_type[idxBlock()];

	if (node_type == SOLID_NODE || node_type == BULK)
		return;

	dfloat rho = state.d_rho[idxBlock()] + RHO_0;

	dfloat ux = state.d_ux[idxBlock()];
	dfloat uy = state.d_uy[idxBlock()];

	dfloat mxx = state.d_mxx[idxBlock()];
	dfloat mxy = state.d_mxy[idxBlock()];
	dfloat myy = state.d_myy[idxBlock()];

	const dfloat xb_diff = static_cast<dfloat>(x) - xc;
	const dfloat yb_diff = static_cast<dfloat>(y) - yc;

	const dfloat rb2 = xb_diff * xb_diff + yb_diff * yb_diff;

	const dfloat rb = dsqrt(rb2);
	const dfloat inv_rb = static_cast<dfloat>(1) / rb;

	const dfloat del_x = dsqrt(static_cast<dfloat>(2));

	const dfloat unit_nx = xb_diff * inv_rb;
	const dfloat unit_ny = yb_diff * inv_rb;

	dfloat ux_boundary;
	dfloat uy_boundary;

	dfloat mxx_prime;
	dfloat myy_prime;

	if (node_type & INNER_BOUNDARY)
	{
		const dfloat radii = D_in * static_cast<dfloat>(0.5);
		const dfloat inv_radii = dfloat(1) / radii;

		const dfloat xw = xc + radii * unit_nx;
		const dfloat yw = yc + radii * unit_ny;

		const dfloat xwb_diff = xw - static_cast<dfloat>(x);
		const dfloat ywb_diff = yw - static_cast<dfloat>(y);

		const dfloat dr2 = xwb_diff * xwb_diff + ywb_diff * ywb_diff;
		const dfloat dr = dsqrt(dr2);

		const dfloat x1 = xw + del_x * unit_nx;
		const dfloat y1 = yw + del_x * unit_ny;

		const int int_x1 = int(x1);
		const int int_y1 = int(y1);

		const dfloat x2 = xw + static_cast<dfloat>(2) * del_x * unit_nx;
		const dfloat y2 = yw + static_cast<dfloat>(2) * del_x * unit_ny;

		const int int_x2 = int(x2);
		const int int_y2 = int(y2);

		dfloat ux1 = bilinear_velocity_interpolation(x1, y1, int_x1, int_y1, state.d_ux);
		dfloat ux2 = bilinear_velocity_interpolation(x2, y2, int_x2, int_y2, state.d_ux);

		dfloat uy1 = bilinear_velocity_interpolation(x1, y1, int_x1, int_y1, state.d_uy);
		dfloat uy2 = bilinear_velocity_interpolation(x2, y2, int_x2, int_y2, state.d_uy);

		dfloat mxx1, mxx2;
		dfloat myy1, myy2;

		bilinear_moment_interpolation(x1, y1, int_x1, int_y1, state, mxx1, myy1);
		bilinear_moment_interpolation(x2, y2, int_x2, int_y2, state, mxx2, myy2);

		const dfloat dx = dsqrt(static_cast<dfloat>(2.0));
		const dfloat dx2 = dx * dx;

		const dfloat inv_dr2 = dfloat(1) / dx2;
		const dfloat inv_2dr2 = dfloat(0.5) / dx2;

		const dfloat ux_wall = -U_MAX * (yw - yc) * inv_radii;
		const dfloat uy_wall = U_MAX * (xw - xc) * inv_radii;

		const dfloat mxx_wall = ux_wall * ux_wall;
		const dfloat myy_wall = uy_wall * uy_wall;

		const dfloat wall_term = (dfloat(2) * dx2 - dr2 + dfloat(3) * dr * dx) * inv_2dr2;
		const dfloat one_term = dr * (dr - dfloat(2) * dx) * inv_dr2;
		const dfloat two_term = dr * (dr - dx) * inv_2dr2;

		ux_boundary = wall_term * ux_wall + one_term * ux1 - two_term * ux2;
		uy_boundary = wall_term * uy_wall + one_term * uy1 - two_term * uy2;

		mxx_prime = wall_term * mxx_wall + one_term * mxx1 - two_term * mxx2;
		myy_prime = wall_term * myy_wall + one_term * myy1 - two_term * myy2;
	}
	else
	{
		const dfloat radii = D_out * static_cast<dfloat>(0.5);

		const dfloat xw = xc + radii * unit_nx;
		const dfloat yw = yc + radii * unit_ny;

		const dfloat xwb_diff = xw - static_cast<dfloat>(x);
		const dfloat ywb_diff = yw - static_cast<dfloat>(y);

		const dfloat dr2 = xwb_diff * xwb_diff + ywb_diff * ywb_diff;
		const dfloat dr = dsqrt(dr2);

		const dfloat x1 = xw - del_x * unit_nx;
		const dfloat y1 = yw - del_x * unit_ny;

		const int int_x1 = int(x1);
		const int int_y1 = int(y1);

		const dfloat x2 = xw - static_cast<dfloat>(2) * del_x * unit_nx;
		const dfloat y2 = yw - static_cast<dfloat>(2) * del_x * unit_ny;

		const int int_x2 = int(x2);
		const int int_y2 = int(y2);

		dfloat ux1 = bilinear_velocity_interpolation(x1, y1, int_x1, int_y1, state.d_ux);
		dfloat ux2 = bilinear_velocity_interpolation(x2, y2, int_x2, int_y2, state.d_ux);

		dfloat uy1 = bilinear_velocity_interpolation(x1, y1, int_x1, int_y1, state.d_uy);
		dfloat uy2 = bilinear_velocity_interpolation(x2, y2, int_x2, int_y2, state.d_uy);

		dfloat mxx1, mxx2;
		dfloat myy1, myy2;

		bilinear_moment_interpolation(x1, y1, int_x1, int_y1, state, mxx1, myy1);
		bilinear_moment_interpolation(x2, y2, int_x2, int_y2, state, mxx2, myy2);

		const dfloat dx = dsqrt(static_cast<dfloat>(2.0));
		const dfloat dx2 = dx * dx;

		const dfloat inv_dr2 = dfloat(1) / dx2;
		const dfloat inv_2dr2 = dfloat(0.5) / dx2;

		const dfloat one_term = dr * (dr - dfloat(2) * dx) * inv_dr2;
		const dfloat two_term = dr * (dr - dx) * inv_2dr2;

		ux_boundary = one_term * ux1 - two_term * ux2;
		uy_boundary = one_term * uy1 - two_term * uy2;

		mxx_prime = one_term * mxx1 - two_term * mxx2;
		myy_prime = one_term * myy1 - two_term * myy2;
	}

	boundary::numerical_solution(node_type, OMEGA,
								 mxx_prime, myy_prime,
								 rho, ux_boundary, uy_boundary,
								 mxx, mxy, myy);

	state.d_rho[idxBlock()] = rho - RHO_0;

	state.d_ux[idxBlock()] = ux_boundary;
	state.d_uy[idxBlock()] = uy_boundary;

	state.d_mxx[idxBlock()] = mxx;
	state.d_mxy[idxBlock()] = mxy;
	state.d_myy[idxBlock()] = myy;
}

__global__ void boundaryAndCollision(LBMState state, dfloat OMEGA, ghostInterfaceData ghostInterface)
{
	const int x = threadIdx.x + blockDim.x * blockIdx.x;
	const int y = threadIdx.y + blockDim.y * blockIdx.y;

	if (x >= NX || y >= NY)
		return;

	dfloat pop[Q];

	unsigned int nodeType = state.d_node_type[idxBlock()];

	if (nodeType == SOLID_NODE)
		return;

	dfloat rho = state.d_rho[idxBlock()] + RHO_0;
	dfloat ux = state.d_ux[idxBlock()];
	dfloat uy = state.d_uy[idxBlock()];
	dfloat mxx = state.d_mxx[idxBlock()];
	dfloat mxy = state.d_mxy[idxBlock()];
	dfloat myy = state.d_myy[idxBlock()];

	ux *= F_M_I_SCALE;
	uy *= F_M_I_SCALE;

	mxx *= F_M_II_SCALE;
	mxy *= F_M_IJ_SCALE;
	myy *= F_M_II_SCALE;

	moment_collision(ux, uy, mxx, mxy, myy, OMEGA);

	second::pop_reconstruction(pop, rho, ux, uy, mxx, mxy, myy);

	state.d_rho[idxBlock()] = rho - RHO_0;

	state.d_ux[idxBlock()] = ux;
	state.d_uy[idxBlock()] = uy;

	state.d_mxx[idxBlock()] = mxx;
	state.d_mxy[idxBlock()] = mxy;
	state.d_myy[idxBlock()] = myy;

	pop_save(ghostInterface, threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y, x, y, pop);
}

__global__ void mlbmKernel(LBMState state, dfloat OMEGA, ghostInterfaceData ghostInterface)
{
	const int x = threadIdx.x + blockDim.x * blockIdx.x;
	const int y = threadIdx.y + blockDim.y * blockIdx.y;

	if (x >= NX || y >= NY)
		return;

	dfloat pop[Q];

	__shared__ dfloat s_pop[BLOCK_LBM_SIZE * (Q - 1)];

	unsigned int nodeType = state.d_node_type[idxBlock()];

	dfloat rho = state.d_rho[idxBlock()] + RHO_0;
	dfloat ux = state.d_ux[idxBlock()];
	dfloat uy = state.d_uy[idxBlock()];
	dfloat mxx = state.d_mxx[idxBlock()];
	dfloat mxy = state.d_mxy[idxBlock()];
	dfloat myy = state.d_myy[idxBlock()];

	second::pop_reconstruction(pop, rho, ux, uy, mxx, mxy, myy);

	const unsigned short int xp1 = (threadIdx.x + 1 + BLOCK_NX) % BLOCK_NX;
	const unsigned short int xm1 = (threadIdx.x - 1 + BLOCK_NX) % BLOCK_NX;

	const unsigned short int yp1 = (threadIdx.y + 1 + BLOCK_NY) % BLOCK_NY;
	const unsigned short int ym1 = (threadIdx.y - 1 + BLOCK_NY) % BLOCK_NY;

	// save populations in shared memory
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 0)] = pop[1];
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 1)] = pop[2];
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 2)] = pop[3];
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 3)] = pop[4];
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 4)] = pop[5];
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 5)] = pop[6];
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 6)] = pop[7];
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 7)] = pop[8];

	// sync threads of the block so all populations are saved
	__syncthreads();

	pop[1] = s_pop[idxPopBlock(xm1, threadIdx.y, 0)];
	pop[2] = s_pop[idxPopBlock(threadIdx.x, ym1, 1)];
	pop[3] = s_pop[idxPopBlock(xp1, threadIdx.y, 2)];
	pop[4] = s_pop[idxPopBlock(threadIdx.x, yp1, 3)];
	pop[5] = s_pop[idxPopBlock(xm1, ym1, 4)];
	pop[6] = s_pop[idxPopBlock(xp1, ym1, 5)];
	pop[7] = s_pop[idxPopBlock(xp1, yp1, 6)];
	pop[8] = s_pop[idxPopBlock(xm1, yp1, 7)];

	/* load pop from global in cover nodes */

	pop_load(ghostInterface, threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y, pop);

	dfloat invRho;

	if (nodeType != BULK)
	{
		// if (y == NY - 1)
		// {
		// 	ux = U_MAX;
		// }

		boundary::boundary_calculation(nodeType, rho, ux, uy, mxx, mxy, myy, pop, OMEGA);
	}
	else
	{
		rho = pop[0] + pop[1] + pop[2] + pop[3] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8];
		invRho = static_cast<dfloat>(1) / rho;

		ux = ((pop[1] + pop[5] + pop[8]) - (pop[3] + pop[6] + pop[7])) * invRho;
		uy = ((pop[2] + pop[5] + pop[6]) - (pop[4] + pop[7] + pop[8])) * invRho;

		mxx = (pop[1] + pop[3] + pop[5] + pop[6] + pop[7] + pop[8]) * invRho - cs2;
		mxy = ((pop[5] + pop[7]) - (pop[6] + pop[8])) * invRho;
		myy = (pop[2] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8]) * invRho - cs2;
	}

	ux *= F_M_I_SCALE;
	uy *= F_M_I_SCALE;

	mxx *= F_M_II_SCALE;
	mxy *= F_M_IJ_SCALE;
	myy *= F_M_II_SCALE;

	moment_collision(ux, uy, mxx, mxy, myy, OMEGA);

	second::pop_reconstruction(pop, rho, ux, uy, mxx, mxy, myy);

	state.d_rho[idxBlock()] = rho - RHO_0;

	state.d_ux[idxBlock()] = U_MAX;
	state.d_uy[idxBlock()] = uy;

	state.d_mxx[idxBlock()] = mxx;
	state.d_mxy[idxBlock()] = mxy;
	state.d_myy[idxBlock()] = myy;

	pop_save(ghostInterface, threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y, x, y, pop);
}