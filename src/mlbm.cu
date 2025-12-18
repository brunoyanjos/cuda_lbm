#include "mlbm.cuh"
#include "globalStructs.h"
#include "globalFunctions.h"
#include "interpolation_utilities.cuh"
#include "includeFiles/interface_handling.cuh"

__global__ void mlbmKernel(LBMState state, ghostInterfaceData ghostInterface)
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

	dfloat rho = RHO_0 + state.d_rho[idxBlock()];
	dfloat ux = state.d_ux[idxBlock()];
	dfloat uy = state.d_uy[idxBlock()];
	dfloat mxx = state.d_mxx[idxBlock()];
	dfloat mxy = state.d_mxy[idxBlock()];
	dfloat myy = state.d_myy[idxBlock()];

	pop_reconstruction(rho, ux, uy, mxx, myy, mxy, pop);

	const unsigned short int xp1 = (threadIdx.x + 1 + BLOCK_NX) % BLOCK_NX;
	const unsigned short int xm1 = (threadIdx.x - 1 + BLOCK_NX) % BLOCK_NX;

	const unsigned short int yp1 = (threadIdx.y + 1 + BLOCK_NY) % BLOCK_NY;
	const unsigned short int ym1 = (threadIdx.y - 1 + BLOCK_NY) % BLOCK_NY;

	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 0)] = pop[1];
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 1)] = pop[2];
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 2)] = pop[3];
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 3)] = pop[4];
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 4)] = pop[5];
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 5)] = pop[6];
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 6)] = pop[7];
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 7)] = pop[8];

	__syncthreads();

	pop[1] = s_pop[idxPopBlock(xm1, threadIdx.y, 0)];
	pop[2] = s_pop[idxPopBlock(threadIdx.x, ym1, 1)];
	pop[3] = s_pop[idxPopBlock(xp1, threadIdx.y, 2)];
	pop[4] = s_pop[idxPopBlock(threadIdx.x, yp1, 3)];
	pop[5] = s_pop[idxPopBlock(xm1, ym1, 4)];
	pop[6] = s_pop[idxPopBlock(xp1, ym1, 5)];
	pop[7] = s_pop[idxPopBlock(xp1, yp1, 6)];
	pop[8] = s_pop[idxPopBlock(xm1, yp1, 7)];

	pop_load(ghostInterface, threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y, pop);

	dfloat invRho;

	if (nodeType != BULK)
	{
		boundary_calculation(nodeType, pop, rho,
							 ux, uy,
							 mxx, mxy, myy,
							 OMEGA);
	}
	else
	{
		rho = pop[0] + pop[1] + pop[2] + pop[3] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8];
		invRho = 1 / rho;

		ux = ((pop[1] + pop[5] + pop[8]) - (pop[3] + pop[6] + pop[7])) * invRho;
		uy = ((pop[2] + pop[5] + pop[6]) - (pop[4] + pop[7] + pop[8])) * invRho;

		mxx = (pop[1] + pop[3] + pop[5] + pop[6] + pop[7] + pop[8]) * invRho - cs2;
		mxy = ((pop[5] + pop[7]) - (pop[6] + pop[8])) * invRho;
		myy = (pop[2] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8]) * invRho - cs2;
	}

	ux = F_M_I_SCALE * ux;
	uy = F_M_I_SCALE * uy;

	mxx = F_M_II_SCALE * mxx;
	mxy = F_M_IJ_SCALE * mxy;
	myy = F_M_II_SCALE * myy;

	moment_collision(ux, uy, mxx, myy, mxy, OMEGA);

	pop_reconstruction(rho, ux, uy, mxx, myy, mxy, pop);

	state.d_rho[idxBlock()] = rho - RHO_0;

	state.d_ux[idxBlock()] = ux;
	state.d_uy[idxBlock()] = uy;

	state.d_mxx[idxBlock()] = mxx;
	state.d_mxy[idxBlock()] = mxy;
	state.d_myy[idxBlock()] = myy;

	pop_save(ghostInterface, threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y, x, y, pop);
}