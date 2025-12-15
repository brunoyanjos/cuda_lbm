#include "mlbm.cuh"
#include "globalStructs.h"
#include "globalFunctions.h"

__global__ void gpuMomCollisionStream(
	LBMState state, unsigned int *dNodeType, ghostInterfaceData ghostInterface, unsigned int step)
{
	const int x = threadIdx.x + blockDim.x * blockIdx.x;
	const int y = threadIdx.y + blockDim.y * blockIdx.y;

	if (x >= NX || y >= NY)
		return;
	dfloat pop[Q];
	__shared__ dfloat s_pop[BLOCK_LBM_SIZE * (Q - 1)];

	// Load moments from global memory

	// rho'
	unsigned int nodeType = dNodeType[idxScalarBlock(threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y)];
	if (nodeType == 0b11111111)
		return;
	dfloat rho = RHO_0 + state.d_rho[idxScalarBlock(threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y)];
	dfloat ux = state.d_ux[idxScalarBlock(threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y)];
	dfloat uy = state.d_uy[idxScalarBlock(threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y)];
	dfloat mxx = state.d_mxx[idxScalarBlock(threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y)];
	dfloat mxy = state.d_mxy[idxScalarBlock(threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y)];
	dfloat myy = state.d_myy[idxScalarBlock(threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y)];

	pop_reconstruction(rho, ux, uy, mxx, myy, mxy, pop);

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

	/* pull */

	pop[1] = s_pop[idxPopBlock(xm1, threadIdx.y, 0)];
	pop[2] = s_pop[idxPopBlock(threadIdx.x, ym1, 1)];
	pop[3] = s_pop[idxPopBlock(xp1, threadIdx.y, 2)];
	pop[4] = s_pop[idxPopBlock(threadIdx.x, yp1, 3)];
	pop[5] = s_pop[idxPopBlock(xm1, ym1, 4)];
	pop[6] = s_pop[idxPopBlock(xp1, ym1, 5)];
	pop[7] = s_pop[idxPopBlock(xp1, yp1, 6)];
	pop[8] = s_pop[idxPopBlock(xm1, yp1, 7)];

	const int tx = threadIdx.x;
	const int ty = threadIdx.y;

	const int bx = blockIdx.x;
	const int by = blockIdx.y;

	const int txm1 = (tx - 1 + BLOCK_NX) % BLOCK_NX;
	const int txp1 = (tx + 1 + BLOCK_NX) % BLOCK_NX;

	const int tym1 = (ty - 1 + BLOCK_NY) % BLOCK_NY;
	const int typ1 = (ty + 1 + BLOCK_NY) % BLOCK_NY;

	const int bxm1 = (bx - 1 + NUM_BLOCK_X) % NUM_BLOCK_X;
	const int bxp1 = (bx + 1 + NUM_BLOCK_X) % NUM_BLOCK_X;

	const int bym1 = (by - 1 + NUM_BLOCK_Y) % NUM_BLOCK_Y;
	const int byp1 = (by + 1 + NUM_BLOCK_Y) % NUM_BLOCK_Y;

	/* load pop from global in cover nodes */
#include "includeFiles/popLoad.inc"

	dfloat invRho;

	if (nodeType != BULK)
	{
		boundary_calculation(nodeType, rho,
							 ux, uy,
							 mxx, myy, mxy,
							 pop);

		// boundary_calculation_old(nodeType, rho, ux, uy,
		// 						 mxx, myy, mxy, pop);

		invRho = 1.0 / rho;
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

	// COLLIDE
	moment_collision(ux, uy, mxx, myy, mxy);

	// calculate post collision populations
	pop_reconstruction(rho, ux, uy, mxx, myy, mxy, pop);

	/* write to global mom */

	state.d_rho[idxScalarBlock(threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y)] = rho - RHO_0;

	state.d_ux[idxScalarBlock(threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y)] = ux;
	state.d_uy[idxScalarBlock(threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y)] = uy;

	state.d_mxx[idxScalarBlock(threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y)] = mxx;
	state.d_mxy[idxScalarBlock(threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y)] = mxy;
	state.d_myy[idxScalarBlock(threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y)] = myy;

#include "includeFiles/popSave.inc"
}
