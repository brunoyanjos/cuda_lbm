#include "mlbm.cuh"
#include "globalStructs.h"
#include "globalFunctions.h"

__global__ void gpuMomCollisionStream(
	latticeNode* nodes, ghostInterfaceData ghostInterface, unsigned int step
)
{
	const int x = threadIdx.x + blockDim.x * blockIdx.x;
	const int y = threadIdx.y + blockDim.y * blockIdx.y;

	if (x >= NX || y >= NY)
		return;
	dfloat pop[Q];
	__shared__ dfloat s_pop[BLOCK_LBM_SIZE * (Q - 1)];

	std::size_t idx = idxScalarBlock(threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y);

	// Load moments from global memory

	// rho'
	unsigned int nodeType = nodes[idx].node_type;
	if (nodeType == 0b11111111)
		return;
	dfloat rhoVar = RHO_0 + nodes[idx].rho;
	dfloat ux_t30 = nodes[idx].ux;
	dfloat uy_t30 = nodes[idx].uy;
	dfloat m_xx_t45 = nodes[idx].mxx;
	dfloat m_xy_t90 = nodes[idx].mxy;
	dfloat m_yy_t45 = nodes[idx].myy;

	pop_reconstruction(rhoVar, ux_t30, uy_t30, m_xx_t45, m_yy_t45, m_xy_t90, pop);

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
		boundary_calculation(nodeType, &rhoVar, &ux_t30, &uy_t30, &m_xx_t45, &m_yy_t45, &m_xy_t90, pop, nodes, x, y);

		invRho = 1.0 / rhoVar;
	}
	else
	{
		rhoVar = pop[0] + pop[1] + pop[2] + pop[3] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8];
		invRho = 1 / rhoVar;

		ux_t30 = ((pop[1] + pop[5] + pop[8]) - (pop[3] + pop[6] + pop[7])) * invRho;
		uy_t30 = ((pop[2] + pop[5] + pop[6]) - (pop[4] + pop[7] + pop[8])) * invRho;

		m_xx_t45 = (pop[1] + pop[3] + pop[5] + pop[6] + pop[7] + pop[8]) * invRho - cs2;
		m_xy_t90 = ((pop[5] + pop[7]) - (pop[6] + pop[8])) * invRho;
		m_yy_t45 = (pop[2] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8]) * invRho - cs2;
	}

	ux_t30 = F_M_I_SCALE * ux_t30;
	uy_t30 = F_M_I_SCALE * uy_t30;

	m_xx_t45 = F_M_II_SCALE * (m_xx_t45);
	m_xy_t90 = F_M_IJ_SCALE * (m_xy_t90);
	m_yy_t45 = F_M_II_SCALE * (m_yy_t45);

	// COLLIDE
	moment_collision(ux_t30, uy_t30, &m_xx_t45, &m_yy_t45, &m_xy_t90);

	// calculate post collision populations
	pop_reconstruction(rhoVar, ux_t30, uy_t30, m_xx_t45, m_yy_t45, m_xy_t90, pop);

	/* write to global mom */

	nodes[idx].rho = rhoVar - RHO_0;

	nodes[idx].ux = ux_t30;
	nodes[idx].uy = uy_t30;

	nodes[idx].mxx = m_xx_t45;
	nodes[idx].mxy = m_xy_t90;
	nodes[idx].myy = m_yy_t45;

#include "includeFiles/popSave.inc"
}
