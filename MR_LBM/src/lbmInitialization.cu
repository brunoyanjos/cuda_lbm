#include "lbmInitialization.cuh"
#include <cmath>

__global__ void gpuInitialization_nodes(latticeNode *nodes)
{
	int x = threadIdx.x + blockDim.x * blockIdx.x;
	int y = threadIdx.y + blockDim.y * blockIdx.y;

	if (x >= NX || y >= NY)
		return;

	std::size_t idx = idxScalarBlock(threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y);

	// first moments
	dfloat rho, ux, uy;

	rho = RHO_0;
	ux = U_0_X;
	uy = U_0_Y;

	// zeroth moment
	nodes[idx].rho = rho - RHO_0;
	nodes[idx].ux = F_M_I_SCALE * ux;
	nodes[idx].uy = F_M_I_SCALE * uy;

	// second moments
	// define equilibrium populations
	dfloat pop[Q];
	for (int i = 0; i < Q; i++)
	{
		pop[i] = w[i] * RHO_0 * (1.0 + 3.0 * (ux * cx[i] + uy * cy[i]) + 4.5 * (ux * ux * (cx[i] * cx[i] - cs2) + uy * uy * (cx[i] * cx[i] - cs2)) + 9 * ux * uy * cx[i] * cy[i]);
	}

	dfloat invRho = 1.0 / rho;
	dfloat pixx = (pop[1] + pop[3] + pop[5] + pop[6] + pop[7] + pop[8]) * invRho - cs2;
	dfloat pixy = ((pop[5] + pop[7]) - (pop[6] + pop[8])) * invRho;
	dfloat piyy = (pop[2] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8]) * invRho - cs2;

	nodes[idx].mxx = F_M_II_SCALE * pixx;
	nodes[idx].mxy = F_M_IJ_SCALE * pixy;
	nodes[idx].myy = F_M_II_SCALE * piyy;
}

__global__ void gpuInitialization_pop(latticeNode *nodes, ghostInterfaceData ghostInterface)
{
	int x = threadIdx.x + blockDim.x * blockIdx.x;
	int y = threadIdx.y + blockDim.y * blockIdx.y;

	if (x >= NX || y >= NY)
		return;

	std::size_t idx = idxScalarBlock(threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y);

	// zeroth moment
	dfloat rhoVar = RHO_0 + nodes[idx].rho;
	dfloat ux_t30 = nodes[idx].ux;
	dfloat uy_t30 = nodes[idx].uy;
	dfloat m_xx_t45 = nodes[idx].mxx;
	dfloat m_xy_t90 = nodes[idx].mxy;
	dfloat m_yy_t45 = nodes[idx].myy;

	dfloat pop[Q];

	pop_reconstruction(rhoVar, ux_t30, uy_t30, m_xx_t45, m_yy_t45, m_xy_t90, pop);

	// thread xyz
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	// block xyz
	int bx = blockIdx.x;
	int by = blockIdx.y;

	if (threadIdx.x == 0)
	{ // w
		ghostInterface.fGhost.X_0[idxPopX(ty, 0, bx, by)] = pop[3];
		ghostInterface.fGhost.X_0[idxPopX(ty, 1, bx, by)] = pop[6];
		ghostInterface.fGhost.X_0[idxPopX(ty, 2, bx, by)] = pop[7];
	}
	else if (threadIdx.x == (BLOCK_NX - 1))
	{
		ghostInterface.fGhost.X_1[idxPopX(ty, 0, bx, by)] = pop[1];
		ghostInterface.fGhost.X_1[idxPopX(ty, 1, bx, by)] = pop[5];
		ghostInterface.fGhost.X_1[idxPopX(ty, 2, bx, by)] = pop[8];
	}

	if (threadIdx.y == 0)
	{ // s
		ghostInterface.fGhost.Y_0[idxPopY(tx, 0, bx, by)] = pop[4];
		ghostInterface.fGhost.Y_0[idxPopY(tx, 1, bx, by)] = pop[7];
		ghostInterface.fGhost.Y_0[idxPopY(tx, 2, bx, by)] = pop[8];
	}
	else if (threadIdx.y == (BLOCK_NY - 1))
	{
		ghostInterface.fGhost.Y_1[idxPopY(tx, 0, bx, by)] = pop[2];
		ghostInterface.fGhost.Y_1[idxPopY(tx, 1, bx, by)] = pop[5];
		ghostInterface.fGhost.Y_1[idxPopY(tx, 2, bx, by)] = pop[6];
	}
}

__global__ void gpuInitialization_nodeType_bulk(latticeNode *nodes)
{
	int x = threadIdx.x + blockDim.x * blockIdx.x;
	int y = threadIdx.y + blockDim.y * blockIdx.y;

	if (x >= NX || y >= NY)
		return;

	std::size_t idx = idxScalarBlock(threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y);

	nodes[idx].node_type = BULK;
}

__global__ void gpuInitialization_nodeType(latticeNode *nodes)
{
	int x = threadIdx.x + blockDim.x * blockIdx.x;
	int y = threadIdx.y + blockDim.y * blockIdx.y;

	if (x >= NX || y >= NY)
		return;

	unsigned int nodeType;

	std::size_t idx = idxScalarBlock(threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y);
	boundary_definition(&nodeType, x, y);

	if (nodeType != BULK)
		nodes[idx].node_type = (unsigned int)nodeType;
}