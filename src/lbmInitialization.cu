#include "lbmInitialization.cuh"
#include <cmath>

__global__ void gpuInitialization_mom(LBMState state)
{
	int x = threadIdx.x + blockDim.x * blockIdx.x;
	int y = threadIdx.y + blockDim.y * blockIdx.y;
	if (x >= NX || y >= NY)
		return;

	// first moments
	dfloat rho, ux, uy;

	rho = RHO_0;
	ux = U_0_X;
	uy = U_0_Y;

	const uint8_t nodeType = state.d_node_type[idxBlock()];

	if (nodeType == SOLID_NODE)
		rho = static_cast<dfloat>(0);

	// zeroth moment
	state.d_rho[idxBlock()] = rho - RHO_0;
	state.d_ux[idxBlock()] = F_M_I_SCALE * ux;
	state.d_uy[idxBlock()] = F_M_I_SCALE * uy;

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

	state.d_mxx[idxBlock()] = F_M_II_SCALE * pixx;
	state.d_mxy[idxBlock()] = F_M_IJ_SCALE * pixy;
	state.d_myy[idxBlock()] = F_M_II_SCALE * piyy;
}

__global__ void gpuInitialization_pop(
	LBMState state, ghostInterfaceData ghostInterface)
{
	int x = threadIdx.x + blockDim.x * blockIdx.x;
	int y = threadIdx.y + blockDim.y * blockIdx.y;
	if (x >= NX || y >= NY)
		return;

	// zeroth moment
	dfloat rhoVar = RHO_0 + state.d_rho[idxBlock()];
	dfloat ux_t30 = state.d_ux[idxBlock()];
	dfloat uy_t30 = state.d_uy[idxBlock()];
	dfloat m_xx_t45 = state.d_mxx[idxBlock()];
	dfloat m_xy_t90 = state.d_mxy[idxBlock()];
	dfloat m_yy_t45 = state.d_myy[idxBlock()];

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

__global__ void gpuInitialization_nodeType(
	unsigned int *dNodeType)
{
	int x = threadIdx.x + blockDim.x * blockIdx.x;
	int y = threadIdx.y + blockDim.y * blockIdx.y;

	if (x >= NX || y >= NY)
		return;

	dNodeType[idxScalarBlock(threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y)] = boundary_definition(x, y);
}

__host__ void hostInitialization_nodeType_bulk(
	unsigned int *hNodeType)
{
	int x, y;

	for (y = 0; y < NY; y++)
	{
		for (x = 0; x < NX; x++)
		{
			hNodeType[idxScalarBlock(x % BLOCK_NX, y % BLOCK_NY, x / BLOCK_NX, y / BLOCK_NY)] = BULK;
		}
	}
}

__host__ void hostInitialization_nodeType(
	uint8_t *hNodeType)
{
	int x, y;

	for (y = 0; y < NY; y++)
	{
		for (x = 0; x < NX; x++)
		{
			hNodeType[idxScalarBlock(x % BLOCK_NX, y % BLOCK_NY, x / BLOCK_NX, y / BLOCK_NY)] = annular_boundary_definition(x, y);
		}
	}
}

__host__ inline uint8_t getNodeSafe(const uint8_t *node_type,
									int x, int y)
{
	// Se estiver fora do domínio, considere SOLID
	if (x < 0 || x >= NX || y < 0 || y >= NY)
		return SOLID_NODE;

	// caso contrário, acessa normalmente
	return node_type[idxBlockCoord(x, y)];
}

void initialize_boundaries(LBMState &state)
{

	for (int y = 0; y < NY; ++y)
	{
		for (int x = 0; x < NX; ++x)
		{
			uint8_t node_0 = getNodeSafe(state.h_node_type, x, y) == SOLID_NODE;

			uint8_t node_1 = getNodeSafe(state.h_node_type, x + 1, y) == BULK;
			uint8_t node_2 = getNodeSafe(state.h_node_type, x, y + 1) == BULK;
			uint8_t node_3 = getNodeSafe(state.h_node_type, x - 1, y) == BULK;
			uint8_t node_4 = getNodeSafe(state.h_node_type, x, y - 1) == BULK;

			uint8_t node_5 = getNodeSafe(state.h_node_type, x + 1, y + 1) == BULK;
			uint8_t node_6 = getNodeSafe(state.h_node_type, x - 1, y + 1) == BULK;
			uint8_t node_7 = getNodeSafe(state.h_node_type, x - 1, y - 1) == BULK;
			uint8_t node_8 = getNodeSafe(state.h_node_type, x + 1, y - 1) == BULK;

			const bool anyBulk = node_1 || node_2 || node_3 || node_4 || node_5 || node_6 || node_7 || node_8;

			const dfloat inner_radius = dfloat(D) / 2.0;
			const dfloat outer_radius = dfloat(NX - 1) / 2.0;

			const dfloat medium_radius = (inner_radius + outer_radius) / 2;

			const dfloat xc_local = xc;
			const dfloat yc_local = yc;

			const dfloat dx = dfloat(x) - xc_local;
			const dfloat dy = dfloat(y) - yc_local;
			const dfloat dist = dsqrt(dx * dx + dy * dy);

			if (node_0 && anyBulk)
			{
				const uint8_t bit_1 = !node_3 && !node_4 && !node_7 ? 0 : 1;
				const uint8_t bit_2 = !node_1 && !node_4 && !node_8 ? 0 : 1;
				const uint8_t bit_4 = !node_2 && !node_3 && !node_6 ? 0 : 1;
				const uint8_t bit_8 = !node_1 && !node_2 && !node_5 ? 0 : 1;

				uint8_t bc_number = bit_1 * NORTH_EAST +
									bit_2 * NORTH_WEST +
									bit_4 * SOUTH_EAST +
									bit_8 * SOUTH_WEST;

				if (dist < medium_radius)
					bc_number += INNER_BOUNDARY;

				state.h_node_type[idxBlockCoord(x, y)] = bc_number;
			}
		}
	}
}
