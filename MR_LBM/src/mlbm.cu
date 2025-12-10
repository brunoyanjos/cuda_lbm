#include "mlbm.cuh"
#include "globalStructs.h"
#include "globalFunctions.h"
#include "init/state.cuh"
#include "boundaries/node_type.h"
#include "boundaries/boundary_formulation.cuh"
#include "colrec/collision_and_reconstruction.cuh"
#include "interface/interface_handling.cuh"

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

__global__ void updateInnerBoundaries(dfloat *fMom, dfloat OMEGA, unsigned int step)
{
	// cylinderProperties property = cylinder_properties[threadIdx.x];

	// const int xb = (int)property.xb;
	// const int yb = (int)property.yb;

	// const int tx = xb % BLOCK_NX;
	// const int ty = yb % BLOCK_NY;

	// const int bx = xb / BLOCK_NX;
	// const int by = yb / BLOCK_NY;

	// dfloat rhoVar = RHO_0 + fMom[idxMom(tx, ty, M_RHO_INDEX, bx, by)];

	// dfloat ux_t30 = fMom[idxMom(tx, ty, M_UX_INDEX, bx, by)];
	// dfloat uy_t30 = fMom[idxMom(tx, ty, M_UY_INDEX, bx, by)];

	// dfloat m_xx_t45 = fMom[idxMom(tx, ty, M_MXX_INDEX, bx, by)];
	// dfloat m_xy_t90 = fMom[idxMom(tx, ty, M_MXY_INDEX, bx, by)];
	// dfloat m_yy_t45 = fMom[idxMom(tx, ty, M_MYY_INDEX, bx, by)];

	// if (property.isBulk == false)
	// {
	// for first point

	// dfloat ux1;
	// dfloat uy1;

	// bilinear_velocity_interpolation(property.x1, property.y1, int(property.x1), int(property.y1),
	// 								int(property.x1) + 1, int(property.y1) + 1, fMom, M_UX_INDEX, &ux1);
	// bilinear_velocity_interpolation(property.x1, property.y1, int(property.x1), int(property.y1),
	// 								int(property.x1) + 1, int(property.y1) + 1, fMom, M_UY_INDEX, &uy1);

	// // for second point

	// dfloat ux2;
	// dfloat uy2;

	// bilinear_velocity_interpolation(property.x2, property.y2, int(property.x2), int(property.y2),
	// 								int(property.x2) + 1, int(property.y2) + 1, fMom, M_UX_INDEX, &ux2);
	// bilinear_velocity_interpolation(property.x2, property.y2, int(property.x2), int(property.y2),
	// 								int(property.x2) + 1, int(property.y2) + 1, fMom, M_UY_INDEX, &uy2);

	// // moment interpolation to first point
	// dfloat mxx1 = 0.0;
	// dfloat myy1 = 0.0;
	// dfloat mxx2 = 0.0;
	// dfloat myy2 = 0.0;

	// if (ROTATIONAL_COORDINATES)
	// {
	// 	bilinear_moment_interpolation(property.x1, property.y1, int(property.x1), int(property.y1),
	// 								  int(property.x1) + 1, int(property.y1) + 1, fMom, &mxx1, &myy1);
	// 	bilinear_moment_interpolation(property.x2, property.y2, int(property.x2), int(property.y2),
	// 								  int(property.x2) + 1, int(property.y2) + 1, fMom, &mxx2, &myy2);
	// }

	// if (CALCULATE_PRESSURE && step >= STAT_BEGIN_TIME && step <= STAT_END_TIME)
	// {
	// 	dfloat rho1;
	// 	dfloat rho2;
	// 	dfloat rho3;

	// 	bilinear_density_interpolation(property.x1, property.y1,
	// 								   int(property.x1), int(property.y1),
	// 								   int(property.x1) + 1, int(property.y1) + 1,
	// 								   fMom, M_RHO_INDEX, &rho1);
	// 	bilinear_density_interpolation(property.x2, property.y2,
	// 								   int(property.x2), int(property.y2),
	// 								   int(property.x2) + 1, int(property.y2) + 1,
	// 								   fMom, M_RHO_INDEX, &rho2);
	// 	bilinear_density_interpolation(property.x3, property.y3,
	// 								   int(property.x3), int(property.y3),
	// 								   int(property.x3) + 1, int(property.y3) + 1,
	// 								   fMom, M_RHO_INDEX, &rho3);

	// 	pressure_extrapolation_old(property.xw, property.yw,
	// 							   property.x1, property.y1,
	// 							   property.x2, property.y2,
	// 							   property.x3, property.y3,
	// 							   rho1, rho2, rho3, &(cylinder_properties[threadIdx.x].ps));
	// }

	// const dfloat delta = property.dr;

	// ux_t30 = extrapolation(delta, ux1, ux2);
	// uy_t30 = extrapolation(delta, uy1, uy2);

	// const dfloat m_xx_int = extrapolation(delta, mxx1, mxx2);
	// const dfloat m_yy_int = extrapolation(delta, myy1, myy2);

	// if (ROTATIONAL_COORDINATES)
	// {
	// 	if (RHO_STRONG)
	// 	{
	// 		numericalSolution_rotation(&rhoVar, ux_t30, uy_t30, &m_xx_t45, &m_xy_t90, &m_yy_t45, m_xx_int, m_yy_int, property.is, property.os, OMEGA, xb, yb);
	// 	}
	// 	if (RHO_EQ)
	// 	{
	// 		numericalSolution_rotation_rhoeq(&rhoVar, ux_t30, uy_t30, &m_xx_t45, &m_xy_t90, &m_yy_t45, m_xx_int, m_yy_int, property.is, property.os, OMEGA, xb, yb);
	// 	}
	// }
	// else
	// {
	// 	numericalSolution(&rhoVar, ux_t30, uy_t30, &m_xx_t45, &m_xy_t90, &m_yy_t45, property.is, property.os, OMEGA);
	// }
	// }

	// fMom[idxMom(tx, ty, M_RHO_INDEX, bx, by)] = rhoVar - RHO_0;

	// fMom[idxMom(tx, ty, M_UX_INDEX, bx, by)] = ux_t30;
	// fMom[idxMom(tx, ty, M_UY_INDEX, bx, by)] = uy_t30;

	// fMom[idxMom(tx, ty, M_MXX_INDEX, bx, by)] = m_xx_t45;
	// fMom[idxMom(tx, ty, M_MXY_INDEX, bx, by)] = m_xy_t90;
	// fMom[idxMom(tx, ty, M_MYY_INDEX, bx, by)] = m_yy_t45;
}

__global__ void boundaryAndCollision(
	dfloat *fMom, size_t cylinder_count, dfloat OMEGA, unsigned int *dNodeType,
	ghostInterfaceData ghostInterface, unsigned int step)
{
	const int x = threadIdx.x + blockDim.x * blockIdx.x;
	const int y = threadIdx.y + blockDim.y * blockIdx.y;

	if (x >= NX || y >= NY)
		return;
	// dfloat pop[Q];

	// Load moments from global memory

	// rho'
	unsigned int nodeType = dNodeType[idxScalarBlock(threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y)];
	if (nodeType == 0b11111111)
		return;
	dfloat rhoVar = RHO_0 + fMom[idxMom(threadIdx.x, threadIdx.y, M_RHO_INDEX, blockIdx.x, blockIdx.y)];
	dfloat ux_t30 = fMom[idxMom(threadIdx.x, threadIdx.y, M_UX_INDEX, blockIdx.x, blockIdx.y)];
	dfloat uy_t30 = fMom[idxMom(threadIdx.x, threadIdx.y, M_UY_INDEX, blockIdx.x, blockIdx.y)];
	dfloat m_xx_t45 = fMom[idxMom(threadIdx.x, threadIdx.y, M_MXX_INDEX, blockIdx.x, blockIdx.y)];
	dfloat m_xy_t90 = fMom[idxMom(threadIdx.x, threadIdx.y, M_MXY_INDEX, blockIdx.x, blockIdx.y)];
	dfloat m_yy_t45 = fMom[idxMom(threadIdx.x, threadIdx.y, M_MYY_INDEX, blockIdx.x, blockIdx.y)];

	ux_t30 = F_M_I_SCALE * ux_t30;
	uy_t30 = F_M_I_SCALE * uy_t30;

	m_xx_t45 = F_M_II_SCALE * (m_xx_t45);
	m_xy_t90 = F_M_IJ_SCALE * (m_xy_t90);
	m_yy_t45 = F_M_II_SCALE * (m_yy_t45);

	// moment_collision(ux_t30, uy_t30, &m_xx_t45, &m_yy_t45, &m_xy_t90, OMEGA);

	// pop_reconstruction(rhoVar, ux_t30, uy_t30, m_xx_t45, m_yy_t45, m_xy_t90, pop);

	// if (nodeType >= 100 && step >= STAT_BEGIN_TIME && step <= STAT_END_TIME && CALCULATE_FORCES)
	// {
	// 	cylinderProperties *bc_property = findCylindeProperty(cylinder_properties, cylinder_count, x, y);

	// 	outgoing_forces(bc_property, cylinder_count, pop);
	// }

	// fMom[idxMom(threadIdx.x, threadIdx.y, M_RHO_INDEX, blockIdx.x, blockIdx.y)] = rhoVar - RHO_0;

	// fMom[idxMom(threadIdx.x, threadIdx.y, M_UX_INDEX, blockIdx.x, blockIdx.y)] = ux_t30;
	// fMom[idxMom(threadIdx.x, threadIdx.y, M_UY_INDEX, blockIdx.x, blockIdx.y)] = uy_t30;

	// fMom[idxMom(threadIdx.x, threadIdx.y, M_MXX_INDEX, blockIdx.x, blockIdx.y)] = m_xx_t45;
	// fMom[idxMom(threadIdx.x, threadIdx.y, M_MXY_INDEX, blockIdx.x, blockIdx.y)] = m_xy_t90;
	// fMom[idxMom(threadIdx.x, threadIdx.y, M_MYY_INDEX, blockIdx.x, blockIdx.y)] = m_yy_t45;

	// pop_save(ghostInterface, threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y, x, y, pop);
}