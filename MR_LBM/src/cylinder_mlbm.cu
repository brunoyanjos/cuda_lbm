#include "mlbm.cuh"
#include "globalStructs.h"
#include "globalFunctions.h"

__global__ void streamingAndMom(
    dfloat *fMom, dfloat OMEGA, unsigned int *dNodeType,
    ghostInterfaceData ghostInterface, cylinderProperties *cylinder_properties, unsigned int step)
{
    const int x = threadIdx.x + blockDim.x * blockIdx.x;
    const int y = threadIdx.y + blockDim.y * blockIdx.y;

    if (x >= NX || y >= NY)
        return;

    dfloat pop[Q];
    dfloat pics2;
    dfloat multiplyTerm;
    __shared__ dfloat s_pop[BLOCK_LBM_SIZE * (Q - 1)];

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

    pop_reconstruction(rhoVar, ux_t30, uy_t30, m_xx_t45, m_xy_t90, m_yy_t45, pop);

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
        boundary_calculation(nodeType, &rhoVar, &ux_t30, &uy_t30, &m_xx_t45, &m_yy_t45, &m_xy_t90, pop, fMom, x, y, OMEGA);
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

    // if (nodeType > 100 && step >= N_STEPS - FORCES_TIME && CALCULATE_FORCES) {
    // 	incoming_forces(nodeType, x, y, cylinder_properties, cylinder_index, pop);
    // }

    fMom[idxMom(threadIdx.x, threadIdx.y, M_RHO_INDEX, blockIdx.x, blockIdx.y)] = rhoVar - RHO_0;

    fMom[idxMom(threadIdx.x, threadIdx.y, M_UX_INDEX, blockIdx.x, blockIdx.y)] = ux_t30;
    fMom[idxMom(threadIdx.x, threadIdx.y, M_UY_INDEX, blockIdx.x, blockIdx.y)] = uy_t30;

    fMom[idxMom(threadIdx.x, threadIdx.y, M_MXX_INDEX, blockIdx.x, blockIdx.y)] = m_xx_t45;
    fMom[idxMom(threadIdx.x, threadIdx.y, M_MXY_INDEX, blockIdx.x, blockIdx.y)] = m_xy_t90;
    fMom[idxMom(threadIdx.x, threadIdx.y, M_MYY_INDEX, blockIdx.x, blockIdx.y)] = m_yy_t45;
}

__global__ void boundaryAndCollision(
    dfloat *fMom, dfloat *fMom_old, dfloat OMEGA, unsigned int *dNodeType,
    ghostInterfaceData ghostInterface, cylinderProperties *cylinder_properties, unsigned int step)
{
    const int x = threadIdx.x + blockDim.x * blockIdx.x;
    const int y = threadIdx.y + blockDim.y * blockIdx.y;

    if (x >= NX || y >= NY)
        return;
    dfloat pop[Q];

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

    // if (nodeType > 100) {
    // 	immersed_boundary_treatment(
    // 		nodeType,
    // 		&rhoVar, &ux_t30, &uy_t30, &m_xx_t45, &m_xy_t90, &m_yy_t45,
    // 		threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y, x, y,
    // 		cylinder_properties,
    // 		fMom_old, OMEGA, step);
    // }

    ux_t30 = F_M_I_SCALE * ux_t30;
    uy_t30 = F_M_I_SCALE * uy_t30;

    m_xx_t45 = F_M_II_SCALE * (m_xx_t45);
    m_xy_t90 = F_M_IJ_SCALE * (m_xy_t90);
    m_yy_t45 = F_M_II_SCALE * (m_yy_t45);

    moment_collision(ux_t30, uy_t30, &m_xx_t45, &m_xy_t90, &m_yy_t45, OMEGA);

    pop_reconstruction(rhoVar, ux_t30, uy_t30, m_xx_t45, m_xy_t90, m_yy_t45, pop);

    // if (nodeType > 100 && step >= N_STEPS - FORCES_TIME && CALCULATE_FORCES) {
    // 	outgoing_forces(nodeType, x, y, cylinder_properties, cylinder_index, pop);
    // }

    fMom[idxMom(threadIdx.x, threadIdx.y, M_RHO_INDEX, blockIdx.x, blockIdx.y)] = rhoVar - RHO_0;

    fMom[idxMom(threadIdx.x, threadIdx.y, M_UX_INDEX, blockIdx.x, blockIdx.y)] = ux_t30;
    fMom[idxMom(threadIdx.x, threadIdx.y, M_UY_INDEX, blockIdx.x, blockIdx.y)] = uy_t30;

    fMom[idxMom(threadIdx.x, threadIdx.y, M_MXX_INDEX, blockIdx.x, blockIdx.y)] = m_xx_t45;
    fMom[idxMom(threadIdx.x, threadIdx.y, M_MXY_INDEX, blockIdx.x, blockIdx.y)] = m_xy_t90;
    fMom[idxMom(threadIdx.x, threadIdx.y, M_MYY_INDEX, blockIdx.x, blockIdx.y)] = m_yy_t45;

    pop_save(ghostInterface, threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y, x, y, pop);
}