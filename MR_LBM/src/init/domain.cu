#include "domain.cuh"
#include "../boundaries/definition.cuh"
#include "../boundaries/node_type.h"
#include "../lbm/equilibrium.cuh"
#include "../colrec/collision_and_reconstruction.cuh"
#include "../globalFunctions.h"
#include "../interface/ghost_interface.cuh"
#include <iostream>

__device__ inline uint8_t getNodeSafe(const uint8_t *node_type,
                                      int x, int y)
{
    // Se estiver fora do domínio, considere SOLID
    if (x < 0 || x >= NX || y < 0 || y >= NY)
        return SOLID_NODE;

    // caso contrário, acessa normalmente
    return node_type[idxBlockCoord(x, y)];
}

__global__ void initialize_bulk(LBMState state)
{
    int x = threadIdx.x + blockDim.x * blockIdx.x;
    int y = threadIdx.y + blockDim.y * blockIdx.y;

    if (x >= NX || y >= NY)
        return;

    state.d_node_type[idxBlock()] = boundary::ldc_definition(x, y);
}

__global__ void initialize_boundaries(LBMState state)
{
    int x = threadIdx.x + blockDim.x * blockIdx.x;
    int y = threadIdx.y + blockDim.y * blockIdx.y;

    if (x >= NX || y >= NY)
        return;

    uint8_t node_0 = getNodeSafe(state.d_node_type, x, y) == SOLID_NODE;

    uint8_t node_1 = getNodeSafe(state.d_node_type, x + 1, y) == BULK;
    uint8_t node_2 = getNodeSafe(state.d_node_type, x, y + 1) == BULK;
    uint8_t node_3 = getNodeSafe(state.d_node_type, x - 1, y) == BULK;
    uint8_t node_4 = getNodeSafe(state.d_node_type, x, y - 1) == BULK;

    uint8_t node_5 = getNodeSafe(state.d_node_type, x + 1, y + 1) == BULK;
    uint8_t node_6 = getNodeSafe(state.d_node_type, x - 1, y + 1) == BULK;
    uint8_t node_7 = getNodeSafe(state.d_node_type, x + 1, y - 1) == BULK;
    uint8_t node_8 = getNodeSafe(state.d_node_type, x - 1, y - 1) == BULK;

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

        state.d_node_type[idxBlockCoord(x, y)] = bc_number;
    }
}

__global__ void initialize_moments(LBMState state)
{
    dfloat rho = RHO_0;
    dfloat inv_rho = 1.0 / rho;

    dfloat ux = dfloat(0.0);
    dfloat uy = dfloat(0.0);

    if (state.d_node_type[idxBlock()] == SOLID_NODE)
        rho -= RHO_0;

    if (state.d_node_type[idxBlock()] == NORTH ||
        state.d_node_type[idxBlock()] == NORTH_WEST ||
        state.d_node_type[idxBlock()] == NORTH_EAST)
    {
        ux = U_MAX;
    }

    dfloat pop[9];

    equilibrium(pop, rho, ux, uy);

    const dfloat mxx = (pop[1] + pop[3] + pop[5] + pop[6] + pop[7] + pop[8]) * inv_rho - cs2;
    const dfloat mxy = ((pop[5] + pop[7]) - (pop[6] + pop[8])) * inv_rho;
    const dfloat myy = (pop[2] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8]) * inv_rho - cs2;

    state.d_rho[idxBlock()] = rho - RHO_0;

    state.d_ux[idxBlock()] = ux * F_M_I_SCALE;
    state.d_uy[idxBlock()] = uy * F_M_I_SCALE;

    state.d_mxx[idxBlock()] = mxx * F_M_II_SCALE;
    state.d_mxy[idxBlock()] = mxy * F_M_IJ_SCALE;
    state.d_myy[idxBlock()] = myy * F_M_II_SCALE;
}

__global__ void gpuInitialization_pop(LBMState state, ghostInterfaceData ghostInterface)
{
    int x = threadIdx.x + blockDim.x * blockIdx.x;
    int y = threadIdx.y + blockDim.y * blockIdx.y;

    if (x >= NX || y >= NY)
        return;

    // zeroth moment
    dfloat rho = RHO_0 + state.d_rho[idxBlock()];

    dfloat ux = state.d_ux[idxBlock()];
    dfloat uy = state.d_uy[idxBlock()];

    dfloat mxx = state.d_mxx[idxBlock()];
    dfloat mxy = state.d_mxy[idxBlock()];
    dfloat myy = state.d_myy[idxBlock()];

    dfloat pop[Q];

    second::pop_reconstruction(pop, rho, ux, uy, mxx, mxy, myy);

    if (y == NY - 1)
    {
        for (int i = 0; y < 9; ++i)
        {
            printf("%f %f %f %f %f %f %f %f %f\n", pop[0],
                   pop[1], pop[2], pop[3],
                   pop[4], pop[5], pop[6],
                   pop[7], pop[8]);
        }
    }

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

__host__ void defining_geometry(LBMState state, dfloat &D_in, dfloat &D_out)
{
    for (int y = 0; y < NY; ++y)
    {
        for (int x = 0; x < NX; ++x)
        {
            uint8_t node_type = state.h_node_type[idxBlockCoord(x, y)];

            if (node_type != BULK && node_type != SOLID_NODE)
            {
                const dfloat inner_radius = dfloat(D) / 2.0;
                const dfloat outer_radius = dfloat(NX - 1) / 2.0;

                const dfloat medium_radius = (inner_radius + outer_radius) / 2;

                const dfloat xc_local = xc;
                const dfloat yc_local = yc;

                const dfloat dx = dfloat(x) - xc_local;
                const dfloat dy = dfloat(y) - yc_local;
                const dfloat dist = dsqrt(dx * dx + dy * dy);

                const dfloat value = 2 * dist;

                if (dist < medium_radius && value > D_in)
                {
                    D_in = value;
                }
                else if (dist > medium_radius && value < D_out)
                {
                    D_out = value;
                }
            }
        }
    }
}

void init_domain(LBMState &state, ghostInterfaceData &ghostInterface, dfloat &D_out, dfloat &D_in)
{
    initialize_bulk<<<gridBlock, threadBlock>>>(state);
    checkCudaErrors(cudaDeviceSynchronize());

    // initialize_boundaries<<<gridBlock, threadBlock>>>(state);
    // checkCudaErrors(cudaDeviceSynchronize());

    checkCudaErrors(cudaMemcpy(state.h_node_type, state.d_node_type, state.bytes_types, cudaMemcpyDeviceToHost));
    // defining_geometry(state, D_in, D_out);

    initialize_moments<<<gridBlock, threadBlock>>>(state);
    checkCudaErrors(cudaDeviceSynchronize());

    gpuInitialization_pop<<<gridBlock, threadBlock>>>(state, ghostInterface);
    checkCudaErrors(cudaDeviceSynchronize());

    interfaceCudaMemcpy(ghostInterface, ghostInterface.gGhost, ghostInterface.fGhost, cudaMemcpyDeviceToDevice, QF);
    checkCudaErrors(cudaDeviceSynchronize());
}