#include "domain.cuh"
#include "../boundaries/definition.cuh"
#include "../boundaries/node_type.h"
#include "../globalFunctions.h"
#include <iostream>

__global__ void initialize_boundaries(LBMState state)
{
    int x = threadIdx.x + blockDim.x * blockIdx.x;
    int y = threadIdx.y + blockDim.y * blockIdx.y;

    state.d_node_type[idxBlock()] = boundary::definition(x, y);
}

__global__ void initialize_moments(LBMState state)
{
    const dfloat rho = RHO_0;

    if (state.d_node_type[idxBlock()] == SOLID_NODE)
    {
        state.d_rho[idxBlock()] = dfloat(0.0);
    }
    else
    {
        state.d_rho[idxBlock()] = rho;
    }

    state.d_rho[idxBlock()] -= RHO_0;
}

void init_domain(LBMState &state)
{
    initialize_boundaries<<<gridBlock, threadBlock>>>(state);
    initialize_moments<<<gridBlock, threadBlock>>>(state);
}