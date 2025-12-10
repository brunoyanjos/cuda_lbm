#pragma once

#include "../var.h"
#include "state.cuh"
#include "../globalStructs.h"

__device__ inline uint8_t getNodeSafe(const uint8_t *node_type, int x, int y);

__host__ void defining_geometry(LBMState state, dfloat &D_in, dfloat &D_out);

__host__ void init_domain(LBMState &state, ghostInterfaceData &ghostInterface, dfloat &D_out, dfloat &D_in);