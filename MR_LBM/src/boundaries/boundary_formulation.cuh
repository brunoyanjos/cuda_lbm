#pragma once
#include <stdint.h>
#include "../var.h"

namespace boundary
{
    __device__ void evaluate_dir(uint8_t node_type, uint32_t &incoming_mask, uint32_t &outgoing_mask);

    __device__ void eval_incoming_properties(uint8_t node_type, dfloat *pop,
                                             dfloat &rho, dfloat &ux, dfloat &uy,
                                             dfloat &mxx, dfloat &mxy, dfloat &myy);
}