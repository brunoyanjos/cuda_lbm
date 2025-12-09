#pragma once

#include <cstdint>
#include <cstddef>

#include "../var.h"

struct LBMState
{
    // DEVICE
    uint8_t *d_node_type;

    dfloat *d_rho;

    dfloat *d_ux;
    dfloat *d_uy;

    dfloat *d_mxx;
    dfloat *d_mxy;
    dfloat *d_myy;

    // HOST
    uint8_t *h_node_type;
    dfloat *h_rho;

    dfloat *h_ux;
    dfloat *h_uy;

    dfloat *h_mxx;
    dfloat *h_mxy;
    dfloat *h_myy;

    // SIZES
    size_t bytes_fields;
    size_t bytes_types;
};

[[nodiscard]] __host__ LBMState init_state();

__host__ void upload_state_to_host(LBMState &state);

__host__ void free_state(LBMState &state);