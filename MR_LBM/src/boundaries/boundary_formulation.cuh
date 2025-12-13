#pragma once
#include <stdint.h>
#include "../var.h"

namespace boundary
{
    __device__ void evaluate_dir(uint8_t node_type, uint8_t &incoming_mask, uint8_t &outgoing_mask);

    __device__ void eval_incoming_properties(uint8_t node_type, dfloat *pop,
                                             dfloat &rho, dfloat &ux, dfloat &uy,
                                             dfloat &mxx, dfloat &mxy, dfloat &myy);

    __device__ void numerical_solution(uint8_t node_type, dfloat OMEGA,
                                       dfloat mxx_prime, dfloat myy_prime,
                                       dfloat &rho, dfloat ux, dfloat uy,
                                       dfloat &mxx, dfloat &mxy, dfloat &myy);

    __device__ void eval_boundary(uint8_t node_type, dfloat OMEGA, dfloat *pop,
                                  dfloat &rho, dfloat ux, dfloat uy,
                                  dfloat &mxx, dfloat &mxy, dfloat &myy);

    __device__ void boundary_calculation(unsigned int nodeType, dfloat &rho, dfloat &ux,
                                         dfloat &uy, dfloat &mxx, dfloat &myy,
                                         dfloat &mxy, dfloat *pop,
                                         dfloat OMEGA);
}