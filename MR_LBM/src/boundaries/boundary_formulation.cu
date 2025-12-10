#include "boundary_formulation.cuh"

#include <cstdio>

namespace boundary
{
    __device__ void evaluate_dir(uint8_t node_type,
                                 uint8_t &incoming_mask,
                                 uint8_t &outgoing_mask)
    {
        const unsigned int out_dirs[4][3] = {
            {3, 4, 7}, // bit_1_dir
            {1, 4, 8}, // bit_2_dir
            {2, 3, 6}, // bit_4_dir
            {1, 2, 5}, // bit_8_dir
        };

        const unsigned int in_dirs[4][3] = {
            {1, 2, 5}, // bit_1_dir
            {2, 3, 6}, // bit_2_dir
            {1, 4, 8}, // bit_4_dir
            {3, 4, 7}, // bit_8_dir
        };

#pragma unroll 4
        for (unsigned int j = 1; j <= 8; j <<= 1)
        {
            if (node_type & j)
            {
                unsigned int idx = __ffs(j) - 1;

#pragma unroll 3
                for (unsigned int k = 0; k < 3; ++k)
                {
                    unsigned int od = out_dirs[idx][k] - 1;
                    unsigned int id = in_dirs[idx][k] - 1;

                    // liga o bit correspondente à direção
                    outgoing_mask |= (1u << od);
                    incoming_mask |= (1u << id);
                }
            }
        }
    }

    __device__ void eval_incoming_properties(uint8_t node_type, dfloat *pop,
                                             dfloat &rho, dfloat &ux, dfloat &uy,
                                             dfloat &mxx, dfloat &mxy, dfloat &myy)
    {
        const unsigned int x = threadIdx.x + blockDim.x * blockIdx.x;
        const unsigned int y = threadIdx.y + blockDim.y * blockIdx.y;

        uint8_t incoming_mask = 0;
        uint8_t outgoing_mask = 0;

        evaluate_dir(node_type, incoming_mask, outgoing_mask);

        dfloat rho_I = pop[0];

        dfloat mxx_I = -cs2 * pop[0];
        dfloat mxy_I = static_cast<dfloat>(0);
        dfloat myy_I = -cs2 * pop[0];

        const dfloat x_diff = static_cast<dfloat>(x) - static_cast<dfloat>(xc);
        const dfloat y_diff = static_cast<dfloat>(y) - static_cast<dfloat>(yc);

        const dfloat radius = dsqrt(x_diff * x_diff + y_diff * y_diff);
        const dfloat inv_radius = static_cast<dfloat>(1) / radius;

        const dfloat cos_theta = x_diff * inv_radius;
        const dfloat sen_theta = y_diff * inv_radius;

#pragma unroll 8
        for (int i = 1; i < 9; ++i)
        {
            const dfloat cx_prime = cx[i] * cos_theta + cy[i] * sen_theta;
            const dfloat cy_prime = cy[i] * cos_theta - cx[i] * sen_theta;

            const dfloat Hxx = cx_prime * cx_prime - cs2;
            const dfloat Hxy = cx_prime * cy_prime;
            const dfloat Hyy = cy_prime * cy_prime - cs2;

            if (incoming_mask & (1u << (i - 1)))
            {
                rho_I += pop[i];

                mxx_I += pop[i] * Hxx;
                mxy_I += pop[i] * Hxy;
                myy_I += pop[i] * Hyy;
            }
        }

        const dfloat inv_rho_I = static_cast<dfloat>(1) / rho_I;

        rho = rho_I;

        mxx = mxx_I * inv_rho_I;
        mxy = mxy_I * inv_rho_I;
        myy = myy_I * inv_rho_I;
    }
}