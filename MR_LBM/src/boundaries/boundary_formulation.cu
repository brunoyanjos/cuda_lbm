#include "boundary_formulation.cuh"
#include "node_type.h"

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

    __device__ void numerical_solution(uint8_t node_type, dfloat OMEGA,
                                       dfloat mxx_prime, dfloat myy_prime,
                                       dfloat &rho, dfloat ux, dfloat uy,
                                       dfloat &mxx, dfloat &mxy, dfloat &myy)
    {
        const unsigned int x = threadIdx.x + blockDim.x * blockIdx.x;
        const unsigned int y = threadIdx.y + blockDim.y * blockIdx.y;

        uint8_t incoming_mask = 0;
        uint8_t outgoing_mask = 0;

        evaluate_dir(node_type, incoming_mask, outgoing_mask);

        const dfloat omega_var = static_cast<dfloat>(1) - OMEGA;

        dfloat Bxx = -static_cast<dfloat>(0.5) * as2 * w[0];
        dfloat Bxy = static_cast<dfloat>(0);
        dfloat Byy = -static_cast<dfloat>(0.5) * as2 * w[0];

        dfloat A = w[0];

        dfloat Bxx_Hxx = static_cast<dfloat>(0.5) * w[0];
        dfloat Bxy_Hxx = static_cast<dfloat>(0);
        dfloat Byy_Hxx = static_cast<dfloat>(0.5) * w[0];

        dfloat Bxx_Hxy = static_cast<dfloat>(0);
        dfloat Bxy_Hxy = static_cast<dfloat>(0);
        dfloat Byy_Hxy = static_cast<dfloat>(0);

        dfloat Bxx_Hyy = static_cast<dfloat>(0.5) * w[0];
        dfloat Bxy_Hyy = static_cast<dfloat>(0);
        dfloat Byy_Hyy = static_cast<dfloat>(0.5) * w[0];

        dfloat A_Hxx = -cs2 * w[0];
        dfloat A_Hxy = static_cast<dfloat>(0);
        dfloat A_Hyy = -cs2 * w[0];

        const dfloat x_diff = static_cast<dfloat>(x) - static_cast<dfloat>(xc);
        const dfloat y_diff = static_cast<dfloat>(y) - static_cast<dfloat>(yc);

        const dfloat radius = dsqrt(x_diff * x_diff + y_diff * y_diff);
        const dfloat inv_radius = static_cast<dfloat>(1) / radius;

        const dfloat cos_theta = x_diff * inv_radius;
        const dfloat sen_theta = y_diff * inv_radius;
        const dfloat sen_two_theta = static_cast<dfloat>(2.0) * cos_theta * sen_theta;
        const dfloat cos_two_theta = cos_theta * cos_theta - sen_theta * sen_theta;

        dfloat ux_prime = ux * cos_theta + uy * sen_theta;
        dfloat uy_prime = uy * cos_theta - ux * sen_theta;

#pragma unroll 8
        for (int i = 1; i < Q; ++i)
        {
            const dfloat cx_prime = cx[i] * cos_theta + cy[i] * sen_theta;
            const dfloat cy_prime = cy[i] * cos_theta - cx[i] * sen_theta;

            const dfloat Hxx = cx_prime * cx_prime - cs2;
            const dfloat Hxy = cx_prime * cy_prime;
            const dfloat Hyy = cy_prime * cy_prime - cs2;

            const dfloat A_i = w[i] * (1 + as2 * ux_prime * cx_prime + as2 * uy_prime * cy_prime);
            const dfloat Bxx_i = w[i] * as4 * static_cast<dfloat>(0.5) * Hxx;
            const dfloat Bxy_i = w[i] * as4 * static_cast<dfloat>(0.5) * Hxy;
            const dfloat Byy_i = w[i] * as4 * static_cast<dfloat>(0.5) * Hyy;

            if (outgoing_mask & (1u << (i - 1)))
            {
                Bxx += Bxx_i;
                Bxy += Bxy_i;
                Byy += Byy_i;

                A += A_i;
            }

            if (incoming_mask & (1u << (i - 1)))
            {
                Bxx_Hxx += Bxx_i * Hxx;
                Bxy_Hxx += Bxy_i * Hxx;
                Byy_Hxx += Byy_i * Hxx;

                Bxx_Hxy += Bxx_i * Hxy;
                Bxy_Hxy += Bxy_i * Hxy;
                Byy_Hxy += Byy_i * Hxy;

                Bxx_Hyy += Bxx_i * Hyy;
                Bxy_Hyy += Bxy_i * Hyy;
                Byy_Hyy += Byy_i * Hyy;

                A_Hxx += A_i * Hxx;
                A_Hxy += A_i * Hxy;
                A_Hyy += A_i * Hyy;
            }
        }

        const dfloat u_sum = ux_prime * ux_prime * Bxx +
                             static_cast<dfloat>(2) * ux_prime * uy_prime * Bxy +
                             uy_prime * uy_prime * Byy;

        const dfloat mxy_denominator = static_cast<dfloat>(2) * (omega_var * Bxy * mxy - Bxy_Hxy);

        const dfloat mxx_xy = omega_var * Bxx * mxy - Bxx_Hxy;
        const dfloat myy_xy = omega_var * Byy * mxy - Byy_Hxy;

        const dfloat xy_trace = mxx_xy * mxx_prime + myy_xy * myy_prime;

        const dfloat mxy_nominator = A_Hxy - (A + OMEGA * u_sum) * mxy - xy_trace;

        const dfloat mxy_prime = mxy_nominator / mxy_denominator;

        const dfloat mxx_factor = mxx_prime * Bxx;
        const dfloat mxy_factor = static_cast<dfloat>(2) * mxy_prime * Bxy;
        const dfloat myy_factor = myy_prime * Byy;

        const dfloat mom_sum = mxx_factor + mxy_factor + myy_factor;

        const dfloat rho_denominator = A + omega_var * mom_sum + OMEGA * u_sum;

        rho = rho / rho_denominator;

        mxx = mxx_prime * cos_theta * cos_theta +
              myy_prime * sen_theta * sen_theta -
              mxy_prime * sen_two_theta;

        mxy = static_cast<dfloat>(0.5) * (mxx_prime - myy_prime) * sen_two_theta +
              mxy_prime * cos_two_theta;

        myy = mxx_prime * sen_theta * sen_theta +
              myy_prime * cos_theta * cos_theta +
              mxy_prime * sen_two_theta;
    }

    __device__ void eval_boundary(uint8_t node_type, dfloat OMEGA, dfloat *pop,
                                  dfloat &rho, dfloat ux, dfloat uy,
                                  dfloat &mxx, dfloat &mxy, dfloat &myy)
    {
        uint8_t incoming_mask = 0;
        uint8_t outgoing_mask = 0;

        evaluate_dir(node_type, incoming_mask, outgoing_mask);

        const dfloat omega_var = static_cast<dfloat>(1) - OMEGA;

        dfloat rho_I = pop[0];

        dfloat mxx_I = -cs2 * pop[0];
        dfloat mxy_I = static_cast<dfloat>(0);
        dfloat myy_I = -cs2 * pop[0];

        dfloat A = w[0];

        dfloat A_Hxx = -cs2 * w[0];
        dfloat A_Hxy = static_cast<dfloat>(0);
        dfloat A_Hyy = -cs2 * w[0];

        dfloat Bxx = -static_cast<dfloat>(0.5) * as2 * w[0];
        dfloat Bxy = static_cast<dfloat>(0);
        dfloat Byy = -static_cast<dfloat>(0.5) * as2 * w[0];

        dfloat Bxx_Hxx = static_cast<dfloat>(0.5) * w[0];
        dfloat Bxy_Hxx = static_cast<dfloat>(0);
        dfloat Byy_Hxx = static_cast<dfloat>(0.5) * w[0];

        dfloat Bxx_Hxy = static_cast<dfloat>(0);
        dfloat Bxy_Hxy = static_cast<dfloat>(0);
        dfloat Byy_Hxy = static_cast<dfloat>(0);

        dfloat Bxx_Hyy = static_cast<dfloat>(0.5) * w[0];
        dfloat Bxy_Hyy = static_cast<dfloat>(0);
        dfloat Byy_Hyy = static_cast<dfloat>(0.5) * w[0];

#pragma unroll 8
        for (int i = 1; i < Q; ++i)
        {
            const dfloat Hxx = cx[i] * cx[i] - cs2;
            const dfloat Hxy = cx[i] * cy[i];
            const dfloat Hyy = cy[i] * cy[i] - cs2;

            const dfloat A_i = w[i] * (1 + as2 * ux * cx[i] + as2 * uy * cy[i]);
            const dfloat Bxx_i = w[i] * as4 * static_cast<dfloat>(0.5) * Hxx;
            const dfloat Bxy_i = w[i] * as4 * static_cast<dfloat>(0.5) * Hxy;
            const dfloat Byy_i = w[i] * as4 * static_cast<dfloat>(0.5) * Hyy;

            if (outgoing_mask & (1u << (i - 1)))
            {
                A += A_i;

                Bxx += Bxx_i;
                Bxy += Bxy_i;
                Byy += Byy_i;
            }

            if (incoming_mask & (1u << (i - 1)))
            {
                rho_I += pop[i];

                mxx_I += pop[i] * Hxx;
                mxy_I += pop[i] * Hxy;
                myy_I += pop[i] * Hyy;

                Bxx_Hxx += Bxx_i * Hxx;
                Bxy_Hxx += Bxy_i * Hxx;
                Byy_Hxx += Byy_i * Hxx;

                Bxx_Hxy += Bxx_i * Hxy;
                Bxy_Hxy += Bxy_i * Hxy;
                Byy_Hxy += Byy_i * Hxy;

                Bxx_Hyy += Bxx_i * Hyy;
                Bxy_Hyy += Bxy_i * Hyy;
                Byy_Hyy += Byy_i * Hyy;

                A_Hxx += A_i * Hxx;
                A_Hxy += A_i * Hxy;
                A_Hyy += A_i * Hyy;
            }
        }

        const dfloat u_sum = ux * ux * Bxx +
                             static_cast<dfloat>(2) * ux * uy * Bxy +
                             uy * uy * Byy;

        mxx = ux * ux;
        myy = uy * uy;

        const dfloat mxy_denominator = static_cast<dfloat>(2) * (omega_var * Bxy * mxy_I - Bxy_Hxy);

        const dfloat mxx_xy = omega_var * Bxx * mxy_I - Bxx_Hxy;
        const dfloat myy_xy = omega_var * Byy * mxy_I - Byy_Hxy;

        const dfloat xy_trace = mxx_xy * mxx + myy_xy * myy;

        const dfloat mxy_nominator = A_Hxy - (A + OMEGA * u_sum) * mxy_I - xy_trace;

        mxy = mxy_nominator / mxy_denominator;

        const dfloat mxx_factor = mxx * Bxx;
        const dfloat mxy_factor = static_cast<dfloat>(2) * mxy * Bxy;
        const dfloat myy_factor = myy * Byy;

        const dfloat mom_sum = mxx_factor + mxy_factor + myy_factor;

        const dfloat rho_denominator = A + omega_var * mom_sum + OMEGA * u_sum;

        rho = rho_I / rho_denominator;
    }

    __device__ void boundary_calculation(unsigned int nodeType, dfloat &rho, dfloat &ux,
                                         dfloat &uy, dfloat &mxx, dfloat &myy,
                                         dfloat &mxy, dfloat *pop,
                                         dfloat OMEGA)
    {
        switch (nodeType)
        {
        case NORTH:
        {
            const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[3] + pop[5] + pop[6];
            const dfloat inv_rhoIn = 1.0f / rhoIn;

            const dfloat mxyIn = (pop[5] - pop[6]) * inv_rhoIn;

            ux = U_MAX;
            uy = 0.0f;

            rho = 6.0f * rhoIn / 5.0f;

            mxx = U_MAX * U_MAX;
            mxy = 5.0f * mxyIn / 3.0f - U_MAX / 3.0f;
            myy = 0.0f;

            break;
        }
        case SOUTH:
        {
            const dfloat rhoIn = pop[0] + pop[1] + pop[3] + pop[4] + pop[7] + pop[8];
            const dfloat inv_rhoIn = 1.0f / rhoIn;

            const dfloat mxyIn = (pop[7] - pop[8]) * inv_rhoIn;

            ux = 0.0f;
            uy = 0.0f;

            rho = 6.0f * rhoIn / 5.0f;

            mxx = 0.0f;
            mxy = 5.0f * mxyIn / 3.0f;
            myy = 0.0f;

            break;
        }
        case WEST:
        {
            const dfloat rhoIn = pop[0] + pop[2] + pop[3] + pop[4] + pop[6] + pop[7];
            const dfloat inv_rhoIn = 1.0f / rhoIn;

            const dfloat mxyIn = (pop[7] - pop[6]) * inv_rhoIn;

            ux = 0.0f;
            uy = 0.0f;

            rho = 6.0f * rhoIn / 5.0f;

            mxx = 0.0f;
            mxy = 5.0f * mxyIn / 3.0f;
            myy = 0.0f;

            break;
        }
        case EAST:
        {
            const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[4] + pop[5] + pop[8];
            const dfloat inv_rhoIn = 1.0f / rhoIn;

            const dfloat mxyIn = (pop[5] - pop[8]) * inv_rhoIn;

            ux = 0.0f;
            uy = 0.0f;

            rho = 6.0f * rhoIn / 5.0f;

            mxx = 0.0f;
            mxy = 5.0f * mxyIn / 3.0f;
            myy = 0.0f;

            break;
        }
        case SOUTH_WEST:
        {
            const dfloat rhoIn = pop[0] + pop[3] + pop[4] + pop[7];
            const dfloat inv_rhoIn = 1.0f / rhoIn;

            const dfloat mxyIn = pop[7] * inv_rhoIn;

            ux = 0.0f;
            uy = 0.0f;

            rho = 36.0f * (rhoIn - mxyIn * rhoIn + mxyIn * OMEGA * rhoIn) / (24.0f + OMEGA);

            mxx = 0.0f;
            mxy = (36.0f * mxyIn * rhoIn - (rho)) / (9.0f * (rho));
            myy = 0.0f;

            break;
        }
        case SOUTH_EAST:
        {
            const dfloat rhoIn = pop[0] + pop[1] + pop[4] + pop[8];
            const dfloat inv_rhoIn = 1.0f / rhoIn;

            const dfloat mxyIn = -pop[8] * inv_rhoIn;

            ux = 0.0f;
            uy = 0.0f;

            rho = -36.0f * (mxyIn * OMEGA * rhoIn - rhoIn - mxyIn * rhoIn) / (24 + OMEGA);

            mxx = 0.0f;
            mxy = (36.0f * mxyIn * rhoIn + (rho)) / (9.0f * (rho));
            myy = 0.0f;

            break;
        }
        case NORTH_WEST:
        {
            const dfloat rhoIn = pop[0] + pop[2] + pop[3] + pop[6];
            const dfloat inv_rhoIn = 1.0f / rhoIn;

            const dfloat mxyIn = -pop[6] * inv_rhoIn;

            ux = U_MAX;
            uy = 0.0f;

            rho = -36.0f * (mxyIn * OMEGA * rhoIn - rhoIn - mxyIn * rhoIn) /
                  (24.0f + OMEGA + 18.0f * U_MAX - 3.0f * OMEGA * U_MAX - 18.0f * U_MAX * U_MAX + 3.0f * OMEGA * U_MAX * U_MAX);

            mxx = U_MAX * U_MAX;
            mxy = (36.0f * mxyIn * rhoIn + (rho)-3.0f * U_MAX * (rho) + 3.0f * U_MAX * U_MAX * (rho)) /
                  (9.0f * (rho));
            myy = 0.0f;

            break;
        }
        case NORTH_EAST:
        {
            const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[5];
            const dfloat inv_rhoIn = 1.0f / rhoIn;

            const dfloat mxyIn = pop[5] * inv_rhoIn;

            ux = U_MAX;
            uy = 0.0f;

            rho = 36.0f * (mxyIn * OMEGA * rhoIn + rhoIn - mxyIn * rhoIn) /
                  (24.0f + OMEGA - 18.0f * U_MAX + 3.0f * OMEGA * U_MAX - 18.0f * U_MAX * U_MAX + 3.0f * OMEGA * U_MAX * U_MAX);

            mxx = U_MAX * U_MAX;
            mxy = (36.0f * mxyIn * rhoIn - (rho)-3.0f * U_MAX * (rho)-3.0f * U_MAX * U_MAX * (rho)) /
                  (9.0f * (rho));
            myy = 0.0f;

            break;
        }
        default:
            break;
        }
    }
}