#ifndef LBM_STEPS_CUH
#define LBM_STEPS_CUH

#include "var.h"
#include "newton_raphson.cuh"

__host__ inline void init_pop_eq(dfloat *pop_in,
                                 dfloat rho, dfloat ux, dfloat uy)
{
    dfloat pics2 = 1 - cs2 * (ux * ux + uy * uy);

    dfloat multiplyTerm = W0 * rho;
    pop_in[0] = multiplyTerm * (pics2);

    multiplyTerm = W1 * rho;
    pop_in[1] = multiplyTerm * (pics2 + ux + static_cast<dfloat>(0.5) * ux * ux);
    pop_in[2] = multiplyTerm * (pics2 + uy + static_cast<dfloat>(0.5) * uy * uy);
    pop_in[3] = multiplyTerm * (pics2 - ux + static_cast<dfloat>(0.5) * ux * ux);
    pop_in[4] = multiplyTerm * (pics2 - uy + static_cast<dfloat>(0.5) * uy * uy);

    multiplyTerm = W2 * rho;
    pop_in[5] = multiplyTerm * (pics2 + ux + uy + static_cast<dfloat>(0.5) * ux * ux + static_cast<dfloat>(0.5) * uy * uy + ux * uy);
    pop_in[6] = multiplyTerm * (pics2 - ux + uy + static_cast<dfloat>(0.5) * ux * ux + static_cast<dfloat>(0.5) * uy * uy - ux * uy);
    pop_in[7] = multiplyTerm * (pics2 - ux - uy + static_cast<dfloat>(0.5) * ux * ux + static_cast<dfloat>(0.5) * uy * uy + ux * uy);
    pop_in[8] = multiplyTerm * (pics2 + ux - uy + static_cast<dfloat>(0.5) * ux * ux + static_cast<dfloat>(0.5) * uy * uy - ux * uy);
}

__host__ [[nodiscard]] inline arrayType<9> regularization_mass(const latticeNode &node) noexcept
{
    const dfloat rho = node.rho;
    const dfloat ux = node.ux;
    const dfloat uy = node.uy;
    const dfloat mxx = node.mxx;
    const dfloat mxy = node.mxy;
    const dfloat myy = node.myy;

    arrayType<9> pop;

    const dfloat pics2 = 1 - cs2 * (mxx + myy);
    dfloat multiplyTerm = W0 * rho;

    pop.f[0] = multiplyTerm * (pics2);
    multiplyTerm = W1 * rho;
    pop.f[1] = multiplyTerm * (pics2 + ux + mxx);
    pop.f[2] = multiplyTerm * (pics2 + uy + myy);
    pop.f[3] = multiplyTerm * (pics2 - ux + mxx);
    pop.f[4] = multiplyTerm * (pics2 - uy + myy);

    multiplyTerm = W2 * rho;
    pop.f[5] = multiplyTerm * (pics2 + ux + uy + mxx + myy + mxy);
    pop.f[6] = multiplyTerm * (pics2 - ux + uy + mxx + myy - mxy);
    pop.f[7] = multiplyTerm * (pics2 - ux - uy + mxx + myy + mxy);
    pop.f[8] = multiplyTerm * (pics2 + ux - uy + mxx + myy - mxy);

    return pop;
}

__host__ inline void regularization(
    dfloat *pop_out,
    dfloat rho, dfloat ux, dfloat uy,
    dfloat mxx, dfloat mxy, dfloat myy)
{
    dfloat pics2 = 1 - cs2 * (mxx + myy);

    dfloat multiplyTerm = W0 * rho;
    pop_out[0] = multiplyTerm * (pics2);

    multiplyTerm = W1 * rho;
    pop_out[1] = multiplyTerm * (pics2 + ux + mxx);
    pop_out[2] = multiplyTerm * (pics2 + uy + myy);
    pop_out[3] = multiplyTerm * (pics2 - ux + mxx);
    pop_out[4] = multiplyTerm * (pics2 - uy + myy);

    multiplyTerm = W2 * rho;
    pop_out[5] = multiplyTerm * (pics2 + ux + uy + mxx + myy + mxy);
    pop_out[6] = multiplyTerm * (pics2 - ux + uy + mxx + myy - mxy);
    pop_out[7] = multiplyTerm * (pics2 - ux - uy + mxx + myy + mxy);
    pop_out[8] = multiplyTerm * (pics2 + ux - uy + mxx + myy - mxy);
}

__host__ inline void
boundary_condition(unsigned int node_type, const dfloat *pop,
                   dfloat *rho, dfloat *ux, dfloat *uy,
                   dfloat *mxx, dfloat *mxy, dfloat *myy,
                   dfloat omega)
{
    switch (node_type)
    {
    case NORTH:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[3] + pop[5] + pop[6];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat mxxIn = (pop[1] + pop[3] + pop[5] + pop[6]) * inv_rhoIn - cs2;
        const dfloat mxyIn = (pop[5] - pop[6]) * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[5] + pop[6]) * inv_rhoIn - cs2;

        *ux = U_MAX;
        *uy = 0.0f;

        const dfloat rhoVar = 6.0f * rhoIn / 5.0f;

        *mxx = U_MAX * U_MAX;
        *mxx = 0.0f;
        *mxy = 5.0f * mxyIn / 3.0f - U_MAX / 3.0f;
        *myy = 0.0f;

        *rho = rhoVar;

        break;
    }
    case SOUTH:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[3] + pop[4] + pop[7] + pop[8];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat mxxIn = (pop[1] + pop[3] + pop[7] + pop[8]) * inv_rhoIn - cs2;
        const dfloat mxyIn = (pop[7] - pop[8]) * inv_rhoIn;
        const dfloat myyIn = (pop[3] + pop[7] + pop[8]) * inv_rhoIn - cs2;

        *ux = 0.0f;
        *uy = 0.0f;

        const dfloat rhoVar = 6.0f * rhoIn / 5.0f;

        *mxx = 0.0f;
        *mxy = 5.0f * mxyIn / 3.0f;
        *myy = 0.0f;

        *rho = rhoVar;

        break;
    }
    case WEST:
    {
        const dfloat rhoIn = pop[0] + pop[2] + pop[3] + pop[4] + pop[6] + pop[7];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat mxxIn = (pop[3] + pop[6] + pop[7]) * inv_rhoIn - cs2;
        const dfloat mxyIn = (pop[7] - pop[6]) * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[4] + pop[6] + pop[7]) * inv_rhoIn - cs2;

        *ux = 0.0f;
        *uy = 0.0f;

        const dfloat rhoVar = 6.0f * rhoIn / 5.0f;

        *mxx = 0.0f;
        *mxy = 5.0f * mxyIn / 3.0f;
        *myy = 0.0f;

        *rho = rhoVar;

        break;
    }
    case EAST:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[4] + pop[5] + pop[8];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat mxxIn = (pop[1] + pop[5] + pop[8]) * inv_rhoIn - cs2;
        const dfloat mxyIn = (pop[5] - pop[8]) * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[4] + pop[5] + pop[8]) * inv_rhoIn - cs2;

        *ux = 0.0f;
        *uy = 0.0f;

        const dfloat rhoVar = 6.0f * rhoIn / 5.0f;

        *mxx = 0.0f;
        *mxy = 5.0f * mxyIn / 3.0f;
        *myy = 0.0f;

        *rho = rhoVar;

        break;
    }
    case SOUTH_WEST:
    {
        const dfloat rhoIn = pop[0] + pop[3] + pop[4] + pop[7];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat mxxIn = (pop[3] + pop[7]) * inv_rhoIn - cs2;
        const dfloat mxyIn = pop[7] * inv_rhoIn;
        const dfloat myyIn = (pop[4] + pop[7]) * inv_rhoIn - cs2;

        *ux = 0.0f;
        *uy = 0.0f;

        const dfloat rhoVar = 36.0f * (rhoIn - mxyIn * rhoIn + mxyIn * omega * rhoIn) /
                              (24.0f + omega);

        *mxx = 0.0f;
        *mxy = (36.0f * mxyIn * rhoIn - (rhoVar)) /
               (9.0f * (rhoVar));
        *myy = 0.0f;

        *rho = rhoVar;

        break;
    }
    case SOUTH_EAST:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[4] + pop[8];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat mxxIn = (pop[1] + pop[8]) * inv_rhoIn - cs2;
        const dfloat mxyIn = -pop[8] * inv_rhoIn;
        const dfloat myyIn = (pop[4] + pop[8]) * inv_rhoIn - cs2;

        *uy = 0.0f;
        *ux = 0.0f;

        const dfloat rhoVar = -36.0f * (mxyIn * omega * rhoIn - rhoIn - mxyIn * rhoIn) /
                              (24 + omega);

        *mxx = 0.0f;
        *mxy = (36.0f * mxyIn * rhoIn + (rhoVar)) / (9.0f * (rhoVar));
        *myy = 0.0f;

        *rho = rhoVar;

        break;
    }
    case NORTH_WEST:
    {
        const dfloat rhoIn = pop[0] + pop[2] + pop[3] + pop[6];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat mxxIn = (pop[3] + pop[6]) * inv_rhoIn - cs2;
        const dfloat mxyIn = -pop[6] * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[6]) * inv_rhoIn - cs2;

        *ux = U_MAX;
        *uy = 0.0f;

        const dfloat rhoVar = -36.0f * (mxyIn * omega * rhoIn - rhoIn - mxyIn * rhoIn) /
                              (24.0f + omega + 18.0f * U_MAX - 3.0f * omega * U_MAX - 18.0f * U_MAX * U_MAX + 3.0f * omega * U_MAX * U_MAX);

        *mxx = U_MAX * U_MAX;
        *mxy = (36.0f * mxyIn * rhoIn + (rhoVar)-3.0f * U_MAX * (rhoVar) + 3.0f * U_MAX * U_MAX * (rhoVar)) / (9.0f * (rhoVar));
        *myy = 0.0f;

        *rho = rhoVar;

        break;
    }
    case NORTH_EAST:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[5];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat mxxIn = (pop[1] + pop[5]) * inv_rhoIn - cs2;
        const dfloat mxyIn = pop[5] * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[5]) * inv_rhoIn - cs2;

        *ux = U_MAX;
        *uy = 0.0f;

        const dfloat rhoVar = 36.0f * (mxyIn * omega * rhoIn + rhoIn - mxyIn * rhoIn) /
                              (24.0f + omega - 18.0f * U_MAX + 3.0f * omega * U_MAX - 18.0f * U_MAX * U_MAX + 3.0f * omega * U_MAX * U_MAX);

        *mxx = U_MAX * U_MAX;
        *mxy = (36.0f * mxyIn * rhoIn - (rhoVar)-3.0f * U_MAX * (rhoVar)-3.0f * U_MAX * U_MAX * (rhoVar)) / (9.0f * (rhoVar));
        *myy = 0.0f;

        *rho = rhoVar;

        break;
    }
    case INT_LEFT:
    {
        // const dfloat rhoIn = pop[0] + pop[2] + pop[2] + pop[3] + pop[4] + pop[6] + pop[7];

        // const dfloat rhoUxIn = -(pop[3] + pop[6] + pop[7]);
        // const dfloat rhoUyIn = (pop[2] + pop[6]) - (pop[4] + pop[7]);

        // const dfloat rhoMxxIn = (pop[3] + pop[6] + pop[7]) - rhoIn * cs2;
        // const dfloat rhoMxyIn = pop[7] - pop[6];
        // const dfloat rhoMyyIn = (pop[2] + pop[4] + pop[6] + pop[7]) - rhoIn * cs2;

        // dfloat rho, rhoMxx;

        // const dfloat rhoUy = static_cast<dfloat>(1.5) * (rhoMxyIn + rhoUyIn);
        // const dfloat rhoMxy = static_cast<dfloat>(0.5) * (static_cast<dfloat>(5) * rhoMxyIn + rhoUyIn);
        // const dfloat rhoMyy = static_cast<dfloat>(1.2) * rhoMyyIn;

        // // newton_raphson(rhoIn, rhoUxIn, omega, &rho, rhoUy, &rhoMxx, rhoMyy);

        // const dfloat rhoUx = (static_cast<dfloat>(6) * rhoUxIn + rho + static_cast<dfloat>(3) * rhoMxx) / static_cast<dfloat>(3);

        // *rho = rho;
        // *ux = rhoUx / rho;
        // *uy = rhoUy / rho;
        // *mxx = rhoMxx / rho;
        // *mxy = rhoMxy / rho;
        // *myy = rhoMyy / rho;

        break;
    }
    case INT_TOP:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[3] + pop[5] + pop[6];

        const dfloat rhoUxIn = (pop[1] + pop[5]) - (pop[3] + pop[6]);
        const dfloat rhoUyIn = pop[2] + pop[5] + pop[6];

        const dfloat rhoMxxIn = (pop[1] + pop[3] + pop[5] + pop[6]) - rhoIn * cs2;
        const dfloat rhoMxyIn = (pop[5] - pop[6]);
        const dfloat rhoMyyIn = (pop[2] + pop[5] + pop[6]) - rhoIn * cs2;

        dfloat rhoVar, rhoMyy;

        const dfloat rhoUx = -static_cast<dfloat>(1.5) * (rhoMxyIn - rhoUxIn);
        const dfloat rhoMxx = static_cast<dfloat>(1.2) * rhoMxxIn;
        const dfloat rhoMxy = static_cast<dfloat>(0.5) * (static_cast<dfloat>(5) * rhoMxyIn - rhoUxIn);

        newton_raphson(rhoIn, rhoUyIn, omega, &rhoVar, rhoUx, rhoMxx, &rhoMyy);

        const dfloat rhoUy = (static_cast<dfloat>(6) * rhoUyIn - rhoVar - static_cast<dfloat>(3) * rhoMyy) / static_cast<dfloat>(3);

        *rho = rhoVar;
        *ux = rhoUx / rhoVar;
        *uy = rhoUy / rhoVar;
        *mxx = rhoMxx / rhoVar;
        *mxy = rhoMxy / rhoVar;
        *myy = rhoMyy / rhoVar;
    }
    case INT_TOP_RIGHT:
    {
    }
    case INT_TOP_LEFT:
    {
        const dfloat rhoIn = pop[0] + pop[2] + pop[3] + pop[6];
        const dfloat inv_rho_In = 1.0f / rhoIn;

        const dfloat mxyIn = -pop[6] * inv_rho_In;

        const dfloat rhoVar = -36.0f * (-rhoIn - mxyIn * rhoIn + rhoIn * mxyIn * omega) /
                              (24.0f + 18.0f * U_MAX - 18.0f * U_MAX * U_MAX + omega - 3.0f * U_MAX * omega + 3.0f * U_MAX * U_MAX * omega);
        const dfloat mxyVar = (36.0f * mxyIn * rhoIn + rhoVar - 3.0f * U_MAX * rhoVar + 3.0f * U_MAX * U_MAX * rhoVar) /
                              (9.0f * rhoVar);

        *rho = rhoVar;
        *ux = U_MAX;
        *uy = 0.0f;
        *mxx = U_MAX * U_MAX;
        *mxy = mxyVar;
        *myy = 0.0f;

        break;
    }
    case INT_BOTTOM_LEFT:
    {
        const dfloat rhoIn = pop[0] + pop[3] + pop[4] + pop[7];
        const dfloat inv_rho_In = 1.0f / rhoIn;

        const dfloat mxyIn = pop[7] * inv_rho_In;

        const dfloat rhoVar = 36.0f * (rhoIn - mxyIn * rhoIn + rhoIn * mxyIn * omega) /
                              (24.0f + omega);
        const dfloat mxyVar = (36.0f * mxyIn * rhoIn - rhoVar) /
                              (9.0f * rhoVar);

        *rho = rhoVar;
        *ux = 0.0f;
        *uy = 0.0f;
        *mxx = 0.0f;
        *mxy = mxyVar;
        *myy = 0.0f;

        break;
    }
    default:
        break;
    }
}

__host__ inline void streaming_fine(dfloat *&pop_in, dfloat *&pop_out)
{
    for (size_t y = 0; y < NY_FINE; ++y)
    {
        for (size_t x = 0; x < NX_FINE; ++x)
        {
            size_t xp1 = (x + 1 + NX_FINE) % NX_FINE;
            size_t xm1 = (x - 1 + NX_FINE) % NX_FINE;
            size_t yp1 = (y + 1 + NY_FINE) % NY_FINE;
            size_t ym1 = (y - 1 + NY_FINE) % NY_FINE;

            pop_in[fine_pop_idx(x, y, 0)] = pop_out[fine_pop_idx(x, y, 0)];
            pop_in[fine_pop_idx(xp1, y, 1)] = pop_out[fine_pop_idx(x, y, 1)];
            pop_in[fine_pop_idx(x, yp1, 2)] = pop_out[fine_pop_idx(x, y, 2)];
            pop_in[fine_pop_idx(xm1, y, 3)] = pop_out[fine_pop_idx(x, y, 3)];
            pop_in[fine_pop_idx(x, ym1, 4)] = pop_out[fine_pop_idx(x, y, 4)];
            pop_in[fine_pop_idx(xp1, yp1, 5)] = pop_out[fine_pop_idx(x, y, 5)];
            pop_in[fine_pop_idx(xm1, yp1, 6)] = pop_out[fine_pop_idx(x, y, 6)];
            pop_in[fine_pop_idx(xm1, ym1, 7)] = pop_out[fine_pop_idx(x, y, 7)];
            pop_in[fine_pop_idx(xp1, ym1, 8)] = pop_out[fine_pop_idx(x, y, 8)];
        }
    }
}

__host__ inline void streaming_coarse(dfloat *&pop_in, dfloat *&pop_out)
{
    for (size_t y = 0; y < NY_COARSE + N_OVERLAP_LAYER; ++y)
    {
        for (size_t x = 0; x < NX_COARSE; ++x)
        {
            size_t xp1 = (x + 1 + NX_COARSE) % NX_COARSE;
            size_t xm1 = (x - 1 + NX_COARSE) % NX_COARSE;
            size_t yp1 = (y + 1 + (NY_COARSE + N_OVERLAP_LAYER)) % (NY_COARSE + N_OVERLAP_LAYER);
            size_t ym1 = (y - 1 + (NY_COARSE + N_OVERLAP_LAYER)) % (NY_COARSE + N_OVERLAP_LAYER);

            pop_in[coarse_pop_idx(x, y, 0)] = pop_out[coarse_pop_idx(x, y, 0)];
            pop_in[coarse_pop_idx(xp1, y, 1)] = pop_out[coarse_pop_idx(x, y, 1)];
            pop_in[coarse_pop_idx(x, yp1, 2)] = pop_out[coarse_pop_idx(x, y, 2)];
            pop_in[coarse_pop_idx(xm1, y, 3)] = pop_out[coarse_pop_idx(x, y, 3)];
            pop_in[coarse_pop_idx(x, ym1, 4)] = pop_out[coarse_pop_idx(x, y, 4)];
            pop_in[coarse_pop_idx(xp1, yp1, 5)] = pop_out[coarse_pop_idx(x, y, 5)];
            pop_in[coarse_pop_idx(xm1, yp1, 6)] = pop_out[coarse_pop_idx(x, y, 6)];
            pop_in[coarse_pop_idx(xm1, ym1, 7)] = pop_out[coarse_pop_idx(x, y, 7)];
            pop_in[coarse_pop_idx(xp1, ym1, 8)] = pop_out[coarse_pop_idx(x, y, 8)];
        }
    }
}

__host__ inline void collision(dfloat *mxx, dfloat *mxy, dfloat *myy, dfloat ux, dfloat uy, dfloat omega)
{
    const dfloat omegaVar = omega;
    const dfloat t_omegaVar = 1 - omegaVar;
    const dfloat omegaVar_d2 = omegaVar / 2;

    *mxx = (t_omegaVar * (*mxx) + omegaVar_d2 * ux * ux);
    *myy = (t_omegaVar * (*myy) + omegaVar_d2 * uy * uy);

    *mxy = (t_omegaVar * (*mxy) + omegaVar * ux * uy);
}

__host__ inline void coarse_to_fine(latticeNode *coarse_nodes, latticeNode *fine_nodes)
{
    for (int x = 0; x < NX_COARSE; ++x)
    {
        const size_t fine_id = fine_idx(x * 2, NY_FINE - 1);
        const size_t coarse_id = coarse_idx(x, 0);

        const unsigned temp_type = fine_nodes[fine_id].node_type;

        fine_nodes[fine_id] = coarse_nodes[coarse_id];

        fine_nodes[fine_id].node_type = temp_type;
    }
}

__host__ inline void fine_to_coarse(latticeNode *fine_nodes, latticeNode *coarse_nodes)
{
    for (int x = 0; x < NX_COARSE; ++x)
    {
        const size_t fine_id = fine_idx(x * 2, (NY_FINE - 1) - 2);
        const size_t coarse_id = coarse_idx(x, 0);

        coarse_nodes[coarse_id] = fine_nodes[fine_id];
    }
}

#endif