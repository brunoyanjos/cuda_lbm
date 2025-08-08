#ifndef LBM_STEPS_CUH
#define LBM_STEPS_CUH

#include "var.h"

__host__ inline void init_pop_eq(latticeNode *node)
{
    dfloat rho = (*node).rho;
    dfloat ux = (*node).ux * F_M_I_SCALE;
    dfloat uy = (*node).uy * F_M_I_SCALE;

    dfloat pics2 = 1 - cs2 * (ux * ux + uy * uy);

    dfloat multiplyTerm = W0 * rho;
    (*node).pop_in[0] = multiplyTerm * (pics2);

    multiplyTerm = W1 * rho;
    (*node).pop_in[1] = multiplyTerm * (pics2 + ux + static_cast<dfloat>(0.5) * ux * ux);
    (*node).pop_in[2] = multiplyTerm * (pics2 + uy + static_cast<dfloat>(0.5) * uy * uy);
    (*node).pop_in[3] = multiplyTerm * (pics2 - ux + static_cast<dfloat>(0.5) * ux * ux);
    (*node).pop_in[4] = multiplyTerm * (pics2 - uy + static_cast<dfloat>(0.5) * uy * uy);

    multiplyTerm = W2 * rho;
    (*node).pop_in[5] = multiplyTerm * (pics2 + ux + uy + static_cast<dfloat>(0.5) * ux * ux + static_cast<dfloat>(0.5) * uy * uy + ux * uy);
    (*node).pop_in[6] = multiplyTerm * (pics2 - ux + uy + static_cast<dfloat>(0.5) * ux * ux + static_cast<dfloat>(0.5) * uy * uy - ux * uy);
    (*node).pop_in[7] = multiplyTerm * (pics2 - ux - uy + static_cast<dfloat>(0.5) * ux * ux + static_cast<dfloat>(0.5) * uy * uy + ux * uy);
    (*node).pop_in[8] = multiplyTerm * (pics2 + ux - uy + static_cast<dfloat>(0.5) * ux * ux + static_cast<dfloat>(0.5) * uy * uy - ux * uy);
}

__host__ inline void regularization(latticeNode *node)
{
    dfloat rho = (*node).rho;
    dfloat ux = (*node).ux;
    dfloat uy = (*node).uy;
    dfloat mxx = (*node).mxx;
    dfloat mxy = (*node).mxy;
    dfloat myy = (*node).myy;

    dfloat pics2 = 1 - cs2 * (mxx + myy);

    dfloat multiplyTerm = W0 * rho;
    (*node).pop_out[0] = multiplyTerm * (pics2);

    multiplyTerm = W1 * rho;
    (*node).pop_out[1] = multiplyTerm * (pics2 + ux + mxx);
    (*node).pop_out[2] = multiplyTerm * (pics2 + uy + myy);
    (*node).pop_out[3] = multiplyTerm * (pics2 - ux + mxx);
    (*node).pop_out[4] = multiplyTerm * (pics2 - uy + myy);

    multiplyTerm = W2 * rho;
    (*node).pop_out[5] = multiplyTerm * (pics2 + ux + uy + mxx + myy + mxy);
    (*node).pop_out[6] = multiplyTerm * (pics2 - ux + uy + mxx + myy - mxy);
    (*node).pop_out[7] = multiplyTerm * (pics2 - ux - uy + mxx + myy + mxy);
    (*node).pop_out[8] = multiplyTerm * (pics2 + ux - uy + mxx + myy - mxy);
}

__host__ inline void boundary_condition(latticeNode *node, dfloat omega)
{
    unsigned int nodeType = (*node).node_type;
    const dfloat *pop = (*node).pop_in;

    switch (nodeType)
    {
    case NORTH:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[3] + pop[5] + pop[6];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat mxxIn = (pop[1] + pop[3] + pop[5] + pop[6]) * inv_rhoIn - cs2;
        const dfloat mxyIn = (pop[5] - pop[6]) * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[5] + pop[6]) * inv_rhoIn - cs2;

        (*node).ux = U_MAX;
        (*node).uy = 0.0f;

        const dfloat rhoVar = 6.0f * rhoIn / 5.0f;

        (*node).mxx = U_MAX * U_MAX;
        (*node).mxx = 0.0f;
        (*node).mxy = 5.0f * mxyIn / 3.0f - U_MAX / 3.0f;
        (*node).myy = 0.0f;

        (*node).rho = rhoVar;

        break;
    }
    case SOUTH:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[3] + pop[4] + pop[7] + pop[8];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat mxxIn = (pop[1] + pop[3] + pop[7] + pop[8]) * inv_rhoIn - cs2;
        const dfloat mxyIn = (pop[7] - pop[8]) * inv_rhoIn;
        const dfloat myyIn = (pop[3] + pop[7] + pop[8]) * inv_rhoIn - cs2;

        (*node).ux = 0.0f;
        (*node).uy = 0.0f;

        const dfloat rhoVar = 6.0f * rhoIn / 5.0f;

        (*node).mxx = 0.0f;
        (*node).mxy = 5.0f * mxyIn / 3.0f;
        (*node).myy = 0.0f;

        (*node).rho = rhoVar;

        break;
    }
    case WEST:
    {
        const dfloat rhoIn = pop[0] + pop[2] + pop[3] + pop[4] + pop[6] + pop[7];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat mxxIn = (pop[3] + pop[6] + pop[7]) * inv_rhoIn - cs2;
        const dfloat mxyIn = (pop[7] - pop[6]) * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[4] + pop[6] + pop[7]) * inv_rhoIn - cs2;

        (*node).ux = 0.0f;
        (*node).uy = 0.0f;

        const dfloat rhoVar = 6.0f * rhoIn / 5.0f;

        (*node).mxx = 0.0f;
        (*node).mxy = 5.0f * mxyIn / 3.0f;
        (*node).myy = 0.0f;

        (*node).rho = rhoVar;

        break;
    }
    case EAST:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[4] + pop[5] + pop[8];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat mxxIn = (pop[1] + pop[5] + pop[8]) * inv_rhoIn - cs2;
        const dfloat mxyIn = (pop[5] - pop[8]) * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[4] + pop[5] + pop[8]) * inv_rhoIn - cs2;

        (*node).ux = 0.0f;
        (*node).uy = 0.0f;

        const dfloat rhoVar = 6.0f * rhoIn / 5.0f;

        (*node).mxx = 0.0f;
        (*node).mxy = 5.0f * mxyIn / 3.0f;
        (*node).myy = 0.0f;

        (*node).rho = rhoVar;

        break;
    }
    case SOUTH_WEST:
    {
        const dfloat rhoIn = pop[0] + pop[3] + pop[4] + pop[7];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat mxxIn = (pop[3] + pop[7]) * inv_rhoIn - cs2;
        const dfloat mxyIn = pop[7] * inv_rhoIn;
        const dfloat myyIn = (pop[4] + pop[7]) * inv_rhoIn - cs2;

        (*node).ux = 0.0f;
        (*node).uy = 0.0f;

        const dfloat rhoVar = 36.0f * (rhoIn - mxyIn * rhoIn + mxyIn * omega * rhoIn) /
                              (24.0f + omega);

        (*node).mxx = 0.0f;
        (*node).mxy = (36.0f * mxyIn * rhoIn - (rhoVar)) /
                      (9.0f * (rhoVar));
        (*node).myy = 0.0f;

        (*node).rho = rhoVar;

        break;
    }
    case SOUTH_EAST:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[4] + pop[8];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat mxxIn = (pop[1] + pop[8]) * inv_rhoIn - cs2;
        const dfloat mxyIn = -pop[8] * inv_rhoIn;
        const dfloat myyIn = (pop[4] + pop[8]) * inv_rhoIn - cs2;

        (*node).uy = 0.0f;
        (*node).ux = 0.0f;

        const dfloat rhoVar = -36.0f * (mxyIn * omega * rhoIn - rhoIn - mxyIn * rhoIn) /
                              (24 + omega);

        (*node).mxx = 0.0f;
        (*node).mxy = (36.0f * mxyIn * rhoIn + (rhoVar)) / (9.0f * (rhoVar));
        (*node).myy = 0.0f;

        (*node).rho = rhoVar;

        break;
    }
    case NORTH_WEST:
    {
        const dfloat rhoIn = pop[0] + pop[2] + pop[3] + pop[6];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat mxxIn = (pop[3] + pop[6]) * inv_rhoIn - cs2;
        const dfloat mxyIn = -pop[6] * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[6]) * inv_rhoIn - cs2;

        (*node).ux = U_MAX;
        (*node).uy = 0.0f;

        const dfloat rhoVar = -36.0f * (mxyIn * omega * rhoIn - rhoIn - mxyIn * rhoIn) /
                              (24.0f + omega + 18.0f * U_MAX - 3.0f * omega * U_MAX - 18.0f * U_MAX * U_MAX + 3.0f * omega * U_MAX * U_MAX);

        (*node).mxx = U_MAX * U_MAX;
        (*node).mxy = (36.0f * mxyIn * rhoIn + (rhoVar)-3.0f * U_MAX * (rhoVar) + 3.0f * U_MAX * U_MAX * (rhoVar)) / (9.0f * (rhoVar));
        (*node).myy = 0.0f;

        (*node).rho = rhoVar;

        break;
    }
    case NORTH_EAST:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[5];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat mxxIn = (pop[1] + pop[5]) * inv_rhoIn - cs2;
        const dfloat mxyIn = pop[5] * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[5]) * inv_rhoIn - cs2;

        (*node).ux = U_MAX;
        (*node).uy = 0.0f;

        const dfloat rhoVar = 36.0f * (mxyIn * omega * rhoIn + rhoIn - mxyIn * rhoIn) /
                              (24.0f + omega - 18.0f * U_MAX + 3.0f * omega * U_MAX - 18.0f * U_MAX * U_MAX + 3.0f * omega * U_MAX * U_MAX);

        (*node).mxx = U_MAX * U_MAX;
        (*node).mxy = (36.0f * mxyIn * rhoIn - (rhoVar)-3.0f * U_MAX * (rhoVar)-3.0f * U_MAX * U_MAX * (rhoVar)) / (9.0f * (rhoVar));
        (*node).myy = 0.0f;

        (*node).rho = rhoVar;

        break;
    }
    case INT_LEFT:
    {
        const dfloat rhoIn = pop[0] + pop[2] + pop[3] + pop[4] + pop[6] + pop[7];

        const dfloat rhoUxIn = -(pop[3] + pop[6] + pop[7]);
        const dfloat rhoUyIn = (pop[2] + pop[6]) - (pop[4] + pop[7]);

        const dfloat rhoMxxIn = (pop[3] + pop[6] + pop[7]) - rhoIn * cs2;
        const dfloat rhoMxyIn = pop[7] - pop[6];
        const dfloat rhoMyyIn = (pop[2] + pop[4] + pop[6] + pop[7]) - rhoIn * cs2;

        const dfloat sqrt_term =
            -25.0f * rhoIn * rhoIn -
            225.0f * rhoMxyIn * rhoMxyIn +
            180.0f * rhoIn * rhoMyyIn +
            108.0f * rhoMyyIn * rhoMyyIn -
            250.0f * rhoIn * rhoUxIn -
            540.0f * rhoMyyIn * rhoUxIn +
            575.0f * rhoUxIn * rhoUxIn -
            450.0f * rhoMxyIn * rhoUyIn -
            225.0f * rhoUyIn * rhoUyIn;

        const dfloat rho = (3.0f / 10.0f) *
                           (5.0f * rhoIn + 6.0f * rhoMyyIn - 15.0f * rhoUxIn - (1.0f / std::sqrt(3.0f)) * std::sqrt(sqrt_term));

        const dfloat rhoUx = -rhoIn + rhoUxIn + rho;
        const dfloat rhoUy = (3.0f / 2.0f) * (rhoMxyIn + rhoUyIn);

        const dfloat rhoMxx = (1.0f / 3.0f) * (-3.0f * rhoIn - 3.0f * rhoUxIn + 2.0f * rho);
        const dfloat rhoMxy = (1.0f / 2.0f) * (5.0f * rhoMxyIn + rhoUyIn);
        const dfloat rhoMyy = (6.0f / 5.0f) * rhoMyyIn;

        (*node).rho = rho;
        (*node).ux = rhoUx / rho;
        (*node).uy = rhoUy / rho;
        (*node).mxx = rhoMxx / rho;
        (*node).mxy = rhoMxy / rho;
        (*node).myy = rhoMyy / rho;

        break;
    }
    case INT_TOP_LEFT:
    {
        const dfloat rhoIn = pop[0] + pop[2] + pop[3] + pop[6];
        const dfloat inv_rho_In = 1.0f / rhoIn;

        const dfloat mxyIn = -pop[6] * inv_rho_In;

        const dfloat rho = -36.0f * (-rhoIn - mxyIn * rhoIn + rhoIn * mxyIn * omega) /
                           (24.0f + 18.0f * U_MAX - 18.0f * U_MAX * U_MAX + omega - 3.0f * U_MAX * omega + 3.0f * U_MAX * U_MAX * omega);
        const dfloat mxy = (36.0f * mxyIn * rhoIn + rho - 3.0f * U_MAX * rho + 3.0f * U_MAX * U_MAX * rho) /
                           (9.0f * rho);

        (*node).rho = rho;
        (*node).ux = U_MAX;
        (*node).uy = 0.0f;
        (*node).mxx = U_MAX * U_MAX;
        (*node).mxy = mxy;
        (*node).myy = 0.0f;

        break;
    }
    case INT_BOTTOM_LEFT:
    {
        const dfloat rhoIn = pop[0] + pop[3] + pop[4] + pop[7];
        const dfloat inv_rho_In = 1.0f / rhoIn;

        const dfloat mxyIn = pop[7] * inv_rho_In;

        const dfloat rho = 36.0f * (rhoIn - mxyIn * rhoIn + rhoIn * mxyIn * omega) /
                           (24.0f + omega);
        const dfloat mxy = (36.0f * mxyIn * rhoIn - rho) /
                           (9.0f * rho);

        (*node).rho = rho;
        (*node).ux = 0.0f;
        (*node).uy = 0.0f;
        (*node).mxx = 0.0f;
        (*node).mxy = mxy;
        (*node).myy = 0.0f;

        break;
    }
    default:
        break;
    }
}

__host__ inline void streaming_fine(latticeNode *nodes, size_t x_lattices, size_t y_lattices)
{
    for (size_t y = 0; y < y_lattices; ++y)
    {
        for (size_t x = 0; x < x_lattices; ++x)
        {
            size_t xp1 = (x + 1 + x_lattices) % x_lattices;
            size_t xm1 = (x - 1 + x_lattices) % x_lattices;
            size_t yp1 = (y + 1 + y_lattices) % y_lattices;
            size_t ym1 = (y - 1 + y_lattices) % y_lattices;

            nodes[fine_idx(x, y)].pop_in[0] = nodes[fine_idx(x, y)].pop_out[0];
            nodes[fine_idx(xp1, y)].pop_in[1] = nodes[fine_idx(x, y)].pop_out[1];
            nodes[fine_idx(x, yp1)].pop_in[2] = nodes[fine_idx(x, y)].pop_out[2];
            nodes[fine_idx(xm1, y)].pop_in[3] = nodes[fine_idx(x, y)].pop_out[3];
            nodes[fine_idx(x, ym1)].pop_in[4] = nodes[fine_idx(x, y)].pop_out[4];
            nodes[fine_idx(xp1, yp1)].pop_in[5] = nodes[fine_idx(x, y)].pop_out[5];
            nodes[fine_idx(xm1, yp1)].pop_in[6] = nodes[fine_idx(x, y)].pop_out[6];
            nodes[fine_idx(xm1, ym1)].pop_in[7] = nodes[fine_idx(x, y)].pop_out[7];
            nodes[fine_idx(xp1, ym1)].pop_in[8] = nodes[fine_idx(x, y)].pop_out[8];
        }
    }
}

__host__ inline void streaming_coarse(latticeNode *nodes, size_t x_lattices, size_t y_lattices)
{
    for (size_t y = 0; y < y_lattices; ++y)
    {
        for (size_t x = 0; x < x_lattices; ++x)
        {
            size_t xp1 = (x + 1 + x_lattices) % x_lattices;
            size_t xm1 = (x - 1 + x_lattices) % x_lattices;
            size_t yp1 = (y + 1 + y_lattices) % y_lattices;
            size_t ym1 = (y - 1 + y_lattices) % y_lattices;

            nodes[coarse_idx(x, y)].pop_in[0] = nodes[coarse_idx(x, y)].pop_out[0];
            nodes[coarse_idx(xp1, y)].pop_in[1] = nodes[coarse_idx(x, y)].pop_out[1];
            nodes[coarse_idx(x, yp1)].pop_in[2] = nodes[coarse_idx(x, y)].pop_out[2];
            nodes[coarse_idx(xm1, y)].pop_in[3] = nodes[coarse_idx(x, y)].pop_out[3];
            nodes[coarse_idx(x, ym1)].pop_in[4] = nodes[coarse_idx(x, y)].pop_out[4];
            nodes[coarse_idx(xp1, yp1)].pop_in[5] = nodes[coarse_idx(x, y)].pop_out[5];
            nodes[coarse_idx(xm1, yp1)].pop_in[6] = nodes[coarse_idx(x, y)].pop_out[6];
            nodes[coarse_idx(xm1, ym1)].pop_in[7] = nodes[coarse_idx(x, y)].pop_out[7];
            nodes[coarse_idx(xp1, ym1)].pop_in[8] = nodes[coarse_idx(x, y)].pop_out[8];
        }
    }
}

__host__ inline void collision(latticeNode *node, dfloat omega)
{
    const dfloat omegaVar = omega;
    const dfloat t_omegaVar = 1 - omegaVar;
    const dfloat omegaVar_d2 = omegaVar / 2;

    (*node).mxx = (t_omegaVar * (*node).mxx + omegaVar_d2 * (*node).ux * (*node).ux);
    (*node).myy = (t_omegaVar * (*node).myy + omegaVar_d2 * (*node).uy * (*node).uy);

    (*node).mxy = (t_omegaVar * (*node).mxy + omegaVar * (*node).ux * (*node).uy);
}

__host__ inline void coarse_to_fine(latticeNode *coarse_nodes, latticeNode *fine_nodes)
{
    for (int y = 0; y < NY_COARSE; ++y)
    {
        const unsigned temp_type = fine_nodes[fine_idx(0, y * 2)].node_type;

        fine_nodes[fine_idx(0, y * 2)] = coarse_nodes[coarse_idx(NX_COARSE + N_OVERLAP_LAYER - 1, y)];

        fine_nodes[fine_idx(0, y * 2)].node_type = temp_type;
    }
}

__host__ inline void fine_to_coarse(latticeNode *fine_nodes, latticeNode *coarse_nodes)
{
    for (int y = 0; y < NY_COARSE; ++y)
    {
        coarse_nodes[coarse_idx(NX_COARSE + N_OVERLAP_LAYER - 1, y)] = fine_nodes[fine_idx(2, y * 2)];
    }
}

#endif