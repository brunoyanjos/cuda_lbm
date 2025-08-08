#ifndef METHOD_SOLUTIONS_CUH
#define METHOD_SOLUTIONS_CUH

#include "lbm_steps.cuh"
#include "var.h"
#include "globalStructs.h"
#include "globalFunctions.h"
#include "nodeTypeMap.h"

__host__ inline void fine_grid_solution(latticeNode *nodes, bool isEven)
{
    for (size_t y = 0; y < NY_FINE; ++y)
    {
        for (size_t x = 0; x < NX_FINE; ++x)
        {
            collision(&nodes[fine_idx(x, y)], OMEGA_FINE);
            regularization(&nodes[fine_idx(x, y)]);
        }
    }

    streaming_fine(nodes, NX_FINE, NY_FINE);

    for (size_t y = 0; y < NY_FINE; ++y)
    {
        for (size_t x = 0; x < NX_FINE; ++x)
        {
            unsigned int nodeType = nodes[fine_idx(x, y)].node_type;

            if (nodeType != BULK)
            {
                boundary_condition(&nodes[fine_idx(x, y)], OMEGA_FINE);
            }
            else
            {
                const dfloat *pop = nodes[fine_idx(x, y)].pop_in;

                nodes[fine_idx(x, y)].rho = pop[0] + pop[1] + pop[2] + pop[3] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8];
                const dfloat invRho = 1.0f / nodes[fine_idx(x, y)].rho;

                nodes[fine_idx(x, y)].ux = ((pop[1] + pop[5] + pop[8]) - (pop[3] + pop[6] + pop[7])) * invRho;
                nodes[fine_idx(x, y)].uy = ((pop[2] + pop[5] + pop[6]) - (pop[4] + pop[7] + pop[8])) * invRho;

                nodes[fine_idx(x, y)].mxx = (pop[1] + pop[3] + pop[5] + pop[6] + pop[7] + pop[8]) * invRho - cs2;
                nodes[fine_idx(x, y)].mxy = ((pop[5] + pop[7]) - (pop[6] + pop[8])) * invRho;
                nodes[fine_idx(x, y)].myy = (pop[2] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8]) * invRho - cs2;
            }

            nodes[fine_idx(x, y)].ux = nodes[fine_idx(x, y)].ux * F_M_I_SCALE;
            nodes[fine_idx(x, y)].uy = nodes[fine_idx(x, y)].uy * F_M_I_SCALE;
            nodes[fine_idx(x, y)].mxx = nodes[fine_idx(x, y)].mxx * F_M_II_SCALE;
            nodes[fine_idx(x, y)].mxy = nodes[fine_idx(x, y)].mxy * F_M_IJ_SCALE;
            nodes[fine_idx(x, y)].myy = nodes[fine_idx(x, y)].myy * F_M_II_SCALE;
        }
    }
}

__host__ inline void coarse_grid_solution(latticeNode *nodes)
{
    for (size_t y = 0; y < NY_COARSE; y++)
    {
        for (size_t x = 0; x < NX_COARSE + N_OVERLAP_LAYER; x++)
        {
            collision(&nodes[coarse_idx(x, y)], OMEGA_COARSE);
            regularization(&nodes[coarse_idx(x, y)]);
        }
    }

    streaming_coarse(nodes, NX_COARSE + N_OVERLAP_LAYER, NY_COARSE);

    for (size_t y = 0; y < NY_COARSE; y++)
    {
        for (size_t x = 0; x < NX_COARSE + N_OVERLAP_LAYER; x++)
        {
            unsigned int nodeType = nodes[coarse_idx(x, y)].node_type;

            if (nodeType != BULK)
            {
                boundary_condition(&nodes[coarse_idx(x, y)], OMEGA_COARSE);
            }
            else
            {
                const dfloat *pop = nodes[coarse_idx(x, y)].pop_in;

                nodes[coarse_idx(x, y)].rho = pop[0] + pop[1] + pop[2] + pop[3] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8];
                const dfloat invRho = 1.0f / nodes[coarse_idx(x, y)].rho;

                nodes[coarse_idx(x, y)].ux = ((pop[1] + pop[5] + pop[8]) - (pop[3] + pop[6] + pop[7])) * invRho;
                nodes[coarse_idx(x, y)].uy = ((pop[2] + pop[5] + pop[6]) - (pop[4] + pop[7] + pop[8])) * invRho;

                nodes[coarse_idx(x, y)].mxx = (pop[1] + pop[3] + pop[5] + pop[6] + pop[7] + pop[8]) * invRho - cs2;
                nodes[coarse_idx(x, y)].mxy = ((pop[5] + pop[7]) - (pop[6] + pop[8])) * invRho;
                nodes[coarse_idx(x, y)].myy = (pop[2] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8]) * invRho - cs2;
            }

            nodes[coarse_idx(x, y)].ux = nodes[coarse_idx(x, y)].ux * F_M_I_SCALE;
            nodes[coarse_idx(x, y)].uy = nodes[coarse_idx(x, y)].uy * F_M_I_SCALE;
            nodes[coarse_idx(x, y)].mxx = nodes[coarse_idx(x, y)].mxx * F_M_II_SCALE;
            nodes[coarse_idx(x, y)].mxy = nodes[coarse_idx(x, y)].mxy * F_M_IJ_SCALE;
            nodes[coarse_idx(x, y)].myy = nodes[coarse_idx(x, y)].myy * F_M_II_SCALE;
        }
    }
}

#endif