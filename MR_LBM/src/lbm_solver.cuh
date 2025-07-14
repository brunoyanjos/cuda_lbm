#ifndef METHOD_SOLUTIONS_CUH
#define METHOD_SOLUTIONS_CUH

#include "lbm_steps.cuh"
#include "var.h"
#include "globalStructs.h"
#include "globalFunctions.h"
#include "nodeTypeMap.h"

__host__ inline void fine_grid_solution(latticeNode *nodes)
{
    for (size_t y = 0; y < NY_FINE_GRID; ++y)
    {
        for (size_t x = 0; x < NX_FINE_GRID; ++x)
        {
            unsigned int nodeType = nodes[fine_idx(x, y)].node_type;

            if (nodeType != MISSING_DEFINITION)
            {
                if (nodeType != BULK)
                {
                    if (nodes[fine_idx(x, y)].updated)
                    {
                        const dfloat *pop = nodes[fine_idx(x, y)].pop_in;

                        nodes[fine_idx(x, y)].rho = pop[0] + pop[1] + pop[2] + pop[3] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8];
                        const dfloat invRho = 1.0f / nodes[fine_idx(x, y)].rho;

                        nodes[fine_idx(x, y)].ux = ((pop[1] + pop[5] + pop[8]) - (pop[3] + pop[6] + pop[7])) * invRho;
                        nodes[fine_idx(x, y)].uy = ((pop[2] + pop[5] + pop[6]) - (pop[4] + pop[7] + pop[8])) * invRho;

                        nodes[fine_idx(x, y)].mxx = (pop[1] + pop[3] + pop[5] + pop[6] + pop[7] + pop[8]) * invRho - cs2;
                        nodes[fine_idx(x, y)].mxy = ((pop[5] + pop[7]) - (pop[6] + pop[8])) * invRho;
                        nodes[fine_idx(x, y)].myy = (pop[2] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8]) * invRho - cs2;

                        nodes[fine_idx(x, y)].updated = false;
                    }
                    else
                    {
                        boundary_condition(&nodes[fine_idx(x, y)]);
                    }
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
            }

            nodes[fine_idx(x, y)].ux = nodes[fine_idx(x, y)].ux * F_M_I_SCALE;
            nodes[fine_idx(x, y)].uy = nodes[fine_idx(x, y)].uy * F_M_I_SCALE;
            nodes[fine_idx(x, y)].mxx = nodes[fine_idx(x, y)].mxx * F_M_II_SCALE;
            nodes[fine_idx(x, y)].mxy = nodes[fine_idx(x, y)].mxy * F_M_IJ_SCALE;
            nodes[fine_idx(x, y)].myy = nodes[fine_idx(x, y)].myy * F_M_II_SCALE;

            collision(&nodes[fine_idx(x, y)], OMEGA_FINE);
            regularization(&nodes[fine_idx(x, y)]);
        }
    }

    streaming(nodes, NX_FINE_GRID, NY_FINE_GRID);
}

__host__ inline void coarse_grid_solution(latticeNode *nodes)
{
    for (size_t y = 0; y < NY_COARSE_GRID; y++)
    {
        for (size_t x = 0; x < NX_COARSE_GRID; x++)
        {
            const dfloat *pop = nodes[fine_idx(x, y)].pop_in;

            nodes[fine_idx(x, y)].rho = pop[0] + pop[1] + pop[2] + pop[3] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8];
            const dfloat invRho = 1.0f / nodes[fine_idx(x, y)].rho;

            nodes[fine_idx(x, y)].ux = ((pop[1] + pop[5] + pop[8]) - (pop[3] + pop[6] + pop[7])) * invRho;
            nodes[fine_idx(x, y)].uy = ((pop[2] + pop[5] + pop[6]) - (pop[4] + pop[7] + pop[8])) * invRho;

            nodes[fine_idx(x, y)].mxx = (pop[1] + pop[3] + pop[5] + pop[6] + pop[7] + pop[8]) * invRho - cs2;
            nodes[fine_idx(x, y)].mxy = ((pop[5] + pop[7]) - (pop[6] + pop[8])) * invRho;
            nodes[fine_idx(x, y)].myy = (pop[2] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8]) * invRho - cs2;

            nodes[fine_idx(x, y)].ux = nodes[fine_idx(x, y)].ux * F_M_I_SCALE;
            nodes[fine_idx(x, y)].uy = nodes[fine_idx(x, y)].uy * F_M_I_SCALE;
            nodes[fine_idx(x, y)].mxx = nodes[fine_idx(x, y)].mxx * F_M_II_SCALE;
            nodes[fine_idx(x, y)].mxy = nodes[fine_idx(x, y)].mxy * F_M_IJ_SCALE;
            nodes[fine_idx(x, y)].myy = nodes[fine_idx(x, y)].myy * F_M_II_SCALE;

            collision(&nodes[fine_idx(x, y)], OMEGA_FINE);
            regularization(&nodes[fine_idx(x, y)]);
        }
    }

    streaming(nodes, NX_COARSE_GRID, NY_COARSE_GRID);
}

#endif