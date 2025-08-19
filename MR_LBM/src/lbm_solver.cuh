#ifndef METHOD_SOLUTIONS_CUH
#define METHOD_SOLUTIONS_CUH

#include "lbm_steps.cuh"
#include "var.h"
#include "globalStructs.h"
#include "globalFunctions.h"
#include "nodeTypeMap.h"

__host__ inline void fine_grid_solution(unsigned int *node_type, dfloat *moments, dfloat *pop_in, dfloat *pop_out)
{
    for (size_t y = 0; y < NY_FINE; ++y)
    {
        for (size_t x = 0; x < NX_FINE; ++x)
        {
            dfloat rho = moments[fine_moment_idx(x, y, M_RHO_INDEX)];
            dfloat ux = moments[fine_moment_idx(x, y, M_UX_INDEX)];
            dfloat uy = moments[fine_moment_idx(x, y, M_UY_INDEX)];
            dfloat mxx = moments[fine_moment_idx(x, y, M_MXX_INDEX)];
            dfloat mxy = moments[fine_moment_idx(x, y, M_MXY_INDEX)];
            dfloat myy = moments[fine_moment_idx(x, y, M_MYY_INDEX)];

            collision(&mxx, &mxy, &myy, ux, uy, OMEGA_FINE);
            regularization(pop_out + fine_idx(x, y) * Q, rho, ux, uy, mxx, mxy, myy);
        }
    }

    streaming_fine(pop_in, pop_out);

    for (size_t y = 0; y < NY_FINE; ++y)
    {
        for (size_t x = 0; x < NX_FINE; ++x)
        {
            unsigned int nodeType = node_type[fine_idx(x, y)];

            dfloat rho, ux, uy, mxx, mxy, myy;

            const dfloat *pop = pop_in + fine_idx(x, y) * Q;

            if (nodeType != BULK)
            {
                boundary_condition(nodeType, pop, &rho, &ux, &uy,
                                   &mxx, &mxy, &myy, OMEGA_FINE);
            }
            else
            {
                rho = pop[0] + pop[1] + pop[2] + pop[3] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8];
                const dfloat invRho = 1.0f / rho;

                ux = ((pop[1] + pop[5] + pop[8]) - (pop[3] + pop[6] + pop[7])) * invRho;
                uy = ((pop[2] + pop[5] + pop[6]) - (pop[4] + pop[7] + pop[8])) * invRho;

                mxx = (pop[1] + pop[3] + pop[5] + pop[6] + pop[7] + pop[8]) * invRho - cs2;
                mxy = ((pop[5] + pop[7]) - (pop[6] + pop[8])) * invRho;
                myy = (pop[2] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8]) * invRho - cs2;
            }

            ux = ux * F_M_I_SCALE;
            uy = uy * F_M_I_SCALE;
            mxx = mxx * F_M_II_SCALE;
            mxy = mxy * F_M_IJ_SCALE;
            myy = myy * F_M_II_SCALE;

            moments[fine_moment_idx(x, y, M_RHO_INDEX)] = rho;
            moments[fine_moment_idx(x, y, M_UX_INDEX)] = ux;
            moments[fine_moment_idx(x, y, M_UY_INDEX)] = uy;
            moments[fine_moment_idx(x, y, M_MXX_INDEX)] = mxx;
            moments[fine_moment_idx(x, y, M_MXY_INDEX)] = mxy;
            moments[fine_moment_idx(x, y, M_MYY_INDEX)] = myy;
        }
    }
}

__host__ inline void coarse_grid_solution(unsigned int *node_type, dfloat *moments, dfloat *pop_in, dfloat *pop_out)
{
    for (size_t y = 0; y < NY_COARSE + N_OVERLAP_LAYER; y++)
    {
        for (size_t x = 0; x < NX_COARSE; x++)
        {
            dfloat rho = moments[coarse_moment_idx(x, y, M_RHO_INDEX)];
            dfloat ux = moments[coarse_moment_idx(x, y, M_UX_INDEX)];
            dfloat uy = moments[coarse_moment_idx(x, y, M_UY_INDEX)];
            dfloat mxx = moments[coarse_moment_idx(x, y, M_MXX_INDEX)];
            dfloat mxy = moments[coarse_moment_idx(x, y, M_MXY_INDEX)];
            dfloat myy = moments[coarse_moment_idx(x, y, M_MYY_INDEX)];

            collision(&mxx, &mxy, &myy, ux, uy, OMEGA_COARSE);
            regularization(pop_out + coarse_idx(x, y) * Q, rho, ux, uy, mxx, mxy, myy);
        }
    }

    streaming_coarse(pop_in, pop_out);

    for (size_t y = 0; y < NY_COARSE + N_OVERLAP_LAYER; y++)
    {
        for (size_t x = 0; x < NX_COARSE; x++)
        {
            unsigned int nodeType = node_type[coarse_idx(x, y)];

            dfloat rho, ux, uy, mxx, mxy, myy;

            const dfloat *pop = pop_in + coarse_idx(x, y) * Q;

            if (nodeType != BULK)
            {
                boundary_condition(nodeType, pop, &rho, &ux, &uy,
                                   &mxx, &mxy, &myy, OMEGA_COARSE);
            }
            else
            {
                rho = pop[0] + pop[1] + pop[2] + pop[3] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8];
                const dfloat invRho = 1.0f / rho;

                ux = ((pop[1] + pop[5] + pop[8]) - (pop[3] + pop[6] + pop[7])) * invRho;
                uy = ((pop[2] + pop[5] + pop[6]) - (pop[4] + pop[7] + pop[8])) * invRho;

                mxx = (pop[1] + pop[3] + pop[5] + pop[6] + pop[7] + pop[8]) * invRho - cs2;
                mxy = ((pop[5] + pop[7]) - (pop[6] + pop[8])) * invRho;
                myy = (pop[2] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8]) * invRho - cs2;
            }

            ux = ux * F_M_I_SCALE;
            uy = uy * F_M_I_SCALE;
            mxx = mxx * F_M_II_SCALE;
            mxy = mxy * F_M_IJ_SCALE;
            myy = myy * F_M_II_SCALE;

            moments[coarse_moment_idx(x, y, M_RHO_INDEX)] = rho;
            moments[coarse_moment_idx(x, y, M_UX_INDEX)] = ux;
            moments[coarse_moment_idx(x, y, M_UY_INDEX)] = uy;
            moments[coarse_moment_idx(x, y, M_MXX_INDEX)] = mxx;
            moments[coarse_moment_idx(x, y, M_MXY_INDEX)] = mxy;
            moments[coarse_moment_idx(x, y, M_MYY_INDEX)] = myy;
        }
    }
}

#endif