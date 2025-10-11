#ifndef METHOD_SOLUTIONS_CUH
#define METHOD_SOLUTIONS_CUH

#include "lbm_steps.cuh"
#include "var.h"
#include "globalFunctions.h"
#include "nodeTypeMap.h"
#include CASE_BC
#include EVAL_BC

__host__ inline void fine_grid_solution(unsigned int *node_type, dfloat *moments, dfloat *pop_in, dfloat *pop_out)
{
    for (size_t y = 0; y < NY_FINE; ++y)
    {
        for (size_t x = 0; x < NX_FINE; ++x)
        {
            dfloat rho = moments[idx_mom(x, y, M_RHO_INDEX, NX_FINE)];
            dfloat ux = moments[idx_mom(x, y, M_UX_INDEX, NX_FINE)];
            dfloat uy = moments[idx_mom(x, y, M_UY_INDEX, NX_FINE)];
            dfloat mxx = moments[idx_mom(x, y, M_MXX_INDEX, NX_FINE)];
            dfloat mxy = moments[idx_mom(x, y, M_MXY_INDEX, NX_FINE)];
            dfloat myy = moments[idx_mom(x, y, M_MYY_INDEX, NX_FINE)];

            collision(&mxx, &mxy, &myy, ux, uy, OMEGA_FINE);
            regularization(pop_out + idx_grid(x, y, NX_FINE) * Q, rho, ux, uy, mxx, mxy, myy);
        }
    }

    pop_streaming(pop_in, pop_out, NX_FINE, NY_FINE);

    for (size_t y = 0; y < NY_FINE; ++y)
    {
        for (size_t x = 0; x < NX_FINE; ++x)
        {
            unsigned int nodeType = node_type[idx_grid(x, y, NX_FINE)];

            dfloat rho, ux, uy, mxx, mxy, myy;

            const dfloat *pop = pop_in + idx_grid(x, y, NX_FINE) * Q;

            if (nodeType != BULK)
            {
                evaluate_boundary(nodeType, &rho, &ux, &uy,
                                  &mxx, &mxy, &myy, pop, moments,
                                  x, y, OMEGA_FINE, NX_FINE);
            }
            else
            {
                rho = pop[0] + pop[1] + pop[2] + pop[3] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8];
                const dfloat invRho = static_cast<dfloat>(1) / rho;

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

            moments[idx_mom(x, y, M_RHO_INDEX, NX_FINE)] = rho;
            moments[idx_mom(x, y, M_UX_INDEX, NX_FINE)] = ux;
            moments[idx_mom(x, y, M_UY_INDEX, NX_FINE)] = uy;
            moments[idx_mom(x, y, M_MXX_INDEX, NX_FINE)] = mxx;
            moments[idx_mom(x, y, M_MXY_INDEX, NX_FINE)] = mxy;
            moments[idx_mom(x, y, M_MYY_INDEX, NX_FINE)] = myy;
        }
    }
}

__host__ inline void coarse_grid_solution(unsigned int *node_type, dfloat *moments, dfloat *pop_in, dfloat *pop_out)
{
    for (size_t y = 0; y < NY_COARSE; y++)
    {
        for (size_t x = 0; x < NX_COARSE; x++)
        {
            dfloat rho = moments[idx_mom(x, y, M_RHO_INDEX, NX_COARSE)];
            dfloat ux = moments[idx_mom(x, y, M_UX_INDEX, NX_COARSE)];
            dfloat uy = moments[idx_mom(x, y, M_UY_INDEX, NX_COARSE)];
            dfloat mxx = moments[idx_mom(x, y, M_MXX_INDEX, NX_COARSE)];
            dfloat mxy = moments[idx_mom(x, y, M_MXY_INDEX, NX_COARSE)];
            dfloat myy = moments[idx_mom(x, y, M_MYY_INDEX, NX_COARSE)];

            collision(&mxx, &mxy, &myy, ux, uy, OMEGA_COARSE);
            regularization(pop_out + idx_grid(x, y, NX_COARSE) * Q, rho, ux, uy, mxx, mxy, myy);
        }
    }

    pop_streaming(pop_in, pop_out, NX_COARSE, NY_COARSE);

    for (size_t y = 0; y < NY_COARSE; y++)
    {
        for (size_t x = 0; x < NX_COARSE; x++)
        {
            unsigned int nodeType = node_type[idx_grid(x, y, NX_COARSE)];

            // printf("%02d ", nodeType);

            dfloat rho, ux, uy, mxx, mxy, myy;

            const dfloat *pop = pop_in + idx_grid(x, y, NX_COARSE) * Q;

            if (nodeType != BULK)
            {
                evaluate_boundary(nodeType, &rho, &ux, &uy,
                                  &mxx, &mxy, &myy, pop, moments,
                                  x, y, OMEGA_COARSE, NX_COARSE);
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

            moments[idx_mom(x, y, M_RHO_INDEX, NX_COARSE)] = rho;
            moments[idx_mom(x, y, M_UX_INDEX, NX_COARSE)] = ux;
            moments[idx_mom(x, y, M_UY_INDEX, NX_COARSE)] = uy;
            moments[idx_mom(x, y, M_MXX_INDEX, NX_COARSE)] = mxx;
            moments[idx_mom(x, y, M_MXY_INDEX, NX_COARSE)] = mxy;
            moments[idx_mom(x, y, M_MYY_INDEX, NX_COARSE)] = myy;
        }
    }
}

#endif