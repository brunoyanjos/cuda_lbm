#ifndef METHOD_SOLUTIONS_CUH
#define METHOD_SOLUTIONS_CUH

#include "lbm_steps.cuh"
#include "var.h"
#include "globalFunctions.h"
#include "nodeTypeMap.h"
#include CASE_BC
#include EVAL_BC

__host__ inline void grid_solution(unsigned int *node_type, dfloat *moments, dfloat *pop_in, dfloat *pop_out,
                                   dfloat omega, size_t nx, size_t ny)
{
    for (size_t y = 0; y < ny; ++y)
    {
        for (size_t x = 0; x < nx; ++x)
        {
            const dfloat rho = moments[idx_mom(x, y, M_RHO_INDEX, nx)] + RHO_0;
            const dfloat ux = moments[idx_mom(x, y, M_UX_INDEX, nx)];
            const dfloat uy = moments[idx_mom(x, y, M_UY_INDEX, nx)];
            dfloat mxx = moments[idx_mom(x, y, M_MXX_INDEX, nx)];
            dfloat mxy = moments[idx_mom(x, y, M_MXY_INDEX, nx)];
            dfloat myy = moments[idx_mom(x, y, M_MYY_INDEX, nx)];

            dfloat *pop = pop_out + idx_grid(x, y, nx) * Q;

            regularization(pop, rho, ux, uy, mxx, mxy, myy);
        }
    }

    pop_streaming(pop_in, pop_out, nx, ny);

    for (size_t y = 0; y < ny; ++y)
    {
        for (size_t x = 0; x < nx; ++x)
        {
            unsigned int nodeType = node_type[idx_grid(x, y, nx)];
            const dfloat *pop = pop_in + idx_grid(x, y, nx) * Q;

            dfloat rho, ux, uy, mxx, mxy, myy;

            if (nodeType != BULK)
            {
                evaluate_boundary(nodeType, &rho, &ux, &uy,
                                  &mxx, &mxy, &myy, pop, moments,
                                  x, y, omega, nx);
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

            collision(&mxx, &mxy, &myy, ux, uy, omega);

            moments[idx_mom(x, y, M_RHO_INDEX, nx)] = rho - RHO_0;
            moments[idx_mom(x, y, M_UX_INDEX, nx)] = ux;
            moments[idx_mom(x, y, M_UY_INDEX, nx)] = uy;
            moments[idx_mom(x, y, M_MXX_INDEX, nx)] = mxx;
            moments[idx_mom(x, y, M_MXY_INDEX, nx)] = mxy;
            moments[idx_mom(x, y, M_MYY_INDEX, nx)] = myy;
        }
    }
}

#endif