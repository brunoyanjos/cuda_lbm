#ifndef LBM_STEPS_CUH
#define LBM_STEPS_CUH

#include "var.h"
#include "nodeTypeMap.h"
#include "globalFunctions.h"

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

__host__ inline void pop_streaming(dfloat *&pop_in, dfloat *&pop_out, size_t nx, size_t ny)
{
    for (size_t y = 0; y < ny; ++y)
    {
        for (size_t x = 0; x < nx; ++x)
        {
            size_t xp1 = (x + 1 + nx) % nx;
            size_t xm1 = (x - 1 + nx) % nx;
            size_t yp1 = (y + 1 + ny) % ny;
            size_t ym1 = (y - 1 + ny) % ny;

            pop_in[idx_pop(x, y, 0, nx)] = pop_out[idx_pop(x, y, 0, nx)];
            pop_in[idx_pop(xp1, y, 1, nx)] = pop_out[idx_pop(x, y, 1, nx)];
            pop_in[idx_pop(x, yp1, 2, nx)] = pop_out[idx_pop(x, y, 2, nx)];
            pop_in[idx_pop(xm1, y, 3, nx)] = pop_out[idx_pop(x, y, 3, nx)];
            pop_in[idx_pop(x, ym1, 4, nx)] = pop_out[idx_pop(x, y, 4, nx)];
            pop_in[idx_pop(xp1, yp1, 5, nx)] = pop_out[idx_pop(x, y, 5, nx)];
            pop_in[idx_pop(xm1, yp1, 6, nx)] = pop_out[idx_pop(x, y, 6, nx)];
            pop_in[idx_pop(xm1, ym1, 7, nx)] = pop_out[idx_pop(x, y, 7, nx)];
            pop_in[idx_pop(xp1, ym1, 8, nx)] = pop_out[idx_pop(x, y, 8, nx)];
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

__host__ inline void coarse_to_fine(dfloat *moments_coarse, dfloat *moments_fine, unsigned int *node_type_coarse, unsigned int *node_type_fine)
{
    for (int x = 0; x < NX_COARSE; ++x)
    {
    }
}

__host__ inline void fine_to_coarse(dfloat *moments_fine, dfloat *moments_coarse, unsigned int *node_type_fine, unsigned int *node_type_coarse)
{
    for (int x = 0; x < NX_COARSE; ++x)
    {
    }
}

#endif