#ifndef LBM_STEPS_CUH
#define LBM_STEPS_CUH

#include "var.h"
#include "nodeTypeMap.h"
#include "globalFunctions.h"
#include "interpolation_utilities.cuh"

__host__ inline void eval_pop_eq(dfloat *pop_in,
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
    for (int y = 0; y < NY_FINE; ++y)
    {
        const int xf = 0;
        const int xc = NX_COARSE - 2;

        dfloat pop[9], pop_eq[9], pop_neq[9];

        if (y % 2 == 0)
        {
            const int yc = y / 2;

            const dfloat rho = moments_coarse[idx_mom(xc, yc, M_RHO_INDEX, NX_COARSE)] + RHO_0;

            const dfloat ux = moments_coarse[idx_mom(xc, yc, M_UX_INDEX, NX_COARSE)];
            const dfloat uy = moments_coarse[idx_mom(xc, yc, M_UY_INDEX, NX_COARSE)];

            const dfloat mxx = moments_coarse[idx_mom(xc, yc, M_MXX_INDEX, NX_COARSE)];
            const dfloat mxy = moments_coarse[idx_mom(xc, yc, M_MXY_INDEX, NX_COARSE)];
            const dfloat myy = moments_coarse[idx_mom(xc, yc, M_MYY_INDEX, NX_COARSE)];

            dfloat fi[9], fi_eq[9], fi_neq[9];

            regularization(fi, rho, ux, uy, mxx, mxy, myy);
            eval_pop_eq(fi_eq, rho, ux, uy);

            for (int i = 0; i < 9; ++i)
            {
                fi_neq[i] = fi[i] - fi_eq[i];
            }

            for (int i = 0; i < 9; ++i)
            {
                pop[i] = fi_eq[i] + INV_ALPHA * fi_neq[i];
            }
        }
        else
        {
            if (y == 1)
            {
            }
            else if (y == NY_COARSE - 2)
            {
            }
            else
            {
                const int yc_mh = static_cast<int>(y / 2);
                const int yc_ph = yc_mh + 1;

                const int yc_m3h = yc_mh - 1;
                const int yc_p3h = yc_mh + 2;

                const dfloat rho_ymh = moments_coarse[idx_mom(xc, yc_mh, M_RHO_INDEX, NX_COARSE)];
                const dfloat rho_yph = moments_coarse[idx_mom(xc, yc_ph, M_RHO_INDEX, NX_COARSE)];
                const dfloat rho_ym3h = moments_coarse[idx_mom(xc, yc_m3h, M_RHO_INDEX, NX_COARSE)];
                const dfloat rho_yp3h = moments_coarse[idx_mom(xc, yc_p3h, M_RHO_INDEX, NX_COARSE)];

                const dfloat ux_ymh = moments_coarse[idx_mom(xc, yc_mh, M_UX_INDEX, NX_COARSE)];
                const dfloat ux_yph = moments_coarse[idx_mom(xc, yc_ph, M_UX_INDEX, NX_COARSE)];
                const dfloat ux_ym3h = moments_coarse[idx_mom(xc, yc_m3h, M_UX_INDEX, NX_COARSE)];
                const dfloat ux_yp3h = moments_coarse[idx_mom(xc, yc_p3h, M_UX_INDEX, NX_COARSE)];

                const dfloat uy_ymh = moments_coarse[idx_mom(xc, yc_mh, M_UX_INDEX, NX_COARSE)];
                const dfloat uy_yph = moments_coarse[idx_mom(xc, yc_ph, M_UX_INDEX, NX_COARSE)];
                const dfloat uy_ym3h = moments_coarse[idx_mom(xc, yc_m3h, M_UX_INDEX, NX_COARSE)];
                const dfloat uy_yp3h = moments_coarse[idx_mom(xc, yc_p3h, M_UX_INDEX, NX_COARSE)];

                const dfloat mxx_ymh = moments_coarse[idx_mom(xc, yc_mh, M_MXX_INDEX, NX_COARSE)];
                const dfloat mxx_yph = moments_coarse[idx_mom(xc, yc_ph, M_MXX_INDEX, NX_COARSE)];
                const dfloat mxx_ym3h = moments_coarse[idx_mom(xc, yc_m3h, M_MXX_INDEX, NX_COARSE)];
                const dfloat mxx_yp3h = moments_coarse[idx_mom(xc, yc_p3h, M_MXX_INDEX, NX_COARSE)];

                const dfloat mxy_ymh = moments_coarse[idx_mom(xc, yc_mh, M_MXY_INDEX, NX_COARSE)];
                const dfloat mxy_yph = moments_coarse[idx_mom(xc, yc_ph, M_MXY_INDEX, NX_COARSE)];
                const dfloat mxy_ym3h = moments_coarse[idx_mom(xc, yc_m3h, M_MXY_INDEX, NX_COARSE)];
                const dfloat mxy_yp3h = moments_coarse[idx_mom(xc, yc_p3h, M_MXY_INDEX, NX_COARSE)];

                const dfloat myy_ymh = moments_coarse[idx_mom(xc, yc_mh, M_MYY_INDEX, NX_COARSE)];
                const dfloat myy_yph = moments_coarse[idx_mom(xc, yc_ph, M_MYY_INDEX, NX_COARSE)];
                const dfloat myy_ym3h = moments_coarse[idx_mom(xc, yc_m3h, M_MYY_INDEX, NX_COARSE)];
                const dfloat myy_yp3h = moments_coarse[idx_mom(xc, yc_p3h, M_MYY_INDEX, NX_COARSE)];
            }
        }

        const dfloat rho = pop[0] + pop[1] + pop[2] + pop[3] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8];
        const dfloat invRho = static_cast<dfloat>(1) / rho;

        const dfloat ux = ((pop[1] + pop[5] + pop[8]) - (pop[3] + pop[6] + pop[7])) * invRho;
        const dfloat uy = ((pop[2] + pop[5] + pop[6]) - (pop[4] + pop[7] + pop[8])) * invRho;

        const dfloat mxx = (pop[1] + pop[3] + pop[5] + pop[6] + pop[7] + pop[8]) * invRho - cs2;
        const dfloat mxy = ((pop[5] + pop[7]) - (pop[6] + pop[8])) * invRho;
        const dfloat myy = (pop[2] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8]) * invRho - cs2;

        moments_fine[idx_mom(xf, y, M_RHO_INDEX, NX_COARSE)] = rho - RHO_0;

        moments_fine[idx_mom(xf, y, M_UX_INDEX, NX_COARSE)] = ux;
        moments_fine[idx_mom(xf, y, M_UY_INDEX, NX_COARSE)] = uy;

        moments_fine[idx_mom(xf, y, M_MXX_INDEX, NX_COARSE)] = mxx;
        moments_fine[idx_mom(xf, y, M_MXY_INDEX, NX_COARSE)] = mxy;
        moments_fine[idx_mom(xf, y, M_MYY_INDEX, NX_COARSE)] = myy;
    }
}

__host__ inline void fine_to_coarse(dfloat *moments_fine, dfloat *moments_coarse, unsigned int *node_type_fine, unsigned int *node_type_coarse)
{
    for (int y = 0; y < NY_COARSE; ++y)
    {
        const int xc = NX_COARSE - 1;
        const int xf = 2;

        const int yf = y * 2;

        const dfloat rho_fine = moments_fine[idx_mom(xf, yf, M_RHO_INDEX, NX_FINE)] + RHO_0;

        const dfloat ux_fine = moments_fine[idx_mom(xf, yf, M_UX_INDEX, NX_FINE)];
        const dfloat uy_fine = moments_fine[idx_mom(xf, yf, M_UY_INDEX, NX_FINE)];

        dfloat pop[9];
        dfloat pop_eq[9];
        dfloat pop_neq_filtered[9];

        eval_pop_eq(pop_eq, rho_fine, ux_fine, uy_fine);

        int count = 0;

        for (int i = 0; i < 9; ++i)
        {
            int xn = xf + cx[i];
            int yn = yf + cy[i];

            if (yn >= 0 && yn < NY_FINE)
            {
                dfloat fi[9];
                dfloat fi_eq[9];

                const dfloat rho = moments_fine[idx_mom(xn, yn, M_RHO_INDEX, NX_FINE)] + RHO_0;

                const dfloat ux = moments_fine[idx_mom(xn, yn, M_UX_INDEX, NX_FINE)];
                const dfloat uy = moments_fine[idx_mom(xn, yn, M_UY_INDEX, NX_FINE)];

                const dfloat mxx = moments_fine[idx_mom(xn, yn, M_MXX_INDEX, NX_FINE)];
                const dfloat mxy = moments_fine[idx_mom(xn, yn, M_MXY_INDEX, NX_FINE)];
                const dfloat myy = moments_fine[idx_mom(xn, yn, M_MYY_INDEX, NX_FINE)];

                regularization(fi, rho, ux, uy, mxx, mxy, myy);
                eval_pop_eq(fi_eq, rho, ux, uy);

                for (int j = 0; j < 9; ++j)
                {
                    pop_neq_filtered[j] = fi[j] - fi_eq[j];
                }

                count++;
            }
        }

        for (int i = 0; i < 9; ++i)
        {
            pop_neq_filtered[i] = pop_neq_filtered[i] / count;
        }

        for (int i = 0; i < 9; ++i)
        {
            pop[i] = pop_eq[i] + ALPHA * pop_neq_filtered[i];
        }

        const dfloat rho = pop[0] + pop[1] + pop[2] + pop[3] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8];
        const dfloat invRho = static_cast<dfloat>(1) / rho;

        const dfloat ux = ((pop[1] + pop[5] + pop[8]) - (pop[3] + pop[6] + pop[7])) * invRho;
        const dfloat uy = ((pop[2] + pop[5] + pop[6]) - (pop[4] + pop[7] + pop[8])) * invRho;

        const dfloat mxx = (pop[1] + pop[3] + pop[5] + pop[6] + pop[7] + pop[8]) * invRho - cs2;
        const dfloat mxy = ((pop[5] + pop[7]) - (pop[6] + pop[8])) * invRho;
        const dfloat myy = (pop[2] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8]) * invRho - cs2;

        moments_coarse[idx_mom(xc, y, M_RHO_INDEX, NX_COARSE)] = rho - RHO_0;

        moments_coarse[idx_mom(xc, y, M_UX_INDEX, NX_COARSE)] = ux;
        moments_coarse[idx_mom(xc, y, M_UY_INDEX, NX_COARSE)] = uy;

        moments_coarse[idx_mom(xc, y, M_MXX_INDEX, NX_COARSE)] = mxx;
        moments_coarse[idx_mom(xc, y, M_MXY_INDEX, NX_COARSE)] = mxy;
        moments_coarse[idx_mom(xc, y, M_MYY_INDEX, NX_COARSE)] = myy;
    }
}

#endif