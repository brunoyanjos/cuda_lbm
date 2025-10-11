#ifndef EVAL_BOUNDARIES_CUH
#define EVAL_BOUNDARIES_CUH

#include "../../var.h"
#include "../../nodeTypeMap.h"
#include "../../globalFunctions.h"

__host__ inline void evaluate_boundary(unsigned int nodeType,
                                       dfloat *rhoVar,
                                       dfloat *ux, dfloat *uy,
                                       dfloat *mxx, dfloat *myy, dfloat *mxy,
                                       const dfloat *pop, dfloat *moms, int x, int y,
                                       dfloat omega, size_t nx)
{
    switch (nodeType)
    {
    case NORTH:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[3] + pop[5] + pop[6];
        const dfloat inv_rhoIn = 1.0 / rhoIn;

        const dfloat mxxIn = (pop[1] + pop[3] + pop[5] + pop[6]) * inv_rhoIn - cs2;
        const dfloat mxyIn = (pop[5] - pop[6]) * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[5] + pop[6]) * inv_rhoIn - cs2;

        *ux = 0.0;
        *uy = 0.0;

        *rhoVar = 6.0 * rhoIn / 5.0;

        *mxx = 0.0;
        *mxy = 5.0 * mxyIn / 3.0;
        *myy = 0.0;

        break;
    }
    case SOUTH:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[3] + pop[4] + pop[7] + pop[8];
        const dfloat inv_rhoIn = 1.0 / rhoIn;

        const dfloat mxxIn = (pop[1] + pop[3] + pop[7] + pop[8]) * inv_rhoIn - cs2;
        const dfloat mxyIn = (pop[7] - pop[8]) * inv_rhoIn;
        const dfloat myyIn = (pop[3] + pop[7] + pop[8]) * inv_rhoIn - cs2;

        *ux = 0.0;
        *uy = 0.0;

        *rhoVar = 6.0 * rhoIn / 5.0;

        *mxx = 0.0;
        *mxy = 5.0 * mxyIn / 3.0;
        *myy = 0.0;

        break;
    }
    case WEST:
    {
        const dfloat rhoIn = pop[0] + pop[2] + pop[3] + pop[4] + pop[6] + pop[7];
        const dfloat inv_rhoIn = 1.0 / rhoIn;

        const dfloat mxxIn = (pop[3] + pop[6] + pop[7]) * inv_rhoIn - cs2;
        const dfloat mxyIn = (pop[7] - pop[6]) * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[4] + pop[6] + pop[7]) * inv_rhoIn - cs2;

        // const dfloat H = NY - 1;
        *uy = 0.0;
        // *ux = 4.0 * U_MAX * ((y / H) - (y / H) * (y / H));
        *ux = U_MAX;

        const dfloat rho = (4.0 * rhoIn + 3.0 * rhoIn * mxxIn) / (3.0 - 3.0 * (*ux));

        *mxx = (rho + 9.0 * rhoIn * mxxIn + 3.0 * rho * (*ux)) / (6.0 * rho);
        *mxy = 2.0 * rhoIn * mxyIn / rho;
        *myy = 6.0 * rhoIn * myyIn / (5.0 * rho);

        *rhoVar = rho;

        break;
    }
    case EAST:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[4] + pop[5] + pop[8];
        const dfloat inv_rhoIn = 1.0 / rhoIn;

        const dfloat mxxIn = (pop[1] + pop[5] + pop[8]) * inv_rhoIn - cs2;
        const dfloat mxyIn = (pop[5] - pop[8]) * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[4] + pop[5] + pop[8]) * inv_rhoIn - cs2;

        const dfloat rho = RHO_0 + moms[idx_mom(x - 1, y, M_RHO_INDEX, nx)];
        *ux = moms[idx_mom(x - 1, y, M_UX_INDEX, nx)] / F_M_I_SCALE;
        *uy = moms[idx_mom(x - 1, y, M_UY_INDEX, nx)] / F_M_I_SCALE;

        *rhoVar = RHO_0;

        *mxx = (rho + 9.0 * rhoIn * mxxIn - 3.0 * rho * *ux) / (6.0 * rho);
        *mxy = (6.0 * rhoIn * mxyIn - rho * *uy) / (3.0 * rho);
        *myy = 6.0 * rhoIn * myyIn / (5.0 * rho);

        break;
    }

    case SOUTH_WEST:
    {
        const dfloat rhoIn = pop[0] + pop[3] + pop[4] + pop[7];
        const dfloat inv_rhoIn = 1.0 / rhoIn;

        const dfloat mxxIn = (pop[3] + pop[7]) * inv_rhoIn - cs2;
        const dfloat mxyIn = pop[7] * inv_rhoIn;
        const dfloat myyIn = (pop[4] + pop[7]) * inv_rhoIn - cs2;

        *ux = 0.0;
        *uy = 0.0;

        *rhoVar = 36.0 * (rhoIn - mxyIn * rhoIn + mxyIn * omega * rhoIn) / (24.0 + omega);

        *mxx = 0.0;
        *mxy = (36.0 * mxyIn * rhoIn - (*rhoVar)) / (9.0 * (*rhoVar));
        *myy = 0.0;

        break;
    }
    case SOUTH_EAST:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[4] + pop[8];
        const dfloat inv_rhoIn = 1.0 / rhoIn;

        const dfloat mxxIn = (pop[1] + pop[8]) * inv_rhoIn - cs2;
        const dfloat mxyIn = -pop[8] * inv_rhoIn;
        const dfloat myyIn = (pop[4] + pop[8]) * inv_rhoIn - cs2;

        *ux = 0.0;
        *uy = 0.0;

        *rhoVar = -36.0 * (mxyIn * omega * rhoIn - rhoIn - mxyIn * rhoIn) / (24 + omega);

        *mxx = 0.0;
        *mxy = (36.0 * mxyIn * rhoIn + (*rhoVar)) / (9.0 * (*rhoVar));
        *myy = 0.0;

        break;
    }
    case NORTH_WEST:
    {
        const dfloat rhoIn = pop[0] + pop[2] + pop[3] + pop[6];
        const dfloat inv_rhoIn = 1.0 / rhoIn;

        const dfloat mxxIn = (pop[3] + pop[6]) * inv_rhoIn - cs2;
        const dfloat mxyIn = -pop[6] * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[6]) * inv_rhoIn - cs2;

        *ux = 0.0;
        *uy = 0.0;

        *rhoVar = -36.0 * (mxyIn * omega * rhoIn - rhoIn - mxyIn * rhoIn) / (24 + omega);

        *mxx = 0.0;
        *mxy = (36.0 * mxyIn * rhoIn + (*rhoVar)) / (9.0 * (*rhoVar));
        *myy = 0.0;

        break;
    }
    case NORTH_EAST:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[5];
        const dfloat inv_rhoIn = 1.0 / rhoIn;

        const dfloat mxxIn = (pop[1] + pop[5]) * inv_rhoIn - cs2;
        const dfloat mxyIn = pop[5] * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[5]) * inv_rhoIn - cs2;

        *ux = 0.0;
        *uy = 0.0;

        *rhoVar = 36.0 * (rhoIn - mxyIn * rhoIn + mxyIn * omega * rhoIn) / (24.0 + omega);

        *mxx = 0.0;
        *mxy = (36.0 * mxyIn * rhoIn - (*rhoVar)) / (9.0 * (*rhoVar));
        *myy = 0.0;

        break;
    }
    default:
        break;
    }
}

#endif