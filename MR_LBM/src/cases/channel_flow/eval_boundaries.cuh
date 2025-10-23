#ifndef EVAL_BOUNDARIES_CUH
#define EVAL_BOUNDARIES_CUH

#include "../../var.h"
#include "../../nodeTypeMap.h"
#include "../../globalFunctions.h"

__host__ inline void evaluate_boundary(unsigned int nodeType,
                                       dfloat *rhoVar,
                                       dfloat *ux, dfloat *uy,
                                       dfloat *mxx, dfloat *mxy, dfloat *myy,
                                       const dfloat *pop, dfloat *moms, int x, int y,
                                       dfloat omega, size_t nx)
{
    switch (nodeType)
    {
    case NORTH:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[3] + pop[5] + pop[6];
        const dfloat inv_rhoIn = static_cast<dfloat>(1) / rhoIn;

        const dfloat mxxIn = (pop[1] + pop[3] + pop[5] + pop[6]) * inv_rhoIn - cs2;
        const dfloat mxyIn = (pop[5] - pop[6]) * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[5] + pop[6]) * inv_rhoIn - cs2;

        *ux = static_cast<dfloat>(0);
        *uy = static_cast<dfloat>(0);

        *rhoVar = static_cast<dfloat>(6) * rhoIn / static_cast<dfloat>(5);

        *mxx = static_cast<dfloat>(0);
        *mxy = static_cast<dfloat>(5) * mxyIn / static_cast<dfloat>(3);
        *myy = static_cast<dfloat>(0);

        break;
    }
    case SOUTH:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[3] + pop[4] + pop[7] + pop[8];
        const dfloat inv_rhoIn = static_cast<dfloat>(1) / rhoIn;

        const dfloat mxxIn = (pop[1] + pop[3] + pop[7] + pop[8]) * inv_rhoIn - cs2;
        const dfloat mxyIn = (pop[7] - pop[8]) * inv_rhoIn;
        const dfloat myyIn = (pop[3] + pop[7] + pop[8]) * inv_rhoIn - cs2;

        *ux = static_cast<dfloat>(0);
        *uy = static_cast<dfloat>(0);

        *rhoVar = static_cast<dfloat>(6) * rhoIn / static_cast<dfloat>(5);

        *mxx = static_cast<dfloat>(0);
        *mxy = static_cast<dfloat>(5) * mxyIn / static_cast<dfloat>(3);
        *myy = static_cast<dfloat>(0);

        break;
    }
    case WEST:
    {
        const dfloat rhoIn = pop[0] + pop[2] + pop[3] + pop[4] + pop[6] + pop[7];
        const dfloat inv_rhoIn = static_cast<dfloat>(1) / rhoIn;

        const dfloat mxxIn = (pop[3] + pop[6] + pop[7]) * inv_rhoIn - cs2;
        const dfloat mxyIn = (pop[7] - pop[6]) * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[4] + pop[6] + pop[7]) * inv_rhoIn - cs2;

        *ux = U_MAX;
        *uy = static_cast<dfloat>(0);

        const dfloat rho = (static_cast<dfloat>(4) * rhoIn + static_cast<dfloat>(3) * rhoIn * mxxIn) /
                           (static_cast<dfloat>(3) - static_cast<dfloat>(3) * (*ux));

        *mxx = (rho + static_cast<dfloat>(9) * rhoIn * mxxIn + static_cast<dfloat>(3) * rho * (*ux)) /
               (static_cast<dfloat>(6) * rho);
        *mxy = static_cast<dfloat>(2) * rhoIn * mxyIn / rho;
        *myy = static_cast<dfloat>(6) * rhoIn * myyIn / (static_cast<dfloat>(5) * rho);

        *rhoVar = rho;

        break;
    }
    case EAST:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[4] + pop[5] + pop[8];
        const dfloat inv_rhoIn = static_cast<dfloat>(1) / rhoIn;

        const dfloat mxxIn = (pop[1] + pop[5] + pop[8]) * inv_rhoIn - cs2;
        const dfloat mxyIn = (pop[5] - pop[8]) * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[4] + pop[5] + pop[8]) * inv_rhoIn - cs2;

        const dfloat rho = RHO_0 + moms[idx_mom(x - 1, y, M_RHO_INDEX, nx)];

        *ux = moms[idx_mom(x - 1, y, M_UX_INDEX, nx)] / F_M_I_SCALE;
        *uy = moms[idx_mom(x - 1, y, M_UY_INDEX, nx)] / F_M_I_SCALE;

        *mxx = (rho + static_cast<dfloat>(9) * rhoIn * mxxIn - static_cast<dfloat>(3) * rho * *ux) /
               (static_cast<dfloat>(6) * rho);

        *mxy = (static_cast<dfloat>(6) * rhoIn * mxyIn - rho * *uy) /
               (static_cast<dfloat>(3) * rho);

        *myy = static_cast<dfloat>(6) * rhoIn * myyIn /
               (static_cast<dfloat>(5) * rho);

        *rhoVar = rho;

        break;
    }
    case SOUTH_WEST:
    {
        const dfloat rhoIn = pop[0] + pop[3] + pop[4] + pop[7];
        const dfloat inv_rhoIn = static_cast<dfloat>(1) / rhoIn;

        const dfloat mxxIn = (pop[3] + pop[7]) * inv_rhoIn - cs2;
        const dfloat mxyIn = pop[7] * inv_rhoIn;
        const dfloat myyIn = (pop[4] + pop[7]) * inv_rhoIn - cs2;

        *ux = static_cast<dfloat>(0);
        *uy = static_cast<dfloat>(0);

        *rhoVar = static_cast<dfloat>(36) * (rhoIn - mxyIn * rhoIn + mxyIn * omega * rhoIn) /
                  (static_cast<dfloat>(24) + omega);

        *mxx = static_cast<dfloat>(0);
        *mxy = (static_cast<dfloat>(36) * mxyIn * rhoIn - (*rhoVar)) /
               (static_cast<dfloat>(9) * (*rhoVar));
        *myy = static_cast<dfloat>(0);

        break;
    }
    case SOUTH_EAST:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[4] + pop[8];
        const dfloat inv_rhoIn = static_cast<dfloat>(1) / rhoIn;

        const dfloat mxxIn = (pop[1] + pop[8]) * inv_rhoIn - cs2;
        const dfloat mxyIn = -pop[8] * inv_rhoIn;
        const dfloat myyIn = (pop[4] + pop[8]) * inv_rhoIn - cs2;

        *ux = static_cast<dfloat>(0);
        *uy = static_cast<dfloat>(0);

        *rhoVar = -static_cast<dfloat>(36) * (mxyIn * omega * rhoIn - rhoIn - mxyIn * rhoIn) / (static_cast<dfloat>(24) + omega);

        *mxx = static_cast<dfloat>(0);
        *mxy = (static_cast<dfloat>(36) * mxyIn * rhoIn + (*rhoVar)) / (static_cast<dfloat>(9) * (*rhoVar));
        *myy = static_cast<dfloat>(0);

        break;
    }
    case NORTH_WEST:
    {
        const dfloat rhoIn = pop[0] + pop[2] + pop[3] + pop[6];
        const dfloat inv_rhoIn = static_cast<dfloat>(1) / rhoIn;

        const dfloat mxxIn = (pop[3] + pop[6]) * inv_rhoIn - cs2;
        const dfloat mxyIn = -pop[6] * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[6]) * inv_rhoIn - cs2;

        *ux = static_cast<dfloat>(0);
        *uy = static_cast<dfloat>(0);

        *rhoVar = -static_cast<dfloat>(36) * (mxyIn * omega * rhoIn - rhoIn - mxyIn * rhoIn) / (static_cast<dfloat>(24) + omega);

        *mxx = static_cast<dfloat>(0);
        *mxy = (static_cast<dfloat>(36) * mxyIn * rhoIn + (*rhoVar)) / (static_cast<dfloat>(9) * (*rhoVar));
        *myy = static_cast<dfloat>(0);

        break;
    }
    case NORTH_EAST:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[5];
        const dfloat inv_rhoIn = static_cast<dfloat>(1) / rhoIn;

        const dfloat mxxIn = (pop[1] + pop[5]) * inv_rhoIn - cs2;
        const dfloat mxyIn = pop[5] * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[5]) * inv_rhoIn - cs2;

        *ux = static_cast<dfloat>(0);
        *uy = static_cast<dfloat>(0);

        *rhoVar = static_cast<dfloat>(36) * (rhoIn - mxyIn * rhoIn + mxyIn * omega * rhoIn) /
                  (static_cast<dfloat>(24) + omega);

        *mxx = static_cast<dfloat>(0);
        *mxy = (static_cast<dfloat>(36) * mxyIn * rhoIn - (*rhoVar)) /
               (static_cast<dfloat>(24) * (*rhoVar));
        *myy = static_cast<dfloat>(0);

        break;
    }
    default:
        break;
    }
}

#endif