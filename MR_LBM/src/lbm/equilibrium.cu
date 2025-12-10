#include "equilibrium.cuh"

__device__ void equilibrium(dfloat *pop, const dfloat &rho, const dfloat &ux, const dfloat &uy)
{
    const dfloat u2 = ux * ux + uy * uy;
    const dfloat one_minus_u2 = dfloat(1) - cs2 * u2 / 2;

    dfloat multiplyTerm = W0 * rho;
    pop[0] = multiplyTerm * one_minus_u2;

    multiplyTerm = W1 * rho;
    pop[1] = multiplyTerm * (one_minus_u2 + ux + ux * ux / 2);
    pop[2] = multiplyTerm * (one_minus_u2 + uy + uy * uy / 2);
    pop[3] = multiplyTerm * (one_minus_u2 - ux + ux * ux / 2);
    pop[4] = multiplyTerm * (one_minus_u2 - uy + uy * uy / 2);

    multiplyTerm = W2 * rho;
    pop[5] = multiplyTerm * (one_minus_u2 + ux + uy + u2 / 2 + ux * uy);
    pop[6] = multiplyTerm * (one_minus_u2 - ux + uy + u2 / 2 - ux * uy);
    pop[7] = multiplyTerm * (one_minus_u2 - ux - uy + u2 / 2 + ux * uy);
    pop[8] = multiplyTerm * (one_minus_u2 + ux - uy + u2 / 2 - ux * uy);
}
