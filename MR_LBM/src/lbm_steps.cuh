#ifndef LBM_STEPS_CUH
#define LBM_STEPS_CUH

#include "var.h"

__host__ inline void init_pop_in(latticeNode *node)
{
    dfloat rho = (*node).rho;
    dfloat ux = (*node).ux;
    dfloat uy = (*node).uy;
    dfloat mxx = (*node).mxx;
    dfloat mxy = (*node).mxy;
    dfloat myy = (*node).myy;

    dfloat pics2 = 1 - cs2 * (mxx + myy);

    dfloat multiplyTerm = W0 * rho;
    (*node).pop_in[0] = multiplyTerm * (pics2);

    multiplyTerm = W1 * rho;
    (*node).pop_in[1] = multiplyTerm * (pics2 + ux + mxx);
    (*node).pop_in[2] = multiplyTerm * (pics2 + uy + myy);
    (*node).pop_in[3] = multiplyTerm * (pics2 - ux + mxx);
    (*node).pop_in[4] = multiplyTerm * (pics2 - uy + myy);

    multiplyTerm = W2 * rho;
    (*node).pop_in[5] = multiplyTerm * (pics2 + ux + uy + mxx + myy + mxy);
    (*node).pop_in[6] = multiplyTerm * (pics2 - ux + uy + mxx + myy - mxy);
    (*node).pop_in[7] = multiplyTerm * (pics2 - ux - uy + mxx + myy + mxy);
    (*node).pop_in[8] = multiplyTerm * (pics2 + ux - uy + mxx + myy - mxy);
}

__host__ inline void regularization(latticeNode *node)
{
    dfloat rho = (*node).rho;
    dfloat ux = (*node).ux;
    dfloat uy = (*node).uy;
    dfloat mxx = (*node).mxx;
    dfloat mxy = (*node).mxy;
    dfloat myy = (*node).myy;

    dfloat pics2 = 1 - cs2 * (mxx + myy);

    dfloat multiplyTerm = W0 * rho;
    (*node).pop_out[0] = multiplyTerm * (pics2);

    multiplyTerm = W1 * rho;
    (*node).pop_out[1] = multiplyTerm * (pics2 + ux + mxx);
    (*node).pop_out[2] = multiplyTerm * (pics2 + uy + myy);
    (*node).pop_out[3] = multiplyTerm * (pics2 - ux + mxx);
    (*node).pop_out[4] = multiplyTerm * (pics2 - uy + myy);

    multiplyTerm = W2 * rho;
    (*node).pop_out[5] = multiplyTerm * (pics2 + ux + uy + mxx + myy + mxy);
    (*node).pop_out[6] = multiplyTerm * (pics2 - ux + uy + mxx + myy - mxy);
    (*node).pop_out[7] = multiplyTerm * (pics2 - ux - uy + mxx + myy + mxy);
    (*node).pop_out[8] = multiplyTerm * (pics2 + ux - uy + mxx + myy - mxy);
}

__host__ inline void boundary_condition(latticeNode *node)
{
    unsigned int nodeType = (*node).node_type;
    const dfloat *pop = (*node).pop_in;

    switch (nodeType)
    {
    case NORTH:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[3] + pop[5] + pop[6];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat mxxIn = (pop[1] + pop[3] + pop[5] + pop[6]) * inv_rhoIn - cs2;
        const dfloat mxyIn = (pop[5] - pop[6]) * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[5] + pop[6]) * inv_rhoIn - cs2;

        (*node).ux = U_MAX;
        (*node).uy = 0.0f;

        const dfloat rhoVar = 6.0f * rhoIn / 5.0f;

        (*node).mxx = U_MAX * U_MAX;
        (*node).mxy = 5.0f * mxyIn / 3.0f - U_MAX / 3.0f;
        (*node).myy = 0.0f;

        (*node).rho = rhoVar;

        break;
    }
    case SOUTH:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[3] + pop[4] + pop[7] + pop[8];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat mxxIn = (pop[1] + pop[3] + pop[7] + pop[8]) * inv_rhoIn - cs2;
        const dfloat mxyIn = (pop[7] - pop[8]) * inv_rhoIn;
        const dfloat myyIn = (pop[3] + pop[7] + pop[8]) * inv_rhoIn - cs2;

        (*node).ux = 0.0f;
        (*node).uy = 0.0f;

        const dfloat rhoVar = 6.0f * rhoIn / 5.0f;

        (*node).mxx = 0.0f;
        (*node).mxy = 5.0f * mxyIn / 3.0f;
        (*node).myy = 0.0f;

        (*node).rho = rhoVar;

        break;
    }
    case WEST:
    {
        const dfloat rhoIn = pop[0] + pop[2] + pop[3] + pop[4] + pop[6] + pop[7];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat mxxIn = (pop[3] + pop[6] + pop[7]) * inv_rhoIn - cs2;
        const dfloat mxyIn = (pop[7] - pop[6]) * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[4] + pop[6] + pop[7]) * inv_rhoIn - cs2;

        (*node).ux = 0.0f;
        (*node).uy = 0.0f;

        const dfloat rhoVar = 6.0f * rhoIn / 5.0f;

        (*node).mxx = 0.0f;
        (*node).mxy = 5.0f * mxyIn / 3.0f;
        (*node).myy = 0.0f;

        (*node).rho = rhoVar;

        break;
    }
    case EAST:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[4] + pop[5] + pop[8];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat mxxIn = (pop[1] + pop[5] + pop[8]) * inv_rhoIn - cs2;
        const dfloat mxyIn = (pop[5] - pop[8]) * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[4] + pop[5] + pop[8]) * inv_rhoIn - cs2;

        (*node).ux = 0.0f;
        (*node).uy = 0.0f;

        const dfloat rhoVar = 6.0f * rhoIn / 5.0f;

        (*node).mxx = 0.0f;
        (*node).mxy = 5.0f * mxyIn / 3.0f;
        (*node).myy = 0.0f;

        (*node).rho = rhoVar;

        break;
    }
    case SOUTH_WEST:
    {
        const dfloat rhoIn = pop[0] + pop[3] + pop[4] + pop[7];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat mxxIn = (pop[3] + pop[7]) * inv_rhoIn - cs2;
        const dfloat mxyIn = pop[7] * inv_rhoIn;
        const dfloat myyIn = (pop[4] + pop[7]) * inv_rhoIn - cs2;

        (*node).ux = 0.0f;

        const dfloat rhoVar = 36.0f * (rhoIn - mxyIn * rhoIn + mxyIn * OMEGA * rhoIn) / (24.0f + OMEGA);

        (*node).mxx = 0.0f;
        (*node).mxy = (36.0f * mxyIn * rhoIn - (rhoVar)) / (9.0f * (rhoVar));
        (*node).myy = 0.0f;

        (*node).rho = rhoVar;

        break;
    }
    case SOUTH_EAST:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[4] + pop[8];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat mxxIn = (pop[1] + pop[8]) * inv_rhoIn - cs2;
        const dfloat mxyIn = -pop[8] * inv_rhoIn;
        const dfloat myyIn = (pop[4] + pop[8]) * inv_rhoIn - cs2;

        (*node).uy = 0.0f;

        const dfloat rhoVar = -36.0f * (mxyIn * OMEGA * rhoIn - rhoIn - mxyIn * rhoIn) / (24 + OMEGA);

        (*node).mxx = 0.0f;
        (*node).mxy = (36.0f * mxyIn * rhoIn + (rhoVar)) / (9.0f * (rhoVar));
        (*node).myy = 0.0f;

        (*node).rho = rhoVar;

        break;
    }
    case NORTH_WEST:
    {
        const dfloat rhoIn = pop[0] + pop[2] + pop[3] + pop[6];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat mxxIn = (pop[3] + pop[6]) * inv_rhoIn - cs2;
        const dfloat mxyIn = -pop[6] * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[6]) * inv_rhoIn - cs2;

        (*node).ux = U_MAX;
        (*node).uy = 0.0f;

        const dfloat rhoVar = -36.0f * (mxyIn * OMEGA * rhoIn - rhoIn - mxyIn * rhoIn) / (24.0f + OMEGA + 18.0f * U_MAX - 3.0f * OMEGA * U_MAX - 18.0f * U_MAX * U_MAX + 3.0f * OMEGA * U_MAX * U_MAX);

        (*node).mxx = U_MAX * U_MAX;
        (*node).mxy = (36.0f * mxyIn * rhoIn + (rhoVar)-3.0f * U_MAX * (rhoVar) + 3.0f * U_MAX * U_MAX * (rhoVar)) / (9.0f * (rhoVar));
        (*node).myy = 0.0f;

        (*node).rho = rhoVar;

        break;
    }
    case NORTH_EAST:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[5];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat mxxIn = (pop[1] + pop[5]) * inv_rhoIn - cs2;
        const dfloat mxyIn = pop[5] * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[5]) * inv_rhoIn - cs2;

        (*node).ux = U_MAX;
        (*node).uy = 0.0f;

        const dfloat rhoVar = 36.0f * (mxyIn * OMEGA * rhoIn + rhoIn - mxyIn * rhoIn) / (24.0f + OMEGA - 18.0f * U_MAX + 3.0f * OMEGA * U_MAX - 18.0f * U_MAX * U_MAX + 3.0f * OMEGA * U_MAX * U_MAX);

        (*node).mxx = U_MAX * U_MAX;
        (*node).mxy = (36.0f * mxyIn * rhoIn - (rhoVar)-3.0f * U_MAX * (rhoVar)-3.0f * U_MAX * U_MAX * (rhoVar)) / (9.0f * (rhoVar));
        (*node).myy = 0.0f;

        (*node).rho = rhoVar;

        break;
    }
    case BB_EAST:
    {
        const dfloat rhoIn = pop[0] + pop[2] + pop[3] + pop[4] + pop[6] + pop[7];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat uxIn = -(pop[3] + pop[6] + pop[7]) * inv_rhoIn;
        const dfloat uyIn = ((pop[2] + pop[6]) - (pop[4] + pop[7])) * inv_rhoIn;

        const dfloat mxxIn = (pop[3] + pop[6] + pop[7]) * inv_rhoIn - cs2;
        const dfloat mxyIn = (pop[7] - pop[6]) * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[4] + pop[6] + pop[7]) * inv_rhoIn - cs2;

        // ---------------------------------------------------------
        // Etapa 1: Calcular rho
        // ---------------------------------------------------------
        const dfloat A = -36.0f * rhoIn + 72.0f * mxxIn * rhoIn - 72.0f * myyIn * rhoIn + 252.0f * rhoIn * uxIn;

        const dfloat B =
            60.0f * rhoIn * rhoIn +
            135.0f * mxyIn * mxyIn * rhoIn * rhoIn -
            120.0f * rhoIn * rhoIn * uxIn +
            60.0f * rhoIn * rhoIn * uxIn * uxIn +
            270.0f * mxyIn * rhoIn * rhoIn * uyIn +
            135.0f * rhoIn * rhoIn * uyIn * uyIn;

        const dfloat rhoVar = (1.0f / 40.0f) * (36.0f * rhoIn - 72.0f * mxxIn * rhoIn + 72.0f * myyIn * rhoIn - 252.0f * rhoIn * uxIn +
                                                sqrtf(A * A - 80.0f * B));

        // ---------------------------------------------------------
        // Etapa 2: Calcular ux
        // ---------------------------------------------------------
        const dfloat C =
            135.0f * mxyIn * mxyIn * rhoIn * rhoIn +
            270.0f * mxyIn * rhoIn * uyIn +
            135.0f * rhoIn * rhoIn * uyIn * uyIn +
            72.0f * mxxIn * rhoIn * rhoVar -
            72.0f * myyIn * rhoIn * rhoVar +
            216.0f * rhoIn * uxIn * rhoVar +
            44.0f * rhoVar * rhoVar;

        const dfloat uxVar = (84.0f * rhoVar * rhoVar + sqrtf(7056.0f * pow(rhoVar, 4) - 240.0f * rhoVar * rhoVar * C)) / (120.0f * rhoVar * rhoVar);

        // ---------------------------------------------------------
        // Etapa 3: Calcular uy
        // ---------------------------------------------------------
        const dfloat uyVar = (3.0f * (mxyIn * rhoIn + rhoIn * uyIn)) / (2.0f * rhoVar);

        // ---------------------------------------------------------
        // Etapa 4: Calcular mxx
        // ---------------------------------------------------------
        const dfloat D =
            675.0f * mxyIn * mxyIn * rhoIn * rhoIn +
            1350.0f * mxyIn * rhoIn * uyIn +
            675.0f * rhoIn * rhoIn * uyIn * uyIn +
            360.0f * mxxIn * rhoIn * rhoVar -
            360.0f * myyIn * rhoIn * rhoVar +
            1080.0f * rhoIn * uxIn * rhoVar +
            73.0f * rhoVar * rhoVar;

        const dfloat mxxVar = (-60.0f * rhoIn * uxIn + 11.0f * rhoVar + sqrtf(3.0f) * sqrtf(-rhoVar * rhoVar) * D) / (30.0f * rhoVar);

        // ---------------------------------------------------------
        // Etapa 5: Calcular mxy
        // ---------------------------------------------------------
        const dfloat mxyVar = rhoIn * (5.0f * mxyIn + uyIn) / (2.0f * rhoVar);

        // ---------------------------------------------------------
        // Etapa 6: Calcular myy
        // ---------------------------------------------------------
        const dfloat myyVar = (-90.0f * mxxIn * rhoIn + 90.0f * myyIn * rhoIn - 120.0f * rhoIn * uxIn - 9.0f * rhoVar +
                               sqrtf(3.0f) * sqrtf(-rhoVar * rhoVar) * D) /
                              (75.0f * rhoVar);

        // ---------------------------------------------------------
        // Atribuir ao node
        // ---------------------------------------------------------
        (*node).rho = rhoVar;
        (*node).ux = uxVar;
        (*node).uy = uyVar;
        (*node).mxx = mxxVar;
        (*node).mxy = mxyVar;
        (*node).myy = myyVar;

        break;
    }
    case BB_WEST:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[4] + pop[5] + pop[8];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat uxIn = (pop[1] + pop[5] + pop[8]) * inv_rhoIn;
        const dfloat uyIn = ((pop[2] + pop[5]) - (pop[4] + pop[8])) * inv_rhoIn;

        const dfloat mxxIn = (pop[1] + pop[5] + pop[8]) * inv_rhoIn - cs2;
        const dfloat mxyIn = (pop[5] - pop[8]) * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[4] + pop[5] + pop[8]) * inv_rhoIn - cs2;

        // ---------------------------------------------------------
        // Etapa 1: Calcular rho
        // ---------------------------------------------------------
        const dfloat A = -36.0f * rhoIn + 72.0f * mxxIn * rhoIn - 72.0f * myyIn * rhoIn - 252.0f * rhoIn * uxIn;

        const dfloat B =
            60.0f * rhoIn * rhoIn +
            135.0f * mxyIn * mxyIn * rhoIn * rhoIn +
            120.0f * rhoIn * rhoIn * uxIn +
            60.0f * rhoIn * rhoIn * uxIn * uxIn -
            270.0f * mxyIn * rhoIn * rhoIn * uyIn +
            135.0f * rhoIn * rhoIn * uyIn * uyIn;

        const dfloat rhoVar = (1.0f / 40.0f) * (36.0f * rhoIn - 72.0f * mxxIn * rhoIn + 72.0f * myyIn * rhoIn + 252.0f * rhoIn * uxIn +
                                                sqrtf(A * A - 80.0f * B));

        // ---------------------------------------------------------
        // Etapa 2: Calcular ux
        // ---------------------------------------------------------
        const dfloat C =
            135.0f * mxyIn * mxyIn * rhoIn * rhoIn -
            270.0f * mxyIn * rhoIn * uyIn +
            135.0f * rhoIn * rhoIn * uyIn * uyIn +
            72.0f * mxxIn * rhoIn * rhoVar -
            72.0f * myyIn * rhoIn * rhoVar -
            216.0f * rhoIn * uxIn * rhoVar +
            44.0f * rhoVar * rhoVar;

        const dfloat uxVar = (-84.0f * rhoVar * rhoVar + sqrtf(7056.0f * pow(rhoVar, 4) - 240.0f * rhoVar * rhoVar * C)) / (120.0f * rhoVar * rhoVar);

        // ---------------------------------------------------------
        // Etapa 3: Calcular uy
        // ---------------------------------------------------------
        const dfloat uyVar = (3.0f * (mxyIn * rhoIn - rhoIn * uyIn)) / (2.0f * rhoVar);

        // ---------------------------------------------------------
        // Etapa 4: Calcular mxx
        // ---------------------------------------------------------
        const dfloat D =
            675.0f * mxyIn * mxyIn * rhoIn * rhoIn -
            1350.0f * mxyIn * rhoIn * uyIn +
            675.0f * rhoIn * rhoIn * uyIn * uyIn +
            360.0f * mxxIn * rhoIn * rhoVar -
            360.0f * myyIn * rhoIn * rhoVar -
            1080.0f * rhoIn * uxIn * rhoVar +
            73.0f * rhoVar * rhoVar;

        const dfloat mxxVar = (60.0f * rhoIn * uxIn + 11.0f * rhoVar - sqrtf(3.0f) * sqrtf(-rhoVar * rhoVar) * D) / (30.0f * rhoVar);

        // ---------------------------------------------------------
        // Etapa 5: Calcular mxy
        // ---------------------------------------------------------
        const dfloat mxyVar = rhoIn * (5.0f * mxyIn - uyIn) / (2.0f * rhoVar);

        // ---------------------------------------------------------
        // Etapa 6: Calcular myy
        // ---------------------------------------------------------
        const dfloat myyVar = (-90.0f * mxxIn * rhoIn + 90.0f * myyIn * rhoIn + 120.0f * rhoIn * uxIn - 9.0f * rhoVar -
                               sqrtf(3.0f) * sqrtf(-rhoVar * rhoVar) * D) /
                              (75.0f * rhoVar);

        // ---------------------------------------------------------
        // Armazenar no ponteiro node
        // ---------------------------------------------------------
        (*node).rho = rhoVar;
        (*node).ux = uxVar;
        (*node).uy = uyVar;
        (*node).mxx = mxxVar;
        (*node).mxy = mxyVar;
        (*node).myy = myyVar;

        break;
    }
    case BB_NORTH:
    {
        std::cout << "pop[0]: " << pop[0] << std::endl;
        std::cout << "pop[1]: " << pop[1] << std::endl;
        std::cout << "pop[3]: " << pop[3] << std::endl;
        std::cout << "pop[4]: " << pop[4] << std::endl;
        std::cout << "pop[7]: " << pop[7] << std::endl;
        std::cout << "pop[8]: " << pop[8] << std::endl
                  << std::endl;

        const dfloat rhoIn = pop[0] + pop[1] + pop[3] + pop[4] + pop[7] + pop[8];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat uxIn = ((pop[1] + pop[8]) - (pop[3] + pop[7])) * inv_rhoIn;
        const dfloat uyIn = -(pop[4] + pop[7] + pop[8]) * inv_rhoIn;

        const dfloat mxxIn = (pop[1] + pop[3] + pop[7] + pop[8]) * inv_rhoIn - cs2;
        const dfloat mxyIn = (pop[7] - pop[8]) * inv_rhoIn;
        const dfloat myyIn = (pop[3] + pop[7] + pop[8]) * inv_rhoIn - cs2;

        std::cout << "rhoIn: " << rhoIn << std::endl;
        std::cout << "uxIn: " << uxIn << std::endl;
        std::cout << "uyIn: " << uyIn << std::endl;
        std::cout << "mxxIn: " << mxxIn << std::endl;
        std::cout << "mxyIn: " << mxyIn << std::endl;
        std::cout << "myyIn: " << myyIn << std::endl
                  << std::endl;

        // ---------------------------------------------------------
        // Etapa 1: rho
        // ---------------------------------------------------------
        const dfloat A = -36.0f * rhoIn + 72.0f * mxxIn * rhoIn - 72.0f * myyIn * rhoIn + 252.0f * rhoIn * uyIn;

        const dfloat B =
            60.0f * rhoIn * rhoIn +
            135.0f * mxyIn * mxyIn * rhoIn * rhoIn +
            270.0f * mxyIn * rhoIn * rhoIn * uxIn +
            135.0f * rhoIn * rhoIn * uxIn * uxIn -
            120.0f * rhoIn * rhoIn * uyIn +
            60.0f * rhoIn * rhoIn * uyIn * uyIn;

        const dfloat rhoVar = (1.0f / 40.0f) * (36.0f * rhoIn + 72.0f * mxxIn * rhoIn - 72.0f * myyIn * rhoIn - 252.0f * rhoIn * uyIn +
                                                sqrtf(A * A - 80.0f * B));

        std::cout << "A: " << A << std::endl;
        std::cout << "B: " << B << std::endl;

        std::cout << "sqrt_value: " << A * A - 80.0f * B << std::endl
                  << std::endl;

        // ---------------------------------------------------------
        // Etapa 2: ux
        // ---------------------------------------------------------
        const dfloat uxVar = (3.0f * (mxyIn * rhoIn + rhoIn * uxIn)) / (2.0f * rhoVar);

        // ---------------------------------------------------------
        // Etapa 3: uy
        // ---------------------------------------------------------
        const dfloat C =
            135.0f * mxyIn * mxyIn * rhoIn * rhoIn +
            270.0f * mxyIn * rhoIn * uxIn +
            135.0f * rhoIn * rhoIn * uxIn * uxIn -
            72.0f * mxxIn * rhoIn * rhoVar +
            72.0f * myyIn * rhoIn * rhoVar +
            216.0f * rhoIn * uyIn * rhoVar +
            44.0f * rhoVar * rhoVar;

        const dfloat uyVar = (84.0f * rhoVar * rhoVar + sqrtf(7056.0f * pow(rhoVar, 4) - 240.0f * rhoVar * rhoVar * C)) / (120.0f * rhoVar * rhoVar);

        // ---------------------------------------------------------
        // Etapa 4: mxx
        // ---------------------------------------------------------
        const dfloat D =
            -675.0f * mxyIn * mxyIn * rhoIn * rhoIn -
            1350.0f * mxyIn * rhoIn * uxIn -
            675.0f * rhoIn * rhoIn * uxIn * uxIn +
            360.0f * mxxIn * rhoIn * rhoVar -
            360.0f * myyIn * rhoIn * rhoVar -
            1080.0f * rhoIn * uyIn * rhoVar -
            73.0f * rhoVar * rhoVar;

        const dfloat mxxVar = (90.0f * mxxIn * rhoIn - 90.0f * myyIn * rhoIn - 120.0f * rhoIn * uyIn - 9.0f * rhoVar +
                               sqrtf(3.0f) * sqrtf(rhoVar * rhoVar) * D) /
                              (75.0f * rhoVar);

        // ---------------------------------------------------------
        // Etapa 5: mxy
        // ---------------------------------------------------------
        const dfloat mxyVar = rhoIn * (5.0f * mxyIn + uxIn) / (2.0f * rhoVar);

        // ---------------------------------------------------------
        // Etapa 6: myy
        // ---------------------------------------------------------
        const dfloat E =
            675.0f * mxyIn * mxyIn * rhoIn * rhoIn +
            1350.0f * mxyIn * rhoIn * uxIn +
            675.0f * rhoIn * rhoIn * uxIn * uxIn +
            360.0f * mxxIn * rhoIn * rhoVar -
            360.0f * myyIn * rhoIn * rhoVar +
            1080.0f * rhoIn * uyIn * rhoVar +
            73.0f * rhoVar * rhoVar;

        const dfloat myyVar = (-60.0f * rhoIn * uyIn + 11.0f * rhoVar + sqrtf(3.0f) * sqrtf(-rhoVar * rhoVar) * E) / (30.0f * rhoVar);

        // ---------------------------------------------------------
        // Atribuir ao node
        // ---------------------------------------------------------
        (*node).rho = rhoVar;
        (*node).ux = uxVar;
        (*node).uy = uyVar;
        (*node).mxx = mxxVar;
        (*node).mxy = mxyVar;
        (*node).myy = myyVar;

        std::cout << "rho: " << rhoVar << std::endl;
        std::cout << "ux: " << uxVar << std::endl;
        std::cout << "uy: " << uyVar << std::endl;
        std::cout << "mxx: " << mxxVar << std::endl;
        std::cout << "mxy: " << mxyVar << std::endl;
        std::cout << "myy: " << myyVar << std::endl
                  << std::endl;

        break;
    }
    case BB_SOUTH:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[3] + pop[5] + pop[6];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat uxIn = ((pop[1] + pop[5]) - (pop[3] + pop[6])) * inv_rhoIn;
        const dfloat uyIn = (pop[2] + pop[5] + pop[6]) * inv_rhoIn;

        const dfloat mxxIn = (pop[1] + pop[3] + pop[5] + pop[6]) * inv_rhoIn - cs2;
        const dfloat mxyIn = (pop[5] - pop[6]) * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[5] + pop[6]) * inv_rhoIn - cs2;

        // ---------------------------------------------------------
        // Etapa 1: rho
        // ---------------------------------------------------------
        const dfloat A = -36.0f * rhoIn + 72.0f * mxxIn * rhoIn - 72.0f * myyIn * rhoIn + 252.0f * rhoIn * uyIn;

        const dfloat B =
            60.0f * rhoIn * rhoIn +
            135.0f * mxyIn * mxyIn * rhoIn * rhoIn -
            270.0f * mxyIn * rhoIn * rhoIn * uxIn +
            135.0f * rhoIn * rhoIn * uxIn * uxIn +
            120.0f * rhoIn * rhoIn * uyIn +
            60.0f * rhoIn * rhoIn * uyIn * uyIn;

        const dfloat rhoVar = (1.0f / 40.0f) * (36.0f * rhoIn + 72.0f * mxxIn * rhoIn - 72.0f * myyIn * rhoIn + 252.0f * rhoIn * uyIn +
                                                sqrtf(A * A - 80.0f * B));

        // ---------------------------------------------------------
        // Etapa 2: ux
        // ---------------------------------------------------------
        const dfloat uxVar = -(3.0f * (mxyIn * rhoIn - rhoIn * uxIn)) / (2.0f * rhoVar);

        // ---------------------------------------------------------
        // Etapa 3: uy
        // ---------------------------------------------------------
        const dfloat C =
            135.0f * mxyIn * mxyIn * rhoIn * rhoIn -
            270.0f * mxyIn * rhoIn * uxIn +
            135.0f * rhoIn * rhoIn * uxIn * uxIn -
            72.0f * mxxIn * rhoIn * rhoVar +
            72.0f * myyIn * rhoIn * rhoVar -
            216.0f * rhoIn * uyIn * rhoVar +
            44.0f * rhoVar * rhoVar;

        const dfloat uyVar = (-84.0f * rhoVar * rhoVar + sqrtf(7056.0f * pow(rhoVar, 4) - 240.0f * rhoVar * rhoVar * C)) / (120.0f * rhoVar * rhoVar);

        // ---------------------------------------------------------
        // Etapa 4: mxx
        // ---------------------------------------------------------
        const dfloat D =
            -675.0f * mxyIn * mxyIn * rhoIn * rhoIn +
            1350.0f * mxyIn * rhoIn * uxIn -
            675.0f * rhoIn * rhoIn * uxIn * uxIn +
            360.0f * mxxIn * rhoIn * rhoVar -
            360.0f * myyIn * rhoIn * rhoVar +
            1080.0f * rhoIn * uyIn * rhoVar -
            73.0f * rhoVar * rhoVar;

        const dfloat mxxVar = (90.0f * mxxIn * rhoIn - 90.0f * myyIn * rhoIn + 120.0f * rhoIn * uyIn - 9.0f * rhoVar -
                               sqrtf(3.0f) * sqrtf(rhoVar * rhoVar) * D) /
                              (75.0f * rhoVar);

        // ---------------------------------------------------------
        // Etapa 5: mxy
        // ---------------------------------------------------------
        const dfloat mxyVar = rhoIn * (5.0f * mxyIn - uxIn) / (2.0f * rhoVar);

        // ---------------------------------------------------------
        // Etapa 6: myy
        // ---------------------------------------------------------
        const dfloat E =
            675.0f * mxyIn * mxyIn * rhoIn * rhoIn -
            1350.0f * mxyIn * rhoIn * uxIn +
            675.0f * rhoIn * rhoIn * uxIn * uxIn -
            360.0f * mxxIn * rhoIn * rhoVar +
            360.0f * myyIn * rhoIn * rhoVar -
            1080.0f * rhoIn * uyIn * rhoVar +
            73.0f * rhoVar * rhoVar;

        const dfloat myyVar = (60.0f * rhoIn * uyIn + 11.0f * rhoVar - sqrtf(3.0f) * sqrtf(rhoVar * rhoVar) * E) / (30.0f * rhoVar);

        // ---------------------------------------------------------
        // Atribuir ao node
        // ---------------------------------------------------------
        (*node).rho = rhoVar;
        (*node).ux = uxVar;
        (*node).uy = uyVar;
        (*node).mxx = mxxVar;
        (*node).mxy = mxyVar;
        (*node).myy = myyVar;

        break;
    }
    case BB_NORTH_EAST:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[3] + pop[4] + pop[6] + pop[7] + pop[8];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat uxIn = ((pop[1] + pop[5] + pop[8]) - (pop[3] + pop[6] + pop[7])) * inv_rhoIn;
        const dfloat uyIn = ((pop[2] + pop[5] + pop[6]) - (pop[4] + pop[7] + pop[8])) * inv_rhoIn;

        const dfloat mxxIn = (pop[1] + pop[3] + pop[6] + pop[7] + pop[8]) * inv_rhoIn - cs2;
        const dfloat mxyIn = (pop[7] - (pop[6] + pop[8])) * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[4] + pop[6] + pop[7] + pop[8]) * inv_rhoIn - cs2;

        // ---------------------------------------------------------
        // Etapa 1: rho
        // ---------------------------------------------------------
        const dfloat A = 9.0f * rhoIn + 9.0f * mxyIn * rhoIn + 9.0f * rhoIn * uxIn + 9.0f * rhoIn * uyIn;
        const dfloat B = 6.0f * rhoIn * rhoIn - 6.0f * rhoIn * rhoIn * uxIn + 3.0f * rhoIn * rhoIn * uxIn * uxIn - 6.0f * rhoIn * rhoIn * uyIn + 3.0f * rhoIn * rhoIn * uyIn * uyIn;

        const dfloat rhoVar = (1.0f / 28.0f) * (9.0f * rhoIn + 9.0f * mxyIn * rhoIn + 9.0f * rhoIn * uxIn + 9.0f * rhoIn * uyIn +
                                                sqrtf(A * A + 56.0f * B));

        // ---------------------------------------------------------
        // Etapa 2: ux
        // ---------------------------------------------------------
        const dfloat uxNumerator =
            6.0f * rhoIn * uxIn - 6.0f * rhoIn * uyIn + 21.0f * rhoVar +
            sqrtf(3.0f) * sqrtf(rhoVar * rhoVar) * (12.0f * rhoIn * rhoIn * uxIn * uxIn - 24.0f * rhoIn * rhoIn * uxIn * uyIn + 12.0f * rhoIn * rhoIn * uyIn * uyIn + 72.0f * mxyIn * rhoIn * rhoVar + 108.0f * rhoIn * uxIn * rhoVar + 108.0f * rhoIn * uyIn * rhoVar - 139.0f * rhoVar * rhoVar);

        const dfloat uxVar = uxNumerator / (12.0f * rhoVar);

        // ---------------------------------------------------------
        // Etapa 3: uy
        // ---------------------------------------------------------
        const dfloat uyNumerator =
            -6.0f * rhoIn * uxIn + 6.0f * rhoIn * uyIn + 21.0f * rhoVar +
            sqrtf(3.0f) * sqrtf(rhoVar * rhoVar) * (12.0f * rhoIn * rhoIn * uxIn * uxIn - 24.0f * rhoIn * rhoIn * uxIn * uyIn + 12.0f * rhoIn * rhoIn * uyIn * uyIn + 72.0f * mxyIn * rhoIn * rhoVar + 108.0f * rhoIn * uxIn * rhoVar + 108.0f * rhoIn * uyIn * rhoVar - 139.0f * rhoVar * rhoVar);

        const dfloat uyVar = uyNumerator / (12.0f * rhoVar);

        // ---------------------------------------------------------
        // Etapa 4: mxx
        // ---------------------------------------------------------
        const dfloat mxxNumerator =
            12.0f * mxxIn * rhoIn - 36.0f * mxyIn * rhoIn - 12.0f * myyIn * rhoIn -
            54.0f * rhoIn * uxIn - 54.0f * rhoIn * uyIn +
            7.0f * sqrtf(3.0f) * sqrtf(rhoVar * rhoVar) * (12.0f * rhoIn * rhoIn * uxIn * uxIn - 24.0f * rhoIn * rhoIn * uxIn * uyIn + 12.0f * rhoIn * rhoIn * uyIn * uyIn + 72.0f * mxyIn * rhoIn * rhoVar + 108.0f * rhoIn * uxIn * rhoVar + 108.0f * rhoIn * uyIn * rhoVar - 139.0f * rhoVar * rhoVar);

        const dfloat mxxVar = mxxNumerator / (143.0f * rhoVar);

        // ---------------------------------------------------------
        // Etapa 5: mxy
        // ---------------------------------------------------------
        const dfloat mxyNumerator =
            12.0f * mxyIn * rhoIn - 6.0f * rhoIn * uxIn - 6.0f * rhoIn * uyIn + 21.0f * rhoVar +
            sqrtf(3.0f) * sqrtf(rhoVar * rhoVar) * (12.0f * rhoIn * rhoIn * uxIn * uxIn - 24.0f * rhoIn * rhoIn * uxIn * uyIn + 12.0f * rhoIn * rhoIn * uyIn * uyIn + 72.0f * mxyIn * rhoIn * rhoVar + 108.0f * rhoIn * uxIn * rhoVar + 108.0f * rhoIn * uyIn * rhoVar - 139.0f * rhoVar * rhoVar);

        const dfloat mxyVar = mxyNumerator / (12.0f * rhoVar);

        // ---------------------------------------------------------
        // Etapa 6: myy
        // ---------------------------------------------------------
        const dfloat myyNumerator =
            -12.0f * mxxIn * rhoIn - 36.0f * mxyIn * rhoIn + 12.0f * myyIn * rhoIn -
            54.0f * rhoIn * uxIn - 54.0f * rhoIn * uyIn +
            7.0f * sqrtf(3.0f) * sqrtf(rhoVar * rhoVar) * (12.0f * rhoIn * rhoIn * uxIn * uxIn - 24.0f * rhoIn * rhoIn * uxIn * uyIn + 12.0f * rhoIn * rhoIn * uyIn * uyIn + 72.0f * mxyIn * rhoIn * rhoVar + 108.0f * rhoIn * uxIn * rhoVar + 108.0f * rhoIn * uyIn * rhoVar - 139.0f * rhoVar * rhoVar);

        const dfloat myyVar = myyNumerator / (143.0f * rhoVar);

        // ---------------------------------------------------------
        // Atribuir ao node
        // ---------------------------------------------------------
        (*node).rho = rhoVar;
        (*node).ux = uxVar;
        (*node).uy = uyVar;
        (*node).mxx = mxxVar;
        (*node).mxy = mxyVar;
        (*node).myy = myyVar;

        break;
    }
    case BB_NORTH_WEST:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[3] + pop[4] + pop[5] + pop[7] + pop[8];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat uxIn = ((pop[1] + pop[5] + pop[8]) - (pop[3] + pop[7])) * inv_rhoIn;
        const dfloat uyIn = ((pop[2] + pop[5]) - (pop[4] + pop[7] + pop[8])) * inv_rhoIn;

        const dfloat mxxIn = (pop[1] + pop[3] + pop[5] + pop[7] + pop[8]) * inv_rhoIn - cs2;
        const dfloat mxyIn = ((pop[5] + pop[7]) - pop[8]) * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[4] + pop[5] + pop[7] + pop[8]) * inv_rhoIn - cs2;

        // ---------------------------------------------------------
        // Etapa 1: rho
        // ---------------------------------------------------------
        const dfloat A = 9.0f * rhoIn - 9.0f * mxyIn * rhoIn - 9.0f * rhoIn * uxIn + 9.0f * rhoIn * uyIn;
        const dfloat B = 6.0f * rhoIn * rhoIn + 6.0f * rhoIn * rhoIn * uxIn + 3.0f * rhoIn * rhoIn * uxIn * uxIn - 6.0f * rhoIn * rhoIn * uyIn + 3.0f * rhoIn * rhoIn * uyIn * uyIn;

        const dfloat rhoVar = (1.0f / 28.0f) * (9.0f * rhoIn - 9.0f * mxyIn * rhoIn - 9.0f * rhoIn * uxIn + 9.0f * rhoIn * uyIn +
                                                sqrtf(A * A + 56.0f * B));

        // ---------------------------------------------------------
        // Etapa 2: ux
        // ---------------------------------------------------------
        const dfloat uxNumerator =
            6.0f * rhoIn * uxIn + 6.0f * rhoIn * uyIn + 21.0f * rhoVar -
            sqrtf(3.0f) * sqrtf(rhoVar * rhoVar) * (12.0f * rhoIn * rhoIn * uxIn * uxIn + 24.0f * rhoIn * rhoIn * uxIn * uyIn + 12.0f * rhoIn * rhoIn * uyIn * uyIn - 72.0f * mxyIn * rhoIn * rhoVar + 108.0f * rhoIn * uxIn * rhoVar - 108.0f * rhoIn * uyIn * rhoVar - 139.0f * rhoVar * rhoVar);

        const dfloat uxVar = uxNumerator / (12.0f * rhoVar);

        // ---------------------------------------------------------
        // Etapa 3: uy
        // ---------------------------------------------------------
        const dfloat uyVar = uxVar; // Idêntico ao valor anterior

        // ---------------------------------------------------------
        // Etapa 4: mxx
        // ---------------------------------------------------------
        const dfloat mxxNumerator =
            12.0f * mxxIn * rhoIn + 36.0f * mxyIn * rhoIn - 12.0f * myyIn * rhoIn +
            54.0f * rhoIn * uxIn - 54.0f * rhoIn * uyIn +
            7.0f * sqrtf(3.0f) * sqrtf(rhoVar * rhoVar) * (12.0f * rhoIn * rhoIn * uxIn * uxIn + 24.0f * rhoIn * rhoIn * uxIn * uyIn + 12.0f * rhoIn * rhoIn * uyIn * uyIn - 72.0f * mxyIn * rhoIn * rhoVar + 108.0f * rhoIn * uxIn * rhoVar - 108.0f * rhoIn * uyIn * rhoVar - 139.0f * rhoVar * rhoVar);

        const dfloat mxxVar = mxxNumerator / (143.0f * rhoVar);

        // ---------------------------------------------------------
        // Etapa 5: mxy
        // ---------------------------------------------------------
        const dfloat mxyNumerator =
            12.0f * mxyIn * rhoIn - 6.0f * rhoIn * uxIn + 6.0f * rhoIn * uyIn -
            21.0f * rhoVar +
            sqrtf(3.0f) * sqrtf(rhoVar * rhoVar) * (12.0f * rhoIn * rhoIn * uxIn * uxIn + 24.0f * rhoIn * rhoIn * uxIn * uyIn + 12.0f * rhoIn * rhoIn * uyIn * uyIn - 72.0f * mxyIn * rhoIn * rhoVar + 108.0f * rhoIn * uxIn * rhoVar - 108.0f * rhoIn * uyIn * rhoVar - 139.0f * rhoVar * rhoVar);

        const dfloat mxyVar = mxyNumerator / (12.0f * rhoVar);

        // ---------------------------------------------------------
        // Etapa 6: myy
        // ---------------------------------------------------------
        const dfloat myyNumerator =
            -12.0f * mxxIn * rhoIn + 36.0f * mxyIn * rhoIn + 12.0f * myyIn * rhoIn +
            54.0f * rhoIn * uxIn - 54.0f * rhoIn * uyIn +
            7.0f * sqrtf(3.0f) * sqrtf(rhoVar * rhoVar) * (12.0f * rhoIn * rhoIn * uxIn * uxIn + 24.0f * rhoIn * rhoIn * uxIn * uyIn + 12.0f * rhoIn * rhoIn * uyIn * uyIn - 72.0f * mxyIn * rhoIn * rhoVar + 108.0f * rhoIn * uxIn * rhoVar - 108.0f * rhoIn * uyIn * rhoVar - 139.0f * rhoVar * rhoVar);

        const dfloat myyVar = myyNumerator / (143.0f * rhoVar);

        // ---------------------------------------------------------
        // Atribuir ao node
        // ---------------------------------------------------------
        (*node).rho = rhoVar;
        (*node).ux = uxVar;
        (*node).uy = uyVar;
        (*node).mxx = mxxVar;
        (*node).mxy = mxyVar;
        (*node).myy = myyVar;

        break;
    }
    case BB_SOUTH_EAST:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[3] + pop[4] + pop[5] + pop[6] + pop[7];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat uxIn = ((pop[1] + pop[5]) - (pop[3] + pop[6] + pop[7])) * inv_rhoIn;
        const dfloat uyIn = ((pop[2] + pop[5] + pop[6]) - (pop[4] + pop[7])) * inv_rhoIn;

        const dfloat mxxIn = (pop[1] + pop[3] + pop[5] + pop[6] + pop[7]) * inv_rhoIn - cs2;
        const dfloat mxyIn = ((pop[5] + pop[7]) - pop[6]) * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[4] + pop[5] + pop[6] + pop[7]) * inv_rhoIn - cs2;

        // ---------------------------------------------------------
        // Etapa 1: rho
        // ---------------------------------------------------------
        const dfloat A = 9.0f * rhoIn - 9.0f * mxyIn * rhoIn - 9.0f * rhoIn * uxIn - 9.0f * rhoIn * uyIn;
        const dfloat B = 6.0f * rhoIn * rhoIn - 6.0f * rhoIn * rhoIn * uxIn + 3.0f * rhoIn * rhoIn * uxIn * uxIn + 6.0f * rhoIn * rhoIn * uyIn + 3.0f * rhoIn * rhoIn * uyIn * uyIn;

        const dfloat rhoVar = (1.0f / 28.0f) * (A + sqrtf(A * A + 56.0f * B));

        // ---------------------------------------------------------
        // Etapa 2: ux
        // ---------------------------------------------------------
        const dfloat uxNumerator =
            6.0f * rhoIn * uxIn + 6.0f * rhoIn * uyIn + 21.0f * rhoVar +
            sqrtf(3.0f) * sqrtf(rhoVar * rhoVar) * (-12.0f * rhoIn * rhoIn * uxIn * uxIn + 24.0f * rhoIn * rhoIn * uxIn * uyIn - 12.0f * rhoIn * rhoIn * uyIn * uyIn + 72.0f * mxyIn * rhoIn * rhoVar - 108.0f * rhoIn * uxIn * rhoVar - 108.0f * rhoIn * uyIn * rhoVar - 139.0f * rhoVar * rhoVar);

        const dfloat uxVar = uxNumerator / (12.0f * rhoVar);

        // ---------------------------------------------------------
        // Etapa 3: uy
        // ---------------------------------------------------------
        const dfloat uyNumerator =
            6.0f * rhoIn * uxIn + 6.0f * rhoIn * uyIn - 21.0f * rhoVar -
            sqrtf(3.0f) * sqrtf(rhoVar * rhoVar) * (12.0f * rhoIn * rhoIn * uxIn * uxIn - 24.0f * rhoIn * rhoIn * uxIn * uyIn + 12.0f * rhoIn * rhoIn * uyIn * uyIn + 72.0f * mxyIn * rhoIn * rhoVar - 108.0f * rhoIn * uxIn * rhoVar - 108.0f * rhoIn * uyIn * rhoVar - 139.0f * rhoVar * rhoVar);

        const dfloat uyVar = uyNumerator / (12.0f * rhoVar);

        // ---------------------------------------------------------
        // Etapa 4: mxx
        // ---------------------------------------------------------
        const dfloat mxxNumerator =
            12.0f * mxxIn * rhoIn + 36.0f * mxyIn * rhoIn - 12.0f * myyIn * rhoIn -
            54.0f * rhoIn * uxIn + 54.0f * rhoIn * uyIn +
            7.0f * sqrtf(3.0f) * sqrtf(rhoVar * rhoVar) * (12.0f * rhoIn * rhoIn * uxIn * uxIn - 24.0f * rhoIn * rhoIn * uxIn * uyIn + 12.0f * rhoIn * rhoIn * uyIn * uyIn + 72.0f * mxyIn * rhoIn * rhoVar - 108.0f * rhoIn * uxIn * rhoVar - 108.0f * rhoIn * uyIn * rhoVar - 139.0f * rhoVar * rhoVar);

        const dfloat mxxVar = mxxNumerator / (143.0f * rhoVar);

        // ---------------------------------------------------------
        // Etapa 5: mxy
        // ---------------------------------------------------------
        const dfloat mxyNumerator =
            12.0f * mxyIn * rhoIn + 6.0f * rhoIn * uxIn - 6.0f * rhoIn * uyIn - 21.0f * rhoVar -
            sqrtf(rhoVar * rhoVar) * sqrtf(3.0f) * (12.0f * rhoIn * rhoIn * uxIn * uxIn - 24.0f * rhoIn * rhoIn * uxIn * uyIn + 12.0f * rhoIn * rhoIn * uyIn * uyIn + 72.0f * mxyIn * rhoIn * rhoVar - 108.0f * rhoIn * uxIn * rhoVar - 108.0f * rhoIn * uyIn * rhoVar - 139.0f * rhoVar * rhoVar);

        const dfloat mxyVar = mxyNumerator / (12.0f * rhoVar);

        // ---------------------------------------------------------
        // Etapa 6: myy
        // ---------------------------------------------------------
        const dfloat myyNumerator =
            -12.0f * mxxIn * rhoIn + 36.0f * mxyIn * rhoIn + 12.0f * myyIn * rhoIn -
            54.0f * rhoIn * uxIn + 54.0f * rhoIn * uyIn +
            7.0f * sqrtf(3.0f) * sqrtf(rhoVar * rhoVar) * (12.0f * rhoIn * rhoIn * uxIn * uxIn - 24.0f * rhoIn * rhoIn * uxIn * uyIn + 12.0f * rhoIn * rhoIn * uyIn * uyIn + 72.0f * mxyIn * rhoIn * rhoVar - 108.0f * rhoIn * uxIn * rhoVar - 108.0f * rhoIn * uyIn * rhoVar - 139.0f * rhoVar * rhoVar);

        const dfloat myyVar = myyNumerator / (143.0f * rhoVar);

        // ---------------------------------------------------------
        // Atribuir ao node
        // ---------------------------------------------------------
        (*node).rho = rhoVar;
        (*node).ux = uxVar;
        (*node).uy = uyVar;
        (*node).mxx = mxxVar;
        (*node).mxy = mxyVar;
        (*node).myy = myyVar;

        break;
    }
    case BB_SOUTH_WEST:
    {
        const dfloat rhoIn = pop[0] + pop[1] + pop[2] + pop[3] + pop[4] + pop[5] + pop[6] + pop[8];
        const dfloat inv_rhoIn = 1.0f / rhoIn;

        const dfloat uxIn = ((pop[1] + pop[5] + pop[8]) - (pop[3] + pop[6])) * inv_rhoIn;
        const dfloat uyIn = ((pop[2] + pop[5] + pop[6]) - (pop[4] + pop[8])) * inv_rhoIn;

        const dfloat mxxIn = (pop[1] + pop[3] + pop[5] + pop[6] + pop[8]) * inv_rhoIn - cs2;
        const dfloat mxyIn = (pop[5] - (pop[6] + pop[8])) * inv_rhoIn;
        const dfloat myyIn = (pop[2] + pop[4] + pop[5] + pop[6] + pop[8]) * inv_rhoIn - cs2;

        // ---------------------------------------------------------
        // Etapa 1: rho
        // ---------------------------------------------------------
        const dfloat A = 9.0f * rhoIn + 9.0f * mxyIn * rhoIn - 9.0f * rhoIn * uxIn - 9.0f * rhoIn * uyIn;
        const dfloat B = 6.0f * rhoIn * rhoIn + 6.0f * rhoIn * rhoIn * uxIn + 3.0f * rhoIn * rhoIn * uxIn * uxIn + 6.0f * rhoIn * rhoIn * uyIn + 3.0f * rhoIn * rhoIn * uyIn * uyIn;

        const dfloat rhoVar = (1.0f / 28.0f) * (A + sqrtf(A * A + 56.0f * B));

        // ---------------------------------------------------------
        // Etapa 2: ux
        // ---------------------------------------------------------
        const dfloat uxNumerator =
            6.0f * rhoIn * uxIn - 6.0f * rhoIn * uyIn - 21.0f * rhoVar +
            sqrtf(3.0f) * sqrtf(rhoVar * rhoVar) * (12.0f * rhoIn * rhoIn * uxIn * uxIn - 24.0f * rhoIn * rhoIn * uxIn * uyIn + 12.0f * rhoIn * rhoIn * uyIn * uyIn + 72.0f * mxyIn * rhoIn * rhoVar - 108.0f * rhoIn * uxIn * rhoVar - 108.0f * rhoIn * uyIn * rhoVar - 139.0f * rhoVar * rhoVar);

        const dfloat uxVar = uxNumerator / (12.0f * rhoVar);

        // ---------------------------------------------------------
        // Etapa 3: uy
        // ---------------------------------------------------------
        const dfloat uyNumerator =
            -6.0f * rhoIn * uxIn + 6.0f * rhoIn * uyIn - 21.0f * rhoVar +
            sqrtf(3.0f) * sqrtf(rhoVar * rhoVar) * (12.0f * rhoIn * rhoIn * uxIn * uxIn - 24.0f * rhoIn * rhoIn * uxIn * uyIn + 12.0f * rhoIn * rhoIn * uyIn * uyIn + 72.0f * mxyIn * rhoIn * rhoVar - 108.0f * rhoIn * uxIn * rhoVar - 108.0f * rhoIn * uyIn * rhoVar - 139.0f * rhoVar * rhoVar);

        const dfloat uyVar = uyNumerator / (12.0f * rhoVar);

        // ---------------------------------------------------------
        // Etapa 4: mxx
        // ---------------------------------------------------------
        const dfloat mxxNumerator =
            12.0f * mxxIn * rhoIn - 36.0f * mxyIn * rhoIn - 12.0f * myyIn * rhoIn +
            54.0f * rhoIn * uxIn + 54.0f * rhoIn * uyIn +
            7.0f * sqrtf(3.0f) * sqrtf(rhoVar * rhoVar) * (12.0f * rhoIn * rhoIn * uxIn * uxIn - 24.0f * rhoIn * rhoIn * uxIn * uyIn + 12.0f * rhoIn * rhoIn * uyIn * uyIn + 72.0f * mxyIn * rhoIn * rhoVar - 108.0f * rhoIn * uxIn * rhoVar - 108.0f * rhoIn * uyIn * rhoVar - 139.0f * rhoVar * rhoVar);

        const dfloat mxxVar = mxxNumerator / (143.0f * rhoVar);

        // ---------------------------------------------------------
        // Etapa 5: mxy
        // ---------------------------------------------------------
        const dfloat mxyNumerator =
            12.0f * mxyIn * rhoIn + 6.0f * rhoIn * uxIn + 6.0f * rhoIn * uyIn + 21.0f * rhoVar -
            sqrtf(3.0f) * sqrtf(rhoVar * rhoVar) * (12.0f * rhoIn * rhoIn * uxIn * uxIn - 24.0f * rhoIn * rhoIn * uxIn * uyIn + 12.0f * rhoIn * rhoIn * uyIn * uyIn + 72.0f * mxyIn * rhoIn * rhoVar - 108.0f * rhoIn * uxIn * rhoVar - 108.0f * rhoIn * uyIn * rhoVar - 139.0f * rhoVar * rhoVar);

        const dfloat mxyVar = mxyNumerator / (12.0f * rhoVar);

        // ---------------------------------------------------------
        // Etapa 6: myy
        // ---------------------------------------------------------
        const dfloat myyNumerator =
            -12.0f * mxxIn * rhoIn - 36.0f * mxyIn * rhoIn + 12.0f * myyIn * rhoIn +
            54.0f * rhoIn * uxIn + 54.0f * rhoIn * uyIn +
            7.0f * sqrtf(3.0f) * sqrtf(rhoVar * rhoVar) * (12.0f * rhoIn * rhoIn * uxIn * uxIn - 24.0f * rhoIn * rhoIn * uxIn * uyIn + 12.0f * rhoIn * rhoIn * uyIn * uyIn + 72.0f * mxyIn * rhoIn * rhoVar - 108.0f * rhoIn * uxIn * rhoVar - 108.0f * rhoIn * uyIn * rhoVar - 139.0f * rhoVar * rhoVar);

        const dfloat myyVar = myyNumerator / (143.0f * rhoVar);

        // ---------------------------------------------------------
        // Atribuir ao node
        // ---------------------------------------------------------
        (*node).rho = rhoVar;
        (*node).ux = uxVar;
        (*node).uy = uyVar;
        (*node).mxx = mxxVar;
        (*node).mxy = mxyVar;
        (*node).myy = myyVar;

        break;
    }
    default:
        break;
    }
}

__host__ inline void streaming(latticeNode *nodes, size_t x_lattices, size_t y_lattices)
{
    for (size_t y = 0; y < y_lattices; ++y)
    {
        for (size_t x = 0; x < x_lattices; ++x)
        {
            size_t xp1 = (x + 1 + x_lattices) % x_lattices;
            size_t xm1 = (x - 1 + x_lattices) % x_lattices;
            size_t yp1 = (y + 1 + y_lattices) % y_lattices;
            size_t ym1 = (y - 1 + y_lattices) % y_lattices;

            nodes[fine_idx(xp1, y)].pop_in[1] = nodes[fine_idx(x, y)].pop_out[1];
            nodes[fine_idx(x, yp1)].pop_in[2] = nodes[fine_idx(x, y)].pop_out[2];
            nodes[fine_idx(xm1, y)].pop_in[3] = nodes[fine_idx(x, y)].pop_out[3];
            nodes[fine_idx(x, ym1)].pop_in[4] = nodes[fine_idx(x, y)].pop_out[4];
            nodes[fine_idx(xp1, yp1)].pop_in[5] = nodes[fine_idx(x, y)].pop_out[5];
            nodes[fine_idx(xm1, yp1)].pop_in[6] = nodes[fine_idx(x, y)].pop_out[6];
            nodes[fine_idx(xm1, ym1)].pop_in[7] = nodes[fine_idx(x, y)].pop_out[7];
            nodes[fine_idx(xp1, ym1)].pop_in[8] = nodes[fine_idx(x, y)].pop_out[8];
        }
    }
}

__host__ inline void collision(latticeNode *node, dfloat omega)
{
    const dfloat omegaVar = omega;
    const dfloat t_omegaVar = 1 - omegaVar;
    const dfloat omegaVar_d2 = omegaVar / 2;

    (*node).mxx = (t_omegaVar * (*node).mxx + omegaVar_d2 * (*node).ux * (*node).ux);
    (*node).myy = (t_omegaVar * (*node).myy + omegaVar_d2 * (*node).uy * (*node).uy);

    (*node).mxy = (t_omegaVar * (*node).mxy + omegaVar * (*node).ux * (*node).uy);
}

#endif