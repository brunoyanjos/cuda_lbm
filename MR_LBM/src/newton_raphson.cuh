#ifndef NEWTON_RAPHSON_CUH
#define NEWTON_RAPHSON_CUH

#include "var.h"

__host__ inline dfloat f1(const dfloat rhoIn, const dfloat rhoUyIn, const dfloat omega, const dfloat rho, const dfloat rhoMyy)
{
    const dfloat rhoUy = rho + static_cast<dfloat>(3) * rhoMyy - static_cast<dfloat>(6) * rhoUyIn;
    const dfloat incoming_terms = static_cast<dfloat>(18) * rhoIn + static_cast<dfloat>(18) * rhoUyIn;

    return incoming_terms * rho + omega * rhoUy * rhoUy - static_cast<dfloat>(18) * rho * rho - static_cast<dfloat>(9) * omega * rho * rhoMyy;
}

__host__ inline dfloat df1_drho(const dfloat rhoIn, const dfloat rhoUyIn, const dfloat omega, const dfloat rho, const dfloat rhoMyy)
{
    const dfloat rhoUy = rho + static_cast<dfloat>(3) * rhoMyy - static_cast<dfloat>(6) * rhoUyIn;
    const dfloat incoming_terms = static_cast<dfloat>(18) * rhoIn + static_cast<dfloat>(18) * rhoUyIn;

    return incoming_terms + static_cast<dfloat>(2) * omega * rhoUy - static_cast<dfloat>(36) * rho - static_cast<dfloat>(9) * omega * rhoMyy;
}

__host__ inline dfloat df1_dmyy(const dfloat rhoIn, const dfloat rhoUyIn, const dfloat omega, const dfloat rho, const dfloat rhoMyy)
{
    const dfloat rhoUy = rho + static_cast<dfloat>(3) * rhoMyy - static_cast<dfloat>(6) * rhoUyIn;

    return static_cast<dfloat>(6) * omega * rhoUy - static_cast<dfloat>(9) * omega * rho;
}

__host__ inline dfloat f2(const dfloat rhoUyIn, const dfloat rho, const dfloat rhoUx, const dfloat rhoMxx, const dfloat rhoMyy)
{
    const dfloat rhoUy = rho + static_cast<dfloat>(3) * rhoMyy - static_cast<dfloat>(6) * rhoUyIn;

    return rho * (rhoMxx + rhoMyy) - rhoUy * rhoUy / static_cast<dfloat>(9) - rhoUx * rhoUx;
}

__host__ inline dfloat df2_drho(const dfloat rhoUyIn, const dfloat rho, const dfloat rhoMxx, const dfloat rhoMyy)
{
    const dfloat rhoUy = rho + static_cast<dfloat>(3) * rhoMyy - static_cast<dfloat>(6) * rhoUyIn;

    return rhoMxx + rhoMyy - static_cast<dfloat>(2) * rhoUy / static_cast<dfloat>(9);
}

__host__ inline dfloat df2_dmyy(const dfloat rhoUyIn, const dfloat rho, const dfloat rhoMyy)
{
    const dfloat rhoUy = rho + static_cast<dfloat>(3) * rhoMyy - static_cast<dfloat>(6) * rhoUyIn;

    return rho - static_cast<dfloat>(2) * rhoUy / static_cast<dfloat>(3);
}

__host__ inline void newton_raphson(const dfloat rhoIn, const dfloat rhoUyIn, const dfloat omega, dfloat *const __restrict__ rho, const dfloat rhoUx, const dfloat rhoMxx, dfloat *const __restrict__ rhoMyy)
{
    dfloat rho_error = 100;
    dfloat rho_mxx_error = 100;
    const dfloat stop_error = 0.1;

    dfloat rho_guess = 1.0;
    dfloat rho_myy_guess = 0.0;

    while (rho_error > stop_error || rho_mxx_error > stop_error)
    {
        const dfloat a1 = df1_drho(rhoIn, rhoUyIn, omega, rho_guess, rho_myy_guess);
        const dfloat a2 = df1_dmyy(rhoIn, rhoUyIn, omega, rho_guess, rho_myy_guess);

        const dfloat a3 = df2_drho(rhoUyIn, rho_guess, rhoMxx, rho_myy_guess);
        const dfloat a4 = df2_dmyy(rhoUyIn, rho_guess, rho_myy_guess);

        const dfloat b1 = rho_guess * a1 + rho_myy_guess * a2 - f1(rhoIn, rhoUyIn, omega, rho_guess, rho_myy_guess);
        const dfloat b2 = rho_guess * a3 + rho_myy_guess * a4 - f2(rhoUyIn, rho_guess, rhoUx, rhoMxx, rho_myy_guess);

        const dfloat calculated_rho_myy = (b1 * a3 - b2 * a1) / (a3 * a2 - a1 * a4);
        const dfloat calculated_rho = (b1 - a2 * calculated_rho_myy) / a1;

        rho_error = std::abs(calculated_rho - rho_guess) / std::abs(calculated_rho) * 100;
        rho_mxx_error = std::abs(calculated_rho_myy - rho_myy_guess) / std::abs(calculated_rho_myy) * 100;

        rho_guess = calculated_rho;
        rho_myy_guess = calculated_rho_myy;
    }

    *rho = rho_guess;
    *rhoMyy = rho_myy_guess;
}

#endif