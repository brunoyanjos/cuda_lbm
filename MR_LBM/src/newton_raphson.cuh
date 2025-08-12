#ifndef NEWTON_RAPHSON_CUH
#define NEWTON_RAPHSON_CUH

#include "var.h"

__host__ inline dfloat f1(const dfloat rhoIn, const dfloat rhoUxIn, const dfloat omega, const dfloat rho, const dfloat rhoMxx)
{
    const dfloat rhoUx = static_cast<dfloat>(6) * rhoUxIn + rho + static_cast<dfloat>(3) * rhoMxx;
    const dfloat incoming_terms = static_cast<dfloat>(18) * rhoIn - static_cast<dfloat>(18) * rhoUxIn;

    return incoming_terms * rho + omega * rhoUx * rhoUx - static_cast<dfloat>(18) * rho * rho - static_cast<dfloat>(9) * omega * rho * rhoMxx;
}

__host__ inline dfloat df1_drho(const dfloat rhoIn, const dfloat rhoUxIn, const dfloat omega, const dfloat rho, const dfloat rhoMxx)
{
    const dfloat rhoUx = static_cast<dfloat>(6) * rhoUxIn + rho + static_cast<dfloat>(3) * rhoMxx;
    const dfloat incoming_terms = static_cast<dfloat>(18) * rhoIn - static_cast<dfloat>(18) * rhoUxIn;

    return incoming_terms + static_cast<dfloat>(2) * omega * rhoUx - static_cast<dfloat>(36) * rho - static_cast<dfloat>(9) * omega * rhoMxx;
}

__host__ inline dfloat df1_dmxx(const dfloat rhoIn, const dfloat rhoUxIn, const dfloat omega, const dfloat rho, const dfloat rhoMxx)
{
    const dfloat rhoUx = static_cast<dfloat>(6) * rhoUxIn + rho + static_cast<dfloat>(3) * rhoMxx;

    return static_cast<dfloat>(6) * omega * rhoUx - static_cast<dfloat>(9) * omega * rho;
}

__host__ inline dfloat f2(const dfloat rhoUxIn, const dfloat rho, const dfloat rhoUy, const dfloat rhoMxx, const dfloat rhoMyy)
{
    const dfloat rhoUx = static_cast<dfloat>(6) * rhoUxIn + rho + static_cast<dfloat>(3) * rhoMxx;

    return rho * (rhoMxx + rhoMyy) - rhoUx * rhoUx / static_cast<dfloat>(9) - rhoUy * rhoUy;
}

__host__ inline dfloat df2_drho(const dfloat rhoUxIn, const dfloat rho, const dfloat rhoUy, const dfloat rhoMxx, const dfloat rhoMyy)
{
    const dfloat rhoUx = static_cast<dfloat>(6) * rhoUxIn + rho + static_cast<dfloat>(3) * rhoMxx;

    return rhoMxx + rhoMyy - static_cast<dfloat>(2) * rhoUx / static_cast<dfloat>(9);
}

__host__ inline dfloat df2_dmxx(const dfloat rhoUxIn, const dfloat rho, const dfloat rhoUy, const dfloat rhoMxx, const dfloat rhoMyy)
{
    const dfloat rhoUx = static_cast<dfloat>(6) * rhoUxIn + rho + static_cast<dfloat>(3) * rhoMxx;

    return rho - static_cast<dfloat>(2) * rhoUx / static_cast<dfloat>(3);
}

__host__ inline void newton_raphson(const dfloat rhoIn, const dfloat rhoUxIn, const dfloat omega, dfloat *rho, const dfloat rhoUy, dfloat *rhoMxx, const dfloat rhoMyy)
{
    dfloat rho_error = 100;
    dfloat rho_mxx_error = 100;
    const dfloat stop_error = 0.1;

    dfloat rho_guess = 1.0;
    dfloat rho_mxx_guess = 0.0;

    while (rho_error > stop_error && rho_mxx_error > stop_error)
    {
        const dfloat a1 = df1_drho(rhoIn, rhoUxIn, omega, rho_guess, rho_mxx_guess);
        const dfloat a2 = df1_dmxx(rhoIn, rhoUxIn, omega, rho_guess, rho_mxx_guess);

        const dfloat a3 = df2_drho(rhoUxIn, rho_guess, rhoUy, rho_mxx_guess, rhoMyy);
        const dfloat a4 = df2_dmxx(rhoUxIn, rho_guess, rhoUy, rho_mxx_guess, rhoMyy);

        const dfloat b1 = rho_guess * a1 + rho_mxx_guess * a2 - f1(rhoIn, rhoUxIn, omega, rho_guess, rho_mxx_guess);
        const dfloat b2 = rho_guess * a3 + rho_mxx_guess * a4 - f2(rhoUxIn, rho_guess, rhoUy, rho_mxx_guess, rhoMyy);

        const dfloat calculated_rho_mxx = (b1 * a3 - b2 * a1) / (a3 * a2 - a1 * a4);
        const dfloat calculated_rho = (b1 - a2 * calculated_rho_mxx) / a1;

        rho_error = std::abs(calculated_rho - rho_guess) / std::abs(calculated_rho) * 100;
        rho_mxx_error = std::abs(calculated_rho_mxx - rho_mxx_guess) / std::abs(calculated_rho_mxx) * 100;

        rho_guess = calculated_rho;
        rho_mxx_guess = calculated_rho_mxx;
    }

    *rho = rho_guess;
    *rhoMxx = rho_mxx_guess;
}

#endif