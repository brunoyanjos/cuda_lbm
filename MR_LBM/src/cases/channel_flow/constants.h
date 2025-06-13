#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "../../var.h"

constexpr dfloat RE = 100;

constexpr int SCALE = 1;
constexpr int N_STEPS = 100000;

constexpr int MACR_SAVE = 1000;

constexpr int N = 64 * SCALE;
constexpr int NX = 10 * N;               // size x of the grid
constexpr int NY = N;               // size y of the grid

constexpr dfloat U_MAX = 0.1;
constexpr dfloat L = N;

constexpr dfloat VISC = U_MAX * (NY-1) / RE;
constexpr dfloat TAU = 0.5 + 3.0 * VISC; // relaxation time

constexpr dfloat OMEGA = 1.0 / TAU;            // (tau)^-1
constexpr dfloat OMEGAd2 = OMEGA / 2.0;        // OMEGA/2
constexpr dfloat OMEGAd9 = OMEGA / 9.0;        // OMEGA/9
constexpr dfloat T_OMEGA = 1.0 - OMEGA;        // 1-OMEGA
constexpr dfloat TT_OMEGA = 1.0 - 0.5 * OMEGA; // 1.0 - OMEGA/2
constexpr dfloat OMEGA_P1 = 1.0 + OMEGA;       // 1+ OMEGA
constexpr dfloat TT_OMEGA_T3 = TT_OMEGA * 3.0; // 3*(1-0.5*OMEGA)

// value for the velocity initial condition in the domain
constexpr dfloat U_0_X = 0.0;
constexpr dfloat U_0_Y = 0.0;
constexpr dfloat U_0_Z = 0.0;
constexpr dfloat RHO_0 = 1.0;

constexpr dfloat MACH_NUMBER = U_MAX / 0.57735026918962;

/* --------------------- INITIALIZATION LOADING DEFINES -------------------- */
constexpr int INI_STEP = 0; // initial simulation step (0 default)

#define BC_X_WALL
#define BC_Y_WALL
// #define BC_Y_PERIODIC

constexpr bool IRBC = false;

#endif // !CONSTANTS_H