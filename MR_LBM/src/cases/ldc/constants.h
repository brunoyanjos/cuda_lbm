#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "../../var.h"

constexpr dfloat RE = 1000;

constexpr int SCALE = 1;

constexpr int MACR_SAVE = 20000;

constexpr int N = 128 * SCALE;
constexpr int NX = N; // size x of the grid
constexpr int NY = N; // size y of the grid

constexpr dfloat U_MAX = 0.0256;
constexpr dfloat L = N;

constexpr dfloat VISC = U_MAX * NX / RE;
constexpr dfloat TAU = 0.5 + 3.0 * VISC; // relaxation time

constexpr dfloat OMEGA = 1.0 / TAU; // (tau)^-1

// value for the velocity initial condition in the domain
constexpr dfloat U_0_X = 0.0;
constexpr dfloat U_0_Y = 0.0;
constexpr dfloat U_0_Z = 0.0;
constexpr dfloat RHO_0 = 1.0;

constexpr dfloat MACH_NUMBER = U_MAX / 0.57735026918962;

/* --------------------- INITIALIZATION LOADING DEFINES -------------------- */
constexpr int INI_STEP = 0; // initial simulation step (0 default)

constexpr int T_STAR_FINAL = 1200;
constexpr int N_STEPS = T_STAR_FINAL * NX / U_MAX;

#define BC_X_WALL
#define BC_Y_WALL

constexpr bool IRBC = true;

#endif // !CONSTANTS_H
