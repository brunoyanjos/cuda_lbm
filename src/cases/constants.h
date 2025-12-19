#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "../var.h"

constexpr dfloat RE = 100;

constexpr int SCALE = 1;
constexpr int D = 128;

constexpr int N = 2 * D;
constexpr int NX = N; // size x of the grid
constexpr int NY = N; // size y of the grid

constexpr dfloat xc = dfloat(NX - 1) / 2;
constexpr dfloat yc = dfloat(NY - 1) / 2;

constexpr dfloat U_MAX = 0.0256;
constexpr dfloat L = N;

// value for the velocity initial condition in the domain

constexpr dfloat U_0_X = 0.0;
constexpr dfloat U_0_Y = 0.0;
constexpr dfloat U_0_Z = 0.0;
constexpr dfloat RHO_0 = 1.0;

constexpr dfloat MACH_NUMBER = U_MAX / 0.57735026918962;

constexpr dfloat R2 = N / 2;
constexpr dfloat R1 = D / 2;

constexpr dfloat VISC = U_MAX * (R2 - R1) / RE;
constexpr dfloat TAU = 0.5 + 3.0 * VISC;

constexpr dfloat OMEGA = 1.0 / TAU;

/* --------------------- INITIALIZATION LOADING DEFINES -------------------- */
constexpr int INI_STEP = 0; // initial simulation step (0 default)

constexpr int t_star = 2000;

constexpr int N_STEPS = t_star * (N - D) / U_MAX;
constexpr int MACR_SAVE = (N - D) / U_MAX;

#define BC_X_WALL
#define BC_Y_WALL

#endif // !CONSTANTS_H
