#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "../var.h"

constexpr dfloat RE = 100;

constexpr int SCALE = 1;

constexpr int MACR_SAVE = 1000;
constexpr int D = 128;

constexpr int N = 2 * D;
constexpr int NX = N; // size x of the grid
constexpr int NY = N; // size y of the grid

constexpr dfloat xc = dfloat(NX - 1) / 2;
constexpr dfloat yc = dfloat(NY - 1) / 2;

constexpr dfloat U_IN = 0.0256;
constexpr dfloat U_OUT = 0.0;
constexpr dfloat L = N;

// value for the velocity initial condition in the domain
constexpr dfloat U_0_X = 0.0;
constexpr dfloat U_0_Y = 0.0;
constexpr dfloat U_0_Z = 0.0;
constexpr dfloat RHO_0 = 1.0;

/* --------------------- INITIALIZATION LOADING DEFINES -------------------- */
constexpr int INI_STEP = 0; // initial simulation step (0 default)

constexpr int N_STEPS = 1000000;

#define BC_X_WALL
#define BC_Y_WALL

#endif // !CONSTANTS_H
