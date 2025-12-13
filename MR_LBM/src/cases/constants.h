#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "../var.h"

constexpr dfloat RE = 1000;
constexpr int D = 256;

constexpr int NX = 64;
constexpr int NY = NX;

constexpr dfloat xc = dfloat(NX - 1) / 2;
constexpr dfloat yc = dfloat(NY - 1) / 2;

constexpr dfloat U_MAX = 0.0256;

constexpr int MACR_SAVE = 10000;
constexpr int tstar = 100;
constexpr int stat_period = 100;

constexpr int N_STEPS = 100000;

constexpr dfloat RHO_0 = 1.0;

constexpr dfloat MACH_NUMBER = U_MAX / 0.57735026918962;

/* --------------------- INITIALIZATION LOADING DEFINES -------------------- */
constexpr int INI_STEP = 0;

#define BC_X_WALL
#define BC_Y_WALL

#endif // !CONSTANTS_H