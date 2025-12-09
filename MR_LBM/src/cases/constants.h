#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "../var.h"

constexpr dfloat RE = 100;
constexpr int D = 32;

constexpr int NX = 2 * D;
constexpr int NY = NX;

constexpr dfloat xc = dfloat(NX - 1) / 2;
constexpr dfloat yc = dfloat(NY - 1) / 2;

constexpr dfloat U_MAX = 0.1;

constexpr int MACR_SAVE = 10 * D / U_MAX;
constexpr int tstar = 100;
constexpr int stat_period = 100;

constexpr int N_STEPS = (tstar + stat_period) * D / U_MAX;

constexpr dfloat RHO_0 = 1.0;

constexpr dfloat MACH_NUMBER = U_MAX / 0.57735026918962;

/* --------------------- INITIALIZATION LOADING DEFINES -------------------- */
constexpr int INI_STEP = 0;

#define BC_X_WALL
#define BC_Y_PERIODIC

#endif // !CONSTANTS_H