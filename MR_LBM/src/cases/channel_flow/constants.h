#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "../../var.h"

constexpr dfloat RE = 100;

constexpr int SCALE = 1;

constexpr int N = 31 * SCALE;
constexpr int NX = 4 * N + 1;
constexpr int NY = N;

constexpr int N_OVERLAP_LAYER = 1;

constexpr int GRID_RATIO = 2;
constexpr dfloat U_MAX = 0.01;

constexpr size_t NX_COARSE = static_cast<size_t>(NX / 2) + 1 + N_OVERLAP_LAYER;
constexpr size_t NY_COARSE = NY;

constexpr size_t NX_FINE = NX;
constexpr size_t NY_FINE = NY * GRID_RATIO - 1;

constexpr size_t NUMBER_OF_COARSE_NODES = NX_COARSE * NY_COARSE;
constexpr size_t NUMBER_OF_FINE_NODES = NX_FINE * NY_FINE;

constexpr dfloat VISC_FINE = U_MAX * (NY_FINE - 1) / RE;
constexpr dfloat VISC_COARSE = U_MAX * (NY_COARSE - 1) / RE;

constexpr dfloat TAU_FINE = 0.5 + 3.0 * VISC_FINE;
constexpr dfloat TAU_COARSE = 0.5 + 3.0 * VISC_COARSE;

constexpr dfloat OMEGA_FINE = 1.0 / TAU_FINE;
constexpr dfloat OMEGA_COARSE = 1.0 / TAU_COARSE;

constexpr dfloat ALPHA = GRID_RATIO * OMEGA_FINE / OMEGA_COARSE;
constexpr dfloat INV_ALPHA = static_cast<dfloat>(1) / ALPHA;

// value for the velocity initial condition in the domain
constexpr dfloat U_0_X = 0.0;
constexpr dfloat U_0_Y = 0.0;
constexpr dfloat U_0_Z = 0.0;
constexpr dfloat RHO_0 = 1.0;

constexpr dfloat MACH_NUMBER = U_MAX / 0.57735026918962;

/* --------------------- INITIALIZATION LOADING DEFINES -------------------- */
constexpr int N_STEPS = 1;
constexpr int MACR_SAVE = 100;

#define BC_X_WALL
#define BC_Y_WALL

constexpr bool IRBC = true;

#endif // !CONSTANTS_H