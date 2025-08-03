#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "../../var.h"

constexpr dfloat RE = 1000;

constexpr int SCALE = 1;

constexpr int MACR_SAVE = 20000;

constexpr int N = 5 * SCALE;
constexpr int NX = N; // size x of the grid
constexpr int NY = N; // size y of the grid

constexpr int N_OVERLAP_LAYER = 1;

constexpr int GRID_RATIO = 2;

constexpr dfloat U_MAX = 0.0256;

constexpr size_t NX_COARSE = static_cast<int>(N / 2) + 1;
constexpr size_t NY_COARSE = NY;

constexpr size_t NX_FINE = NX;
constexpr size_t NY_FINE = NY * GRID_RATIO - 1;

constexpr size_t NUMBER_OF_COARSE_NODES = (NX_COARSE + N_OVERLAP_LAYER) * NY_COARSE;
constexpr size_t NUMBER_OF_FINE_NODES = NX_FINE * NY_FINE;

constexpr dfloat VISC_FINE = U_MAX * (N * 2 - 1) / RE;
constexpr dfloat VISC_COARSE = U_MAX * (NY_COARSE - 1) / RE;
constexpr dfloat TAU_FINE = 0.5 + 3.0 * VISC_FINE;     // relaxation time
constexpr dfloat TAU_COARSE = 0.5 + 3.0 * VISC_COARSE; // relaxation time

constexpr dfloat OMEGA = 1.0;

constexpr dfloat OMEGA_FINE = 1.0 / TAU_FINE;     // (tau)^-1
constexpr dfloat OMEGA_COARSE = 1.0 / TAU_COARSE; // (tau)^-1

// value for the velocity initial condition in the domain
constexpr dfloat U_0_X = 0.0;
constexpr dfloat U_0_Y = 0.0;
constexpr dfloat U_0_Z = 0.0;
constexpr dfloat RHO_0 = 1.0;

constexpr dfloat MACH_NUMBER = U_MAX / 0.57735026918962;

/* --------------------- INITIALIZATION LOADING DEFINES -------------------- */
constexpr int INI_STEP = 0; // initial simulation step (0 default)

// constexpr int T_STAR_FINAL = 1200;
constexpr int N_STEPS = 100;

#define BC_X_WALL
#define BC_Y_WALL

constexpr bool IRBC = true;

#endif // !CONSTANTS_H
