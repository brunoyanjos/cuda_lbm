#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "../../var.h"

constexpr dfloat RE = 1000;

constexpr int SCALE = 1;

constexpr int MACR_SAVE = 20000;

constexpr int N_COARSE = 3 * SCALE;
constexpr int NX = N_COARSE; // size x of the grid
constexpr int NY = N_COARSE; // size y of the grid

constexpr int N_OVERLAP_LAYER = 1;
constexpr int N_EXTRA_LAYER = 1;

constexpr int GRID_RATIO = 2;

constexpr int NX_TOTAL_SIZE = NX + 2 * N_OVERLAP_LAYER + 2 * N_EXTRA_LAYER;
constexpr int NY_TOTAL_SIZE = NY + 2 * N_OVERLAP_LAYER + 2 * N_EXTRA_LAYER;

constexpr int NX_FINE_GRID = (NX_TOTAL_SIZE - 1) * GRID_RATIO + 1;
constexpr int NY_FINE_GRID = (NY_TOTAL_SIZE - 1) * GRID_RATIO + 1;

constexpr int NX_COARSE_GRID = NX + 2 * N_OVERLAP_LAYER;
constexpr int NY_COARSE_GRID = NY + 2 * N_OVERLAP_LAYER;

constexpr int FINE_WIDTH = (N_OVERLAP_LAYER + N_EXTRA_LAYER) * GRID_RATIO + 1;

constexpr int NUMBER_OF_COARSE_NODES = NX_COARSE_GRID * NY_COARSE_GRID;
constexpr int NUMBER_OF_FINE_NODES = NX_FINE_GRID * NY_FINE_GRID;
// constexpr int NUMBER_OF_FINE_NODES = 2 * (FINE_WIDTH * NX_FINE_GRID + FINE_WIDTH * (NY_FINE_GRID - 2 * FINE_WIDTH));

constexpr dfloat U_MAX = 0.0256;
constexpr dfloat L = N_COARSE;

constexpr dfloat VISC_FINE = U_MAX * (NX_FINE_GRID - 1) / RE;
constexpr dfloat VISC_COARSE = U_MAX * (NX_COARSE_GRID - 1) / RE;
constexpr dfloat TAU_FINE = 0.5 + 3.0 * VISC_FINE; // relaxation time
constexpr dfloat TAU_COARSE = 0.5 + 3.0 * VISC_COARSE; // relaxation time

constexpr dfloat OMEGA = 1.0;

constexpr dfloat OMEGA_FINE = 1.0 / TAU_FINE; // (tau)^-1
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
constexpr int N_STEPS = 3;

#define BC_X_WALL
#define BC_Y_WALL

constexpr bool IRBC = true;

#endif // !CONSTANTS_H
