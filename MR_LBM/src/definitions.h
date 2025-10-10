#ifndef __DEFINITIONS_H
#define __DEFINITIONS_H

#include "var.h"

/* ------------------------------ VELOCITY SET ------------------------------ */
constexpr unsigned char Q = 9;	// number of velocities
constexpr unsigned char QF = 3; // number of velocities on each face
constexpr dfloat W0 = 4.0 / 9;	// population 0 weight (0, 0, 0)
constexpr dfloat W1 = 1.0 / 9;	// adjacent populations (1, 0, 0)
constexpr dfloat W2 = 1.0 / 36; // diagonal populations (1, 1, 0)

// velocities weight vector
__device__ const dfloat w[Q] = {W0,
								W1, W1, W1, W1,
								W2, W2, W2, W2};

constexpr dfloat as2 = 3.0;
constexpr dfloat cs2 = 1.0 / as2;

// populations velocities      0  1  2  3  4  5  6  7  8
__device__ constexpr dfloat cx[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
__device__ constexpr dfloat cy[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};

constexpr dfloat F_M_0_SCALE = 1.0;
constexpr dfloat F_M_I_SCALE = as2;
constexpr dfloat F_M_II_SCALE = as2 * as2 / 2;
constexpr dfloat F_M_IJ_SCALE = as2 * as2;

/* ------------------------------ FINE SIZE ------------------------------ */
#include "arrayIndex.h"

/* ------------------------------ MEMORY SIZE ------------------------------ */

const size_t NUMBER_LBM_NODES = NUMBER_OF_COARSE_NODES * NUMBER_OF_FINE_NODES;

const size_t MEM_SIZE_COARSE_NODES = sizeof(dfloat) * NUMBER_OF_COARSE_NODES;
const size_t MEM_SIZE_FINE_NODES = sizeof(dfloat) * NUMBER_OF_FINE_NODES;

#endif //!__DEFINITIONS_H