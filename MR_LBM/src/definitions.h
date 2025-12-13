#ifndef __DEFINITIONS_H
#define __DEFINITIONS_H

#include "var.h"

/* ------------------------------ VELOCITY SET ------------------------------ */
constexpr unsigned char Q = 9;
constexpr unsigned char QF = 3;
constexpr dfloat W0 = 4.0 / 9;
constexpr dfloat W1 = 1.0 / 9;
constexpr dfloat W2 = 1.0 / 36;

// velocities weight vector
__device__ const dfloat w[Q] = {W0,
								W1, W1, W1, W1,
								W2, W2, W2, W2};

constexpr dfloat as2 = 3.0;
constexpr dfloat as4 = 9.0;
constexpr dfloat cs2 = 1.0 / as2;

// populations velocities      0  1  2  3  4  5  6  7  8
__device__ constexpr dfloat cx[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
__device__ constexpr dfloat cy[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};

constexpr dfloat F_M_0_SCALE = 1.0;
constexpr dfloat F_M_I_SCALE = as2;
constexpr dfloat F_M_II_SCALE = as2 * as2 / 2;
constexpr dfloat F_M_IJ_SCALE = as2 * as2;

/* ------------------------------ MEMORY SIZE ------------------------------ */
#include "arrayIndex.h"

constexpr int SHARED_MEMORY_ELEMENT_SIZE = sizeof(dfloat) * (Q - 1);
constexpr int MAX_ELEMENTS_IN_BLOCK = 48128 / SHARED_MEMORY_ELEMENT_SIZE;

constexpr BlockDim optimalBlockDimArray = findOptimalBlockDimensions(MAX_ELEMENTS_IN_BLOCK);

constexpr int BLOCK_NX = optimalBlockDimArray.x;
constexpr int BLOCK_NY = optimalBlockDimArray.y / 2;

#define BLOCK_LBM_SIZE (BLOCK_NX * BLOCK_NY)

constexpr size_t BLOCK_GHOST_SIZE = BLOCK_NX + BLOCK_NY;

constexpr size_t BLOCK_SIZE = BLOCK_LBM_SIZE + BLOCK_GHOST_SIZE;

constexpr size_t NUM_BLOCK_X = NX / BLOCK_NX;
constexpr size_t NUM_BLOCK_Y = NY / BLOCK_NY;

constexpr size_t NUM_BLOCK = NUM_BLOCK_X * NUM_BLOCK_Y;

constexpr size_t NUMBER_LBM_NODES = NUM_BLOCK * BLOCK_LBM_SIZE;
constexpr size_t NUMBER_GHOST_FACE_X = BLOCK_NY * NUM_BLOCK_X * NUM_BLOCK_Y;
constexpr size_t NUMBER_GHOST_FACE_Y = BLOCK_NX * NUM_BLOCK_X * NUM_BLOCK_Y;

constexpr size_t MEM_SIZE_BLOCK_LBM = sizeof(dfloat) * BLOCK_LBM_SIZE * NUMBER_MOMENTS;
constexpr size_t MEM_SIZE_BLOCK_GHOST = sizeof(dfloat) * BLOCK_GHOST_SIZE * Q;
constexpr size_t MEM_SIZE_BLOCK_TOTAL = MEM_SIZE_BLOCK_GHOST + MEM_SIZE_BLOCK_LBM;

constexpr size_t NUMBER_LBM_POP_NODES = NX * NY;

// memory size
constexpr size_t MEM_SIZE_SCALAR = sizeof(dfloat) * NUMBER_LBM_POP_NODES;
constexpr size_t MEM_SIZE_POP = sizeof(dfloat) * NUMBER_LBM_POP_NODES * Q;
constexpr size_t MEM_SIZE_MOM = sizeof(dfloat) * NUMBER_LBM_NODES * NUMBER_MOMENTS;

constexpr dim3 threadBlock(BLOCK_NX, BLOCK_NY);
constexpr dim3 gridBlock(NUM_BLOCK_X, NUM_BLOCK_Y);

#endif //!__DEFINITIONS_H