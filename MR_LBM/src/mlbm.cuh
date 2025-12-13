#ifndef __MLBM_H
#define __MLBM_H

#include <string>
#include <math.h>
#include <cuda.h>
#include <curand.h>
#include <cuda_runtime.h>
#include <builtin_types.h>
#include "var.h"
#include "globalStructs.h"
#include "init/state.cuh"

__global__ void streamingAndMom(LBMState state, dfloat OMEGA, ghostInterfaceData ghostInterface);

__global__ void updateBoundaries(LBMState state, dfloat OMEGA, dfloat D_in, dfloat D_out);

__global__ void boundaryAndCollision(LBMState state, dfloat OMEGA, ghostInterfaceData ghostInterface);

__global__ void mlbmKernel(LBMState state, dfloat OMEGA, ghostInterfaceData ghostInterface);

#endif