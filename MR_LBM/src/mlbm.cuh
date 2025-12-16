#ifndef __MLBM_H
#define __MLBM_H

#include <string>
#include <math.h>
#include <cuda.h>
#include <curand.h>
#include <cuda_runtime.h>
#include <builtin_types.h>
#include "var.h"
#include "includeFiles/interface.h"
#include "boundaryCondition.cuh"

#include COLREC
#include CASE_BC

#include "globalStructs.h"

__global__ void streaming_and_moments(LBMState state, ghostInterfaceData ghostInterface, dfloat OMEGA);

__global__ void boundary_condition_and_interpolation(LBMState state, ghostInterfaceData ghostInterface, dfloat OMEGA);

__global__ void collision_and_interface_saving(LBMState state, ghostInterfaceData ghostInterface, dfloat OMEGA);

#endif