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

__global__ void mlbmKernel(LBMState state, ghostInterfaceData ghostInterface);

#endif