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

__global__ void streamingAndMom(
	dfloat *fMom, dfloat OMEGA, size_t cylinder_counter, unsigned int *dNodeType,
	ghostInterfaceData ghostInterface, unsigned int step);

__global__ void updateInnerBoundaries(dfloat *fMom, dfloat OMEGA, unsigned int step);

__global__ void boundaryAndCollision(
	dfloat *fMom, size_t cylinder_count, dfloat OMEGA, unsigned int *dNodeType,
	ghostInterfaceData ghostInterface, unsigned int step);

#endif