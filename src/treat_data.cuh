#ifndef TREAT_DATA
#define TREAT_DATA

#include <cuda.h>
#include <curand.h>
#include <cuda_runtime.h>
#include <builtin_types.h>

#include "var.h"
#include "globalFunctions.h"

__global__ void velocity_average(dfloat *fMom, dfloat *ux_mean, dfloat *uy_mean, unsigned int step);

#endif