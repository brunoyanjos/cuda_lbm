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

/*
 *   @brief Updates macroscopics and then performs collision and streaming
 *   @param fMom: macroscopics moments
 *   @param ghostInterface interface block transfer information
 *   @param d_mean_rho: mean density, used for density correction
 *   @param d_BC_Fx: boundary condition force x
 *   @param d_BC_Fy: boundary condition force x
 *   @param d_BC_Fz: boundary condition force x
 *   @param step: current time step
 *   @param save: if is necessary save some data
 */
__global__ void gpuMomCollisionStream(
	dfloat* fMom, unsigned int* dNodeType, ghostInterfaceData ghostInterface, unsigned int step
);


#endif __MLBM_H