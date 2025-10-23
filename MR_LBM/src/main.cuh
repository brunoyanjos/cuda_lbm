// main.cuh
#ifndef MAIN_CUH
#define MAIN_CUH

#include <stdio.h>
#include <stdlib.h>

// CUDA INCLUDE
#include <cuda.h>
#include <cuda_runtime.h>
#include "device_launch_parameters.h"

// FILE INCLUDES
#include "var.h"
#include "errorDef.h"
#include "saveData.cuh"
#include "lbm_steps.cuh"
#include "lbm_solver.cuh"
#include "grid_data.cuh"

void initializeCudaEvents(cudaEvent_t &start, cudaEvent_t &stop, cudaEvent_t &start_step, cudaEvent_t &stop_step)
{
	checkCudaErrors(cudaSetDevice(GPU_INDEX));
	checkCudaErrors(cudaEventCreate(&start));
	checkCudaErrors(cudaEventCreate(&stop));
	checkCudaErrors(cudaEventCreate(&start_step));
	checkCudaErrors(cudaEventCreate(&stop_step));

	checkCudaErrors(cudaEventRecord(start, 0));
	checkCudaErrors(cudaEventRecord(start_step, 0));
}

dfloat recordElapsedTime(cudaEvent_t &start_step, cudaEvent_t &stop_step, int step)
{
	checkCudaErrors(cudaEventRecord(stop_step, 0));
	checkCudaErrors(cudaEventSynchronize(stop_step));

	float elapsedTime;
	checkCudaErrors(cudaEventElapsedTime(&elapsedTime, start_step, stop_step));
	elapsedTime *= 0.001;

	size_t nodesUpdatedSync = step * NUMBER_LBM_NODES;
	dfloat MLUPS = (nodesUpdatedSync / 1e6) / elapsedTime;
	return MLUPS;
}

__host__ void allocateHostMemory(unsigned int **node_type_fine, dfloat **moments_fine, dfloat **pop_in_fine, dfloat **pop_out_fine,
								 unsigned int **node_type_coarse, dfloat **moments_coarse, dfloat **pop_in_coarse, dfloat **pop_out_coarse)
{
	checkCudaErrors(cudaMallocHost((void **)node_type_fine, NUMBER_OF_FINE_NODES * sizeof(unsigned int)));
	checkCudaErrors(cudaMallocHost((void **)node_type_coarse, NUMBER_OF_COARSE_NODES * sizeof(unsigned int)));

	checkCudaErrors(cudaMallocHost((void **)moments_fine, NUMBER_OF_FINE_NODES * NUMBER_MOMENTS * sizeof(dfloat)));
	checkCudaErrors(cudaMallocHost((void **)moments_coarse, NUMBER_OF_COARSE_NODES * NUMBER_MOMENTS * sizeof(dfloat)));

	checkCudaErrors(cudaMallocHost((void **)pop_in_fine, NUMBER_OF_FINE_NODES * Q * sizeof(dfloat)));
	checkCudaErrors(cudaMallocHost((void **)pop_in_coarse, NUMBER_OF_COARSE_NODES * Q * sizeof(dfloat)));

	checkCudaErrors(cudaMallocHost((void **)pop_out_fine, NUMBER_OF_FINE_NODES * Q * sizeof(dfloat)));
	checkCudaErrors(cudaMallocHost((void **)pop_out_coarse, NUMBER_OF_COARSE_NODES * Q * sizeof(dfloat)));
}

__host__ inline void initialize_fine_grid(unsigned int *&node_type, dfloat *&moments, dfloat *&pop_in, dfloat *&pop_out)
{
	for (size_t y = 0; y < NY_FINE; ++y)
	{
		for (size_t x = 0; x < NX_FINE; ++x)
		{
			const dfloat rho = RHO_0;

			moments[idx_mom(x, y, M_RHO_INDEX, NX_FINE)] = rho - RHO_0;

			dfloat inv_rho = static_cast<dfloat>(1) / RHO_0;

			node_type[idx_grid(x, y, NX_FINE)] = fine_boundary_definition(x, y);

			moments[idx_mom(x, y, M_UX_INDEX, NX_FINE)] = node_type[idx_grid(x, y, NX_FINE)] == WEST ? U_MAX : static_cast<dfloat>(0);
			moments[idx_mom(x, y, M_UY_INDEX, NX_FINE)] = static_cast<dfloat>(0);

			dfloat pop[9];

			eval_pop_eq(pop, moments[idx_mom(x, y, M_RHO_INDEX, NX_FINE)],
						moments[idx_mom(x, y, M_UX_INDEX, NX_FINE)], moments[idx_mom(x, y, M_UY_INDEX, NX_FINE)]);

			dfloat mxx = (pop[1] + pop[3] + pop[5] + pop[6] + pop[7] + pop[8]) * inv_rho - cs2;
			dfloat mxy = ((pop[5] + pop[7]) - (pop[6] + pop[8])) * inv_rho;
			dfloat myy = (pop[2] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8]) * inv_rho - cs2;

			collision(&mxx, &mxy, &myy,
					  moments[idx_mom(x, y, M_UX_INDEX, NX_FINE)],
					  moments[idx_mom(x, y, M_UY_INDEX, NX_FINE)], OMEGA_COARSE);

			moments[idx_mom(x, y, M_MXX_INDEX, NX_FINE)] = (pop[1] + pop[3] + pop[5] + pop[6] + pop[7] + pop[8]) * inv_rho - cs2;
			moments[idx_mom(x, y, M_MXY_INDEX, NX_FINE)] = ((pop[5] + pop[7]) - (pop[6] + pop[8])) * inv_rho;
			moments[idx_mom(x, y, M_MYY_INDEX, NX_FINE)] = (pop[2] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8]) * inv_rho - cs2;
		}
	}
}

__host__ void initialize_coarse_grid(unsigned int *&node_type, dfloat *&moments, dfloat *&pop_in, dfloat *&pop_out)
{
	for (size_t y = 0; y < NY_COARSE; y++)
	{
		for (size_t x = 0; x < NX_COARSE; x++)
		{
			const dfloat rho = RHO_0;

			moments[idx_mom(x, y, M_RHO_INDEX, NX_COARSE)] = rho - RHO_0;

			dfloat inv_rho = static_cast<dfloat>(1) / RHO_0;

			node_type[idx_grid(x, y, NX_COARSE)] = coarse_boundary_definition(x, y);

			moments[idx_mom(x, y, M_UX_INDEX, NX_COARSE)] = node_type[idx_grid(x, y, NX_COARSE)] == WEST ? U_MAX : static_cast<dfloat>(0);
			moments[idx_mom(x, y, M_UY_INDEX, NX_COARSE)] = static_cast<dfloat>(0);

			dfloat pop[9];

			eval_pop_eq(pop, moments[idx_mom(x, y, M_RHO_INDEX, NX_COARSE)],
						moments[idx_mom(x, y, M_UX_INDEX, NX_COARSE)], moments[idx_mom(x, y, M_UY_INDEX, NX_COARSE)]);

			dfloat mxx = (pop[1] + pop[3] + pop[5] + pop[6] + pop[7] + pop[8]) * inv_rho - cs2;
			dfloat mxy = ((pop[5] + pop[7]) - (pop[6] + pop[8])) * inv_rho;
			dfloat myy = (pop[2] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8]) * inv_rho - cs2;

			collision(&mxx, &mxy, &myy,
					  moments[idx_mom(x, y, M_UX_INDEX, NX_COARSE)],
					  moments[idx_mom(x, y, M_UY_INDEX, NX_COARSE)], OMEGA_COARSE);

			moments[idx_mom(x, y, M_MXX_INDEX, NX_COARSE)] = (pop[1] + pop[3] + pop[5] + pop[6] + pop[7] + pop[8]) * inv_rho - cs2;
			moments[idx_mom(x, y, M_MXY_INDEX, NX_COARSE)] = ((pop[5] + pop[7]) - (pop[6] + pop[8])) * inv_rho;
			moments[idx_mom(x, y, M_MYY_INDEX, NX_COARSE)] = (pop[2] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8]) * inv_rho - cs2;
		}
	}
}

__host__ void initializeDomain(unsigned int *&node_type_fine, dfloat *&moments_fine, dfloat *&pop_in_fine, dfloat *&pop_out_fine,
							   unsigned int *&node_type_coarse, dfloat *&moments_coarse, dfloat *&pop_in_coarse, dfloat *&pop_out_coarse)
{
	initialize_fine_grid(node_type_fine, moments_fine, pop_in_fine, pop_out_fine);
	initialize_coarse_grid(node_type_coarse, moments_coarse, pop_in_coarse, pop_out_coarse);
}

#endif // MAIN_CUH
