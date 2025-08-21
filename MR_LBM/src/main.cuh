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
#include "globalStructs.h"
#include "errorDef.h"
#include "lbmInitialization.cuh"
#include "mlbm.cuh"
#include "saveData.cuh"
#include "lbm_solver.cuh"

#include "grid_data.cuh"

/*
 *   @brief Swaps the pointers of two dfloat variables.
 *   @param pt1: reference to the first dfloat pointer to be swapped
 *   @param pt2: reference to the second dfloat pointer to be swapped
 */
__host__ __device__ void interfaceSwap(dfloat *&pt1, dfloat *&pt2)
{
	dfloat *temp = pt1;
	pt1 = pt2;
	pt2 = temp;
}

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

/*
 *   @brief Frees the memory allocated for the ghost interface data.
 *   @param ghostInterface: reference to the ghost interface data structure
 */
__host__ void interfaceFree(ghostInterfaceData &ghostInterface)
{
	cudaFree(ghostInterface.fGhost.X_0);
	cudaFree(ghostInterface.fGhost.X_1);
	cudaFree(ghostInterface.fGhost.Y_0);
	cudaFree(ghostInterface.fGhost.Y_1);

	cudaFree(ghostInterface.gGhost.X_0);
	cudaFree(ghostInterface.gGhost.X_1);
	cudaFree(ghostInterface.gGhost.Y_0);
	cudaFree(ghostInterface.gGhost.Y_1);
}

/*
 *   @brief Performs a CUDA memory copy for ghost interface data between source and destination.
 *   @param ghostInterface: reference to the ghost interface data structure
 *   @param dst: destination ghost data structure
 *   @param src: source ghost data structure
 *   @param kind: type of memory copy (e.g., cudaMemcpyHostToDevice)
 *   @param Q: number of quantities in the ghost data that are transfered
 */
__host__ void interfaceCudaMemcpy(GhostInterfaceData &ghostInterface, ghostData &dst, const ghostData &src, cudaMemcpyKind kind, int Q)
{
	struct MemcpyPair
	{
		dfloat *dst;
		const dfloat *src;
		size_t size;
	};

	MemcpyPair memcpyPairs[] = {
		{dst.X_0, src.X_0, sizeof(dfloat) * NUMBER_GHOST_FACE_X * Q},
		{dst.X_1, src.X_1, sizeof(dfloat) * NUMBER_GHOST_FACE_X * Q},
		{dst.Y_0, src.Y_0, sizeof(dfloat) * NUMBER_GHOST_FACE_Y * Q},
		{dst.Y_1, src.Y_1, sizeof(dfloat) * NUMBER_GHOST_FACE_Y * Q},
	};

	checkCudaErrors(cudaDeviceSynchronize());
	for (const auto &pair : memcpyPairs)
	{
		checkCudaErrors(cudaMemcpy(pair.dst, pair.src, pair.size, kind));
	}
}

/*
 *   @brief Swaps the ghost interfaces.
 *   @param ghostInterface: reference to the ghost interface data structure
 */
__host__ void swapGhostInterfaces(GhostInterfaceData &ghostInterface)
{
	// Synchronize device before performing swaps
	checkCudaErrors(cudaDeviceSynchronize());

	// Swap interface pointers for fGhost and gGhost
	interfaceSwap(ghostInterface.fGhost.X_0, ghostInterface.gGhost.X_0);
	interfaceSwap(ghostInterface.fGhost.X_1, ghostInterface.gGhost.X_1);
	interfaceSwap(ghostInterface.fGhost.Y_0, ghostInterface.gGhost.Y_0);
	interfaceSwap(ghostInterface.fGhost.Y_1, ghostInterface.gGhost.Y_1);
}

/*
 *   @brief Allocates memory for the ghost interface data.
 *   @param ghostInterface: reference to the ghost interface data structure
 */
__host__ void interfaceMalloc(ghostInterfaceData &ghostInterface)
{
	cudaMalloc((void **)&(ghostInterface.fGhost.X_0), sizeof(dfloat) * NUMBER_GHOST_FACE_X * QF);
	cudaMalloc((void **)&(ghostInterface.fGhost.X_1), sizeof(dfloat) * NUMBER_GHOST_FACE_X * QF);
	cudaMalloc((void **)&(ghostInterface.fGhost.Y_0), sizeof(dfloat) * NUMBER_GHOST_FACE_Y * QF);
	cudaMalloc((void **)&(ghostInterface.fGhost.Y_1), sizeof(dfloat) * NUMBER_GHOST_FACE_Y * QF);

	cudaMalloc((void **)&(ghostInterface.gGhost.X_0), sizeof(dfloat) * NUMBER_GHOST_FACE_X * QF);
	cudaMalloc((void **)&(ghostInterface.gGhost.X_1), sizeof(dfloat) * NUMBER_GHOST_FACE_X * QF);
	cudaMalloc((void **)&(ghostInterface.gGhost.Y_0), sizeof(dfloat) * NUMBER_GHOST_FACE_Y * QF);
	cudaMalloc((void **)&(ghostInterface.gGhost.Y_1), sizeof(dfloat) * NUMBER_GHOST_FACE_Y * QF);
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

__host__ void allocateDeviceMemory(latticeNode **d_coarse_nodes, latticeNode **d_fine_nodes, GhostInterfaceData *ghostInterface)
{
	cudaMalloc((void **)d_coarse_nodes, MEM_SIZE_NODES);
	cudaMalloc((void **)d_fine_nodes, MEM_SIZE_NODES);
	interfaceMalloc(*ghostInterface);
}

__host__ inline void initialize_fine_grid(unsigned int *&node_type, dfloat *&moments, dfloat *&pop_in, dfloat *&pop_out)
{
	for (size_t y = 0; y < NY_FINE; ++y)
	{
		for (size_t x = 0; x < NX_FINE; ++x)
		{
			moments[fine_moment_idx(x, y, M_RHO_INDEX)] = RHO_0;

			dfloat inv_rho = static_cast<dfloat>(1) / RHO_0;

			moments[fine_moment_idx(x, y, M_UX_INDEX)] = 0.0;
			moments[fine_moment_idx(x, y, M_UY_INDEX)] = 0.0;

			node_type[fine_idx(x, y)] = BULK;

			if (y == 0)
			{
				node_type[fine_idx(x, y)] = SOUTH;
			}

			dfloat pop[9];

			init_pop_eq(pop, moments[fine_moment_idx(x, y, M_RHO_INDEX)],
						moments[fine_moment_idx(x, y, M_UX_INDEX)], moments[fine_moment_idx(x, y, M_UY_INDEX)]);

			moments[fine_moment_idx(x, y, M_MXX_INDEX)] = (pop[1] + pop[3] + pop[5] + pop[6] + pop[7] + pop[8]) * inv_rho - cs2;
			moments[fine_moment_idx(x, y, M_MXY_INDEX)] = ((pop[5] + pop[7]) - (pop[6] + pop[8])) * inv_rho;
			moments[fine_moment_idx(x, y, M_MYY_INDEX)] = (pop[2] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8]) * inv_rho - cs2;
		}
	}
}

__host__ void initialize_coarse_grid(unsigned int *&node_type, dfloat *&moments, dfloat *&pop_in, dfloat *&pop_out)
{
	for (size_t y = 0; y < NY_COARSE + N_OVERLAP_LAYER; y++)
	{
		for (size_t x = 0; x < NX_COARSE; x++)
		{
			moments[coarse_moment_idx(x, y, M_RHO_INDEX)] = RHO_0;

			dfloat inv_rho = static_cast<dfloat>(1) / RHO_0;

			moments[coarse_moment_idx(x, y, M_UX_INDEX)] = 0.0;
			moments[coarse_moment_idx(x, y, M_UY_INDEX)] = 0.0;

			node_type[coarse_idx(x, y)] = BULK;

			if (y == NY_COARSE + N_OVERLAP_LAYER - 1)
			{
				node_type[coarse_idx(x, y)] = NORTH;
				moments[coarse_moment_idx(x, y, M_UX_INDEX)] = U_MAX;
			}

			dfloat pop[9];

			init_pop_eq(pop, moments[coarse_moment_idx(x, y, M_RHO_INDEX)],
						moments[coarse_moment_idx(x, y, M_UX_INDEX)], moments[coarse_moment_idx(x, y, M_UY_INDEX)]);

			moments[coarse_moment_idx(x, y, M_MXX_INDEX)] = (pop[1] + pop[3] + pop[5] + pop[6] + pop[7] + pop[8]) * inv_rho - cs2;
			moments[coarse_moment_idx(x, y, M_MXY_INDEX)] = ((pop[5] + pop[7]) - (pop[6] + pop[8])) * inv_rho;
			moments[coarse_moment_idx(x, y, M_MYY_INDEX)] = (pop[2] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8]) * inv_rho - cs2;
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
