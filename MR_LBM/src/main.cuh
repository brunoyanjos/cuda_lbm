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

__host__ void allocateHostMemory(latticeNode **h_coarse_nodes, latticeNode **h_fine_nodes)
{
	checkCudaErrors(cudaMallocHost((void **)h_coarse_nodes, MEM_SIZE_NODES));
	checkCudaErrors(cudaMallocHost((void **)h_fine_nodes, MEM_SIZE_NODES));
}

__host__ void allocateDeviceMemory(latticeNode **d_coarse_nodes, latticeNode **d_fine_nodes, GhostInterfaceData *ghostInterface)
{
	cudaMalloc((void **)d_coarse_nodes, MEM_SIZE_NODES);
	cudaMalloc((void **)d_fine_nodes, MEM_SIZE_NODES);
	interfaceMalloc(*ghostInterface);
}

__host__ void initializeDomain(
	GhostInterfaceData &ghostInterface,
	latticeNode *&d_nodes, latticeNode *&h_nodes,
	int *step,
	dim3 gridBlock, dim3 threadBlock)
{
	// LBM Initialization
	gpuInitialization_nodes<<<gridBlock, threadBlock>>>(d_nodes);
	gpuInitialization_pop<<<gridBlock, threadBlock>>>(d_nodes, ghostInterface);

	gpuInitialization_nodeType_bulk<<<gridBlock, threadBlock>>>(d_nodes);
	gpuInitialization_nodeType<<<gridBlock, threadBlock>>>(d_nodes);

	// Interface population initialization
	interfaceCudaMemcpy(ghostInterface, ghostInterface.gGhost, ghostInterface.fGhost, cudaMemcpyDeviceToDevice, QF);

	// Synchronize after all initializations
	checkCudaErrors(cudaDeviceSynchronize());

	// Synchronize and transfer data back to host if needed
	checkCudaErrors(cudaDeviceSynchronize());
	checkCudaErrors(cudaMemcpy(h_nodes, d_nodes, sizeof(latticeNode) * NUMBER_LBM_NODES, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaDeviceSynchronize());
}

#endif // MAIN_CUH
