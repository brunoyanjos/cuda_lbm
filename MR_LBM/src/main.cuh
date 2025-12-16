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
#include "time_elapsing.cuh"
#include "checkpoint.cuh"
#include "treat_data.cuh"

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

__host__ void allocateHostMemory(
	dfloat **h_fMom, dfloat **rho, dfloat **ux, dfloat **uy, dfloat **probes, dfloat **ux_mean, dfloat **uy_mean)
{
	checkCudaErrors(cudaMallocHost((void **)h_fMom, MEM_SIZE_MOM));
	checkCudaErrors(cudaMallocHost((void **)rho, MEM_SIZE_SCALAR));
	checkCudaErrors(cudaMallocHost((void **)ux, MEM_SIZE_SCALAR));
	checkCudaErrors(cudaMallocHost((void **)uy, MEM_SIZE_SCALAR));
	checkCudaErrors(cudaMallocHost((void **)probes, MEM_SIZE_PROBES));
	checkCudaErrors(cudaMallocHost((void **)ux_mean, MEM_SIZE_UX_AVG));
	checkCudaErrors(cudaMallocHost((void **)uy_mean, MEM_SIZE_UY_AVG));
}

LBMState init_state()
{
	LBMState state;

	state.D_in = 0;
	state.D_out = NX + 2;

	// Defining Variable Sizes
	state.bytes_fields = NUMBER_LBM_NODES * sizeof(dfloat);
	state.bytes_types = NUMBER_LBM_NODES * sizeof(uint8_t);

	// Allocating HOST memory
	state.h_node_type = (uint8_t *)malloc(state.bytes_types);
	state.h_rho = (dfloat *)malloc(state.bytes_fields);

	state.h_ux = (dfloat *)malloc(state.bytes_fields);
	state.h_uy = (dfloat *)malloc(state.bytes_fields);

	state.h_mxx = (dfloat *)malloc(state.bytes_fields);
	state.h_mxy = (dfloat *)malloc(state.bytes_fields);
	state.h_myy = (dfloat *)malloc(state.bytes_fields);

	// Allocating DEVICE mesmory
	checkCudaErrors(cudaMalloc(&state.d_node_type, state.bytes_types));
	checkCudaErrors(cudaMalloc(&state.d_rho, state.bytes_fields));

	checkCudaErrors(cudaMalloc(&state.d_ux, state.bytes_fields));
	checkCudaErrors(cudaMalloc(&state.d_uy, state.bytes_fields));

	checkCudaErrors(cudaMalloc(&state.d_mxx, state.bytes_fields));
	checkCudaErrors(cudaMalloc(&state.d_mxy, state.bytes_fields));
	checkCudaErrors(cudaMalloc(&state.d_myy, state.bytes_fields));

	return state;
}

__host__ void allocateDeviceMemory(
	unsigned int **dNodeType, GhostInterfaceData *ghostInterface)
{
	cudaMalloc((void **)dNodeType, sizeof(int) * NUMBER_LBM_NODES);

	interfaceMalloc(*ghostInterface);
}

void upload_state_to_host(LBMState &state)
{
	checkCudaErrors(cudaMemcpy(state.h_rho, state.d_rho, state.bytes_fields, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(state.h_ux, state.d_ux, state.bytes_fields, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(state.h_uy, state.d_uy, state.bytes_fields, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(state.h_mxx, state.d_mxx, state.bytes_fields, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(state.h_mxy, state.d_mxy, state.bytes_fields, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(state.h_myy, state.d_myy, state.bytes_fields, cudaMemcpyDeviceToHost));
}

void free_state(LBMState &state)
{
	free(state.h_node_type);
	free(state.h_rho);
	free(state.h_ux);
	free(state.h_uy);
	free(state.h_mxx);
	free(state.h_mxy);
	free(state.h_myy);

	checkCudaErrors(cudaFree(state.d_node_type));
	checkCudaErrors(cudaFree(state.d_rho));
	checkCudaErrors(cudaFree(state.d_ux));
	checkCudaErrors(cudaFree(state.d_uy));
	checkCudaErrors(cudaFree(state.d_mxx));
	checkCudaErrors(cudaFree(state.d_mxy));
	checkCudaErrors(cudaFree(state.d_myy));
}

__host__ bool initializeDomain(
	GhostInterfaceData &ghostInterface,
	LBMState &state,
	dim3 gridBlock, dim3 threadBlock)
{
	// Node type initialization
	hostInitialization_nodeType(state.h_node_type);
	initialize_boundaries(state);

	checkCudaErrors(cudaMemcpy(state.d_node_type, state.h_node_type, state.bytes_types, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaDeviceSynchronize());

	// LBM Initialization
	gpuInitialization_mom<<<gridBlock, threadBlock>>>(state);
	gpuInitialization_pop<<<gridBlock, threadBlock>>>(state, ghostInterface);

	// Interface population initialization
	interfaceCudaMemcpy(ghostInterface, ghostInterface.gGhost, ghostInterface.fGhost, cudaMemcpyDeviceToDevice, QF);

	// Synchronize and transfer data back to host if needed
	checkCudaErrors(cudaDeviceSynchronize());
	upload_state_to_host(state);
	checkCudaErrors(cudaDeviceSynchronize());

	return true;
}

#endif // MAIN_CUH
