
#include "main.cuh"
#include <iostream>
#include <chrono>
#include "saveData.cuh"

using namespace std;

int main()
{
	printf("BLOCK_NX: %d, BLOCK_NY: %d\n", BLOCK_NX, BLOCK_NY);

	folderSetup();

	// set cuda device
	checkCudaErrors(cudaSetDevice(GPU_INDEX));

	// variable declaration
	latticeNode *d_nodes;
	latticeNode *h_nodes;
	ghostInterfaceData ghostInterface;

	/* ----------------- GRID AND THREADS DEFINITION FOR LBM ---------------- */
	dim3 threadBlock(BLOCK_NX, BLOCK_NY);
	dim3 gridBlock(NUM_BLOCK_X, NUM_BLOCK_Y);

	/* ------------------------- ALLOCATION FOR CPU ------------------------- */
	int step = 0;

	allocateHostMemory(&h_nodes);

	/* -------------- ALLOCATION FOR GPU ------------- */
	allocateDeviceMemory(&d_nodes, &ghostInterface);

	// Setup Streams
	cudaStream_t streamsLBM[1];
	checkCudaErrors(cudaSetDevice(GPU_INDEX));
	checkCudaErrors(cudaStreamCreate(&streamsLBM[0]));
	checkCudaErrors(cudaDeviceSynchronize());

	initializeDomain(ghostInterface, d_nodes, h_nodes, &step, gridBlock, threadBlock);

	/* ------------------------------ TIMER EVENTS  ------------------------------ */
	checkCudaErrors(cudaSetDevice(GPU_INDEX));
	cudaEvent_t start, stop, start_step, stop_step;
	initializeCudaEvents(start, stop, start_step, stop_step);
	/* ------------------------------ LBM LOOP ------------------------------ */
	saveSimInfo(step, 0.0);

	/* --------------------------------------------------------------------- */
	/* ---------------------------- BEGIN LOOP ------------------------------ */
	/* --------------------------------------------------------------------- */
	for (step = INI_STEP; step < N_STEPS; step++)
	{
		gpuMomCollisionStream<<<gridBlock, threadBlock>>>(d_nodes, ghostInterface, step);

		// swap interface pointers
		swapGhostInterfaces(ghostInterface);

		if (MACR_SAVE != 0 && step % MACR_SAVE == 0)
		{
			printf("\n----------------------------------- %d -----------------------------------\n", step);

			checkCudaErrors(cudaDeviceSynchronize());
			checkCudaErrors(cudaMemcpy(h_nodes, d_nodes, sizeof(latticeNode) * NUMBER_LBM_NODES, cudaMemcpyDeviceToHost));

			kinetic_energy(h_nodes, step);
			//saveMacr(h_nodes, step);
		}
	}

	/* --------------------------------------------------------------------- */
	/* ------------------------------ END LOOP ------------------------------ */
	/* --------------------------------------------------------------------- */

	checkCudaErrors(cudaDeviceSynchronize());

	// Calculate MLUPS

	dfloat MLUPS = recordElapsedTime(start_step, stop_step, step);
	printf("\n--------------------------- Last Time Step %06d ---------------------------\n", step);
	printf("MLUPS: %f\n", MLUPS);

	/* ------------------------------ POST ------------------------------ */
	checkCudaErrors(cudaMemcpy(h_nodes, d_nodes, sizeof(latticeNode) * NUMBER_LBM_NODES, cudaMemcpyDeviceToHost));
	// save info file
	saveSimInfo(step, MLUPS);
	velocity_profiles(h_nodes, step);

	/* ------------------------------ FREE ------------------------------ */
	cudaFree(d_nodes);
	cudaFree(h_nodes);
	interfaceFree(ghostInterface);
	return 0;
}