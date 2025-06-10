
#include "main.cuh"
#include <iostream>
#include <chrono>
#include "saveData.cuh"

using namespace std;

std::map<unsigned int, unsigned int> cylinder_index;

int main()
{
	printf("BLOCK_NX: %d, BLOCK_NY: %d\n", BLOCK_NX, BLOCK_NY);

	folderSetup();

	// set cuda device
	checkCudaErrors(cudaSetDevice(GPU_INDEX));

	// variable declaration
	dfloat* d_fMom;
	ghostInterfaceData ghostInterface;

	unsigned int* dNodeType;
	unsigned int* hNodeType;

	dfloat* h_fMom;

	dfloat* rho;

	dfloat* ux;
	dfloat* uy;

	dfloat* mxx;
	dfloat* mxy;
	dfloat* myy;

	/* ----------------- GRID AND THREADS DEFINITION FOR LBM ---------------- */
	dim3 threadBlock(BLOCK_NX, BLOCK_NY);
	dim3 gridBlock(NUM_BLOCK_X, NUM_BLOCK_Y);

	/* ------------------------- ALLOCATION FOR CPU ------------------------- */
	int step = 0;

	allocateHostMemory(&h_fMom, &rho, &ux, &uy, &mxx, &mxy, &myy);

	/* -------------- ALLOCATION FOR GPU ------------- */
	allocateDeviceMemory(&d_fMom, &dNodeType, &ghostInterface);

	// Setup Streams
	cudaStream_t streamsLBM[1];
	checkCudaErrors(cudaSetDevice(GPU_INDEX));
	checkCudaErrors(cudaStreamCreate(&streamsLBM[0]));
	checkCudaErrors(cudaDeviceSynchronize());

	initializeDomain(ghostInterface, d_fMom, h_fMom, hNodeType, dNodeType, &step, gridBlock, threadBlock); 

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
		gpuMomCollisionStream << <gridBlock, threadBlock >> > (d_fMom, dNodeType, ghostInterface, step);

		// swap interface pointers
		swapGhostInterfaces(ghostInterface);

		if (MACR_SAVE != 0 && step % MACR_SAVE == 0) {
			printf("\n----------------------------------- %d -----------------------------------\n", step);

			checkCudaErrors(cudaDeviceSynchronize());
			checkCudaErrors(cudaMemcpy(h_fMom, d_fMom, sizeof(dfloat) * NUMBER_LBM_NODES * NUMBER_MOMENTS, cudaMemcpyDeviceToHost));
			saveMacr(h_fMom, rho, ux, uy, step);
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
	checkCudaErrors(cudaMemcpy(h_fMom, d_fMom, sizeof(dfloat) * NUMBER_LBM_NODES * NUMBER_MOMENTS, cudaMemcpyDeviceToHost));
	// save info file
	saveSimInfo(step, MLUPS);

	/* ------------------------------ FREE ------------------------------ */
	cudaFree(d_fMom);
	cudaFree(dNodeType);
	cudaFree(hNodeType);
	cudaFree(hNodeType);
	cudaFree(h_fMom);
	cudaFree(rho);
	cudaFree(ux);
	cudaFree(uy);
	cudaFree(mxx);
	cudaFree(mxy);
	cudaFree(myy);
	interfaceFree(ghostInterface);
	return 0;
}