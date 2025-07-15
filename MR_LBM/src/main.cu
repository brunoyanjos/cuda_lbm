
#include "main.cuh"
#include "saveData.cuh"

using namespace std;

int main()
{
	printf("BLOCK_NX: %d, BLOCK_NY: %d\n", BLOCK_NX, BLOCK_NY);

	folderSetup();

	// set cuda device
	checkCudaErrors(cudaSetDevice(GPU_INDEX));

	// variable declaration
	dfloat *d_fMom;
	ghostInterfaceData ghostInterface;

	unsigned int *dNodeType;
	unsigned int *hNodeType;

	dfloat *h_fMom;
	dfloat *rho;

	dfloat *ux;
	dfloat *uy;

	dfloat *probes;

	/* ----------------- GRID AND THREADS DEFINITION FOR LBM ---------------- */
	dim3 threadBlock(BLOCK_NX, BLOCK_NY);
	dim3 gridBlock(NUM_BLOCK_X, NUM_BLOCK_Y);

	/* ------------------------- ALLOCATION FOR CPU ------------------------- */
	int step = 0;
	int init_step = 0;

	allocateHostMemory(&h_fMom, &rho, &ux, &uy, &probes);

	/* -------------- ALLOCATION FOR GPU ------------- */
	allocateDeviceMemory(&d_fMom, &dNodeType, &ghostInterface);

	// Setup Streams
	cudaStream_t streamsLBM[1];
	checkCudaErrors(cudaSetDevice(GPU_INDEX));
	checkCudaErrors(cudaStreamCreate(&streamsLBM[0]));
	checkCudaErrors(cudaDeviceSynchronize());

	bool success = initializeDomain(ghostInterface, d_fMom, h_fMom, hNodeType, dNodeType, &init_step, gridBlock, threadBlock);

	if (!success)
	{
		return 0;
	}

	/* ------------------------------ TIMER EVENTS  ------------------------------ */
	checkCudaErrors(cudaSetDevice(GPU_INDEX));
	cudaEvent_t start, stop, start_step, stop_step;
	initializeCudaEvents(start, stop, start_step, stop_step);
	/* ------------------------------ LBM LOOP ------------------------------ */
	saveSimInfo(step, 0.0);

	timestep sim_start_time = std::chrono::high_resolution_clock::now();
	timestep step_start = std::chrono::high_resolution_clock::now();
	timestep step_end;

	/* --------------------------------------------------------------------- */
	/* ---------------------------- BEGIN LOOP ------------------------------ */
	/* --------------------------------------------------------------------- */
	for (step = init_step; step <= N_STEPS; ++step)
	{
		gpuMomCollisionStream<<<gridBlock, threadBlock>>>(d_fMom, dNodeType, ghostInterface, step);

		// swap interface pointers
		swapGhostInterfaces(ghostInterface);

		if (step % CHECKPOINT_STEP == 0 && step != init_step)
		{
			checkCudaErrors(cudaDeviceSynchronize());
			checkCudaErrors(cudaMemcpy(h_fMom, d_fMom, sizeof(dfloat) * NUMBER_LBM_NODES * NUMBER_MOMENTS, cudaMemcpyDeviceToHost));

			save_checkpoint(step, h_fMom);
		}

		if (MACR_SAVE != 0 && step % MACR_SAVE == 0)
		{
			printf("\n----------------------------------- (%d/%d) %.2f%% -----------------------------------\n", step, N_STEPS, static_cast<float>(step) / static_cast<float>(N_STEPS) * 100.0f);
			if (step != 0) time_elapsing_count(step_end, step_start, step);
			
			checkCudaErrors(cudaDeviceSynchronize());
			checkCudaErrors(cudaMemcpy(h_fMom, d_fMom, sizeof(dfloat) * NUMBER_LBM_NODES * NUMBER_MOMENTS, cudaMemcpyDeviceToHost));

			kinetic_energy(h_fMom, step);
			saving_probes(h_fMom, probes, step);
			// saveMacr(h_fMom, rho, ux, uy, step);
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

	velocity_profiles(h_fMom, step);

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
	interfaceFree(ghostInterface);
	return 0;
}