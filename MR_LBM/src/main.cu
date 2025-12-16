
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
	auto state = init_state();
	ghostInterfaceData ghostInterface;

	/* ----------------- GRID AND THREADS DEFINITION FOR LBM ---------------- */
	dim3 threadBlock(BLOCK_NX, BLOCK_NY);
	dim3 gridBlock(NUM_BLOCK_X, NUM_BLOCK_Y);

	/* ------------------------- ALLOCATION FOR CPU ------------------------- */
	int step = 0;
	int init_step = 0;

	/* -------------- ALLOCATION FOR GPU ------------- */
	interfaceMalloc(ghostInterface);

	// Setup Streams
	cudaStream_t streamsLBM[1];
	checkCudaErrors(cudaSetDevice(GPU_INDEX));
	checkCudaErrors(cudaStreamCreate(&streamsLBM[0]));
	checkCudaErrors(cudaDeviceSynchronize());

	bool success = initializeDomain(ghostInterface, state, gridBlock, threadBlock);

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

	dfloat VISC = U_MAX * (state.D_out - state.D_in) / RE;
	dfloat TAU = 0.5 + 3.0 * VISC; // relaxation time

	dfloat OMEGA = 1.0 / TAU; // (tau)^-1

	/* --------------------------------------------------------------------- */
	/* ---------------------------- BEGIN LOOP ------------------------------ */
	/* --------------------------------------------------------------------- */
	for (step = init_step; step <= N_STEPS; ++step)
	{
		streaming_and_moments<<<gridBlock, threadBlock>>>(state, ghostInterface, OMEGA);
		checkCudaErrors(cudaDeviceSynchronize());

		boundary_condition_and_interpolation<<<gridBlock, threadBlock>>>(state, ghostInterface, OMEGA);
		checkCudaErrors(cudaDeviceSynchronize());

		collision_and_interface_saving<<<gridBlock, threadBlock>>>(state, ghostInterface, OMEGA);
		checkCudaErrors(cudaDeviceSynchronize());

		swapGhostInterfaces(ghostInterface);
		checkCudaErrors(cudaDeviceSynchronize());

		if (step != 0 && step % MACR_SAVE == 0)
		{
			printf("\n----------------------------------- (%d/%d) %.2f%% -----------------------------------\n",
				   step, N_STEPS, static_cast<float>(step) / static_cast<float>(N_STEPS) * 100.0f);

			if (step != 0)
				time_elapsing_count(step_end, step_start, step);

			checkCudaErrors(cudaDeviceSynchronize());
			upload_state_to_host(state);

			create_vtk(state, step);
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
	// save info file
	saveSimInfo(step, MLUPS);

	/* ------------------------------ FREE ------------------------------ */

	interfaceFree(ghostInterface);
	free_state(state);

	return 0;
}