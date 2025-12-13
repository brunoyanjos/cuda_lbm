#include "saveData.cuh"
#include "interface/ghost_interface.cuh"
#include "mlbm.cuh"
#include "init/state.cuh"
#include "globalFunctions.h"
#include "init/domain.cuh"
#include "time_events.cuh"
#include "post_processing.cuh"
#include "boundaries/node_type.h"

int main()
{
	folderSetup();

	// set cuda device
	checkCudaErrors(cudaSetDevice(GPU_INDEX));

	// variable declaration

	ghostInterfaceData ghostInterface;

	/* ------------------------- ALLOCATION FOR CPU ------------------------- */
	dfloat D_out = D * 2, D_in = 0;

	int step = 0;
	auto state = init_state();

	interfaceMalloc(ghostInterface);

	init_domain(state, ghostInterface, D_out, D_in);

	const dfloat VISC = U_MAX * NX / RE;
	const dfloat TAU = 0.5 + 3.0 * VISC;
	const dfloat OMEGA = 1.0 / TAU;

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
		mlbmKernel<<<gridBlock, threadBlock>>>(state, OMEGA, ghostInterface);
		checkCudaErrors(cudaDeviceSynchronize());

		// streamingAndMom<<<gridBlock, threadBlock>>>(state, OMEGA, ghostInterface);
		// checkCudaErrors(cudaDeviceSynchronize());

		// updateBoundaries<<<gridBlock, threadBlock>>>(state, OMEGA, D_in, D_out);
		// checkCudaErrors(cudaDeviceSynchronize());

		// boundaryAndCollision<<<gridBlock, threadBlock>>>(state, OMEGA, ghostInterface);
		// checkCudaErrors(cudaDeviceSynchronize());

		swapGhostInterfaces(ghostInterface);
		checkCudaErrors(cudaDeviceSynchronize());

		if (step % MACR_SAVE == 0)
		{
			printf("\n----------------------------------- (%d/%d) %.2f%% -----------------------------------\n",
				   step, N_STEPS, static_cast<float>(step) / static_cast<float>(N_STEPS) * 100.0f);

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

	upload_state_to_host(state);
	write_average_velocity_profile(state.h_ux);
	create_vtk(state, step);

	/* ------------------------------ POST ------------------------------ */
	// save info file
	saveSimInfo(step, MLUPS);

	/* ------------------------------ FREE ------------------------------ */
	free_state(state);
	interfaceFree(ghostInterface);
	return 0;
}