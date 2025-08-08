
#include "main.cuh"
#include <iostream>
#include <chrono>
#include "saveData.cuh"

using namespace std;

int main()
{
	std::cout << "TAU_FINE: " << TAU_FINE << std::endl;
	std::cout << "TAU_COARSE: " << TAU_COARSE << std::endl;

	folderSetup();

	// variable declaration
	latticeNode *fine_nodes;
	latticeNode *coarse_nodes;

	/* ------------------------- ALLOCATION FOR CPU ------------------------- */
	size_t step = 0;
	allocateHostMemory(&coarse_nodes, &fine_nodes);
	initializeDomain(fine_nodes, coarse_nodes);

	/* ------------------------------ TIMER EVENTS  ------------------------------ */
	checkCudaErrors(cudaSetDevice(GPU_INDEX));
	cudaEvent_t start, stop, start_step, stop_step;
	initializeCudaEvents(start, stop, start_step, stop_step);
	/* ------------------------------ LBM LOOP ------------------------------ */
	saveSimInfo(step, 0.0);

	/* ---------------------------------------------------------------------- */
	/* ---------------------------- BEGIN LOOP ------------------------------ */
	/* ---------------------------------------------------------------------- */
	for (step = INI_STEP; step < N_STEPS; step++)
	{
		coarse_to_fine(coarse_nodes, fine_nodes);
		fine_to_coarse(fine_nodes, coarse_nodes);

		coarse_grid_solution(coarse_nodes);

		for (size_t fine_step = 0; fine_step < GRID_RATIO; ++fine_step)
		{
			fine_grid_solution(fine_nodes, fine_step == 0);
		}

		if (step % MACR_SAVE == 0)
		{
			saveMacr_coarse(coarse_nodes, step, "002");
			saveMacr_fine(fine_nodes, step, "003");
		}
	}

	/* ---------------------------------------------------------------------- */
	/* ------------------------------ END LOOP ------------------------------ */
	/* ---------------------------------------------------------------------- */

	// Calculate MLUPS
	dfloat MLUPS = recordElapsedTime(start_step, stop_step, step);
	printf("\n--------------------------- Last Time Step %06zu ---------------------------\n", step);
	printf("MLUPS: %f\n", MLUPS);

	/* ------------------------------ POST ------------------------------ */
	saveSimInfo(step, MLUPS);

	cudaFreeHost(fine_nodes);
	cudaFreeHost(coarse_nodes);

	return 0;
}