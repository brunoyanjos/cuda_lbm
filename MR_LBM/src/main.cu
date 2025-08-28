
#include "main.cuh"
#include <iostream>
#include <chrono>
#include "saveData.cuh"

using namespace std;

int main()
{
	// return 0;
	folderSetup();

	// variable declaration
	unsigned int *node_type_fine, *node_type_coarse;
	dfloat *moments_fine, *moments_coarse;
	dfloat *pop_in_fine, *pop_out_fine, *pop_in_coarse, *pop_out_coarse;

	/* ------------------------- ALLOCATION FOR CPU ------------------------- */
	size_t step = 0;
	allocateHostMemory(&node_type_fine, &moments_fine, &pop_in_fine, &pop_out_fine,
					   &node_type_coarse, &moments_coarse, &pop_in_coarse, &pop_out_coarse);

	initializeDomain(node_type_fine, moments_fine, pop_in_fine, pop_out_fine,
					 node_type_coarse, moments_coarse, pop_in_coarse, pop_out_coarse);

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
		coarse_grid_solution(node_type_coarse, moments_coarse, pop_in_coarse, pop_out_coarse);

		// std::cout << std::endl;

		for (size_t fine_step = 0; fine_step < GRID_RATIO; ++fine_step)
		{
			fine_grid_solution(node_type_fine, moments_fine, pop_in_fine, pop_out_fine);

			// std::cout << std::endl;
		}

		// coarse_to_fine(moments_coarse, moments_fine, node_type_coarse, node_type_fine);
		// std::cout << "-------------------------------------------------------------------" << std::endl
		// 		  << std::endl;
		fine_to_coarse(moments_fine, moments_coarse, node_type_fine, node_type_coarse);
		// std::cout << "-------------------------------------------------------------------" << std::endl
		// 		  << std::endl;

		if (step % MACR_SAVE == 0)
		{
			printf("\n--------------------------- Last Time Step %06zu ---------------------------\n", step);

			saveMacr_coarse(moments_coarse, step, "002");
			saveMacr_fine(moments_fine, step, "003");
		}
	}

	/* ---------------------------------------------------------------------- */
	/* ------------------------------ END LOOP ------------------------------ */
	/* ---------------------------------------------------------------------- */

	// Calculate MLUPS
	dfloat MLUPS = recordElapsedTime(start_step, stop_step, step);
	printf("\n--------------------------- Last Time Step %06zu ---------------------------\n", step);
	printf("MLUPS: %f\n", MLUPS);

	coarse_velocity_profile(moments_coarse);
	fine_velocity_profile(moments_fine);

	/* ------------------------------ POST ------------------------------ */
	saveSimInfo(step, MLUPS);

	cudaFreeHost(node_type_fine);
	cudaFreeHost(moments_fine);
	cudaFreeHost(pop_in_fine);
	cudaFreeHost(pop_out_fine);

	cudaFreeHost(node_type_coarse);
	cudaFreeHost(moments_coarse);
	cudaFreeHost(pop_in_coarse);
	cudaFreeHost(pop_out_coarse);

	return 0;
}