
#include "main.cuh"
#include <iostream>
#include <chrono>
#include "saveData.cuh"

using namespace std;

int main()
{
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

	/* --------------------------------------------------------------------- */
	/* ---------------------------- BEGIN LOOP ------------------------------ */
	/* --------------------------------------------------------------------- */

	for (step = INI_STEP; step < N_STEPS; step++)
	{
		coarse_grid_solution(coarse_nodes);

		coarse_to_fine(coarse_nodes, fine_nodes);

		for (size_t fine_step = 0; fine_step < GRID_RATIO; ++fine_step)
		{
			fine_grid_solution(fine_nodes);
		}

		fine_to_coarse(fine_nodes, coarse_nodes);

		string sim_id = ID_SIM;

		std::ofstream coarse_file("GRID/" + sim_id + "/coarse_macr_" + std::to_string(step) + ".dat");
		std::ofstream fine_file("GRID/" + sim_id + "/fine_macr_" + std::to_string(step) + ".dat");

		coarse_file << std::fixed << std::setprecision(12); // formatação com 6 casas decimais
		fine_file << std::fixed << std::setprecision(12);

		for (size_t y = 0; y < NY_COARSE; y++)
		{
			for (size_t x = 0; x < NX_COARSE; x++)
			{
				const dfloat ux = coarse_nodes[coarse_idx(x, y)].ux / F_M_I_SCALE;
				const dfloat uy = coarse_nodes[coarse_idx(x, y)].uy / F_M_I_SCALE;

				const dfloat u2 = ux * ux + uy * uy;

				const dfloat u = std::sqrt(u2);

				coarse_file << u << " ";
			}
			coarse_file << std::endl;
		}

		for (size_t y = 0; y < NY_FINE; ++y)
		{
			for (size_t x = 0; x < NX_FINE; ++x)
			{
				const dfloat ux = fine_nodes[fine_idx(x, y)].ux / F_M_I_SCALE;
				const dfloat uy = fine_nodes[fine_idx(x, y)].uy / F_M_I_SCALE;

				const dfloat u2 = ux * ux + uy * uy;

				const dfloat u = std::sqrt(u2);

				fine_file << u << " ";
			}
			fine_file << std::endl;
		}

		coarse_file.close();
		fine_file.close();
	}

	/* --------------------------------------------------------------------- */
	/* ------------------------------ END LOOP ------------------------------ */
	/* --------------------------------------------------------------------- */

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