
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
		// for (size_t y = 0; y < NY; ++y)
		// {
		// 	unsigned int nodeType_one = fine_nodes[fine_idx(FINE_WIDTH - 1, FINE_WIDTH - 1 + GRID_RATIO * y)].node_type;
		// 	unsigned int nodeType_two = fine_nodes[fine_idx(NX_FINE_GRID - FINE_WIDTH, FINE_WIDTH - 1 + GRID_RATIO * y)].node_type;

		// 	fine_nodes[fine_idx(FINE_WIDTH - 1, FINE_WIDTH - 1 + GRID_RATIO * y)] = coarse_nodes[coarse_idx(N_OVERLAP_LAYER, y + N_OVERLAP_LAYER)];
		// 	fine_nodes[fine_idx(NX_FINE_GRID - FINE_WIDTH, FINE_WIDTH - 1 + GRID_RATIO * y)] = coarse_nodes[coarse_idx(NX_COARSE_GRID - N_OVERLAP_LAYER - 1, y + N_OVERLAP_LAYER)];

		// 	fine_nodes[fine_idx(FINE_WIDTH - 1, FINE_WIDTH - 1 + GRID_RATIO * y)].updated = true;
		// 	fine_nodes[fine_idx(NX_FINE_GRID - FINE_WIDTH, FINE_WIDTH - 1 + GRID_RATIO * y)].updated = true;

		// 	fine_nodes[fine_idx(FINE_WIDTH - 1, FINE_WIDTH - 1 + GRID_RATIO * y)].node_type = nodeType_one;
		// 	fine_nodes[fine_idx(NX_FINE_GRID - FINE_WIDTH, FINE_WIDTH - 1 + GRID_RATIO * y)].node_type = nodeType_two;
		// }

		// for (size_t x = 0; x < NX; ++x)
		// {
		// 	unsigned int nodeType_one = fine_nodes[fine_idx(FINE_WIDTH - 1 + GRID_RATIO * x, FINE_WIDTH - 1)].node_type;
		// 	unsigned int nodeType_two = fine_nodes[fine_idx(FINE_WIDTH - 1 + GRID_RATIO * x, NY_FINE_GRID - FINE_WIDTH)].node_type;

		// 	fine_nodes[fine_idx(FINE_WIDTH - 1 + GRID_RATIO * x, FINE_WIDTH - 1)] = coarse_nodes[coarse_idx(x + N_OVERLAP_LAYER, N_OVERLAP_LAYER)];
		// 	fine_nodes[fine_idx(FINE_WIDTH - 1 + GRID_RATIO * x, NY_FINE_GRID - FINE_WIDTH)] = coarse_nodes[coarse_idx(x + N_OVERLAP_LAYER, NX_COARSE_GRID - N_OVERLAP_LAYER - 1)];

		// 	fine_nodes[fine_idx(FINE_WIDTH - 1 + GRID_RATIO * x, FINE_WIDTH - 1)].updated = true;
		// 	fine_nodes[fine_idx(FINE_WIDTH - 1 + GRID_RATIO * x, NY_FINE_GRID - FINE_WIDTH)].updated = true;

		// 	fine_nodes[fine_idx(FINE_WIDTH - 1 + GRID_RATIO * x, FINE_WIDTH - 1)].node_type = nodeType_one;
		// 	fine_nodes[fine_idx(FINE_WIDTH - 1 + GRID_RATIO * x, NY_FINE_GRID - FINE_WIDTH)].node_type = nodeType_two;
		// }

		for (size_t fine_step = 0; fine_step < GRID_RATIO; ++fine_step)
		{
			fine_grid_solution(fine_nodes);
		}

		// int init_point = N_EXTRA_LAYER * GRID_RATIO;

		// for (size_t y = 0; y < NY_COARSE_GRID; ++y)
		// {
		// 	coarse_nodes[coarse_idx(0, y)] = fine_nodes[fine_idx(init_point, init_point + GRID_RATIO * y)];
		// 	coarse_nodes[coarse_idx(NX_COARSE_GRID - 1, y)] = fine_nodes[fine_idx(NX_FINE_GRID - init_point, init_point + GRID_RATIO * y)];
		// }

		// for (size_t x = 0; x < NX_COARSE_GRID; ++x)
		// {
		// 	coarse_nodes[coarse_idx(x, 0)] = fine_nodes[fine_idx(init_point + GRID_RATIO * x, init_point)];
		// 	coarse_nodes[coarse_idx(x, NY_COARSE_GRID - 1)] = fine_nodes[fine_idx(init_point + GRID_RATIO * x, NY_COARSE_GRID - init_point)];

		// 	coarse_nodes[coarse_idx(x, 0)].node_type = 100;
		// }

		coarse_grid_solution(coarse_nodes);

		std::ofstream file("GRID/001/macr_" + std::to_string(step) + ".dat");

		if (!file.is_open())
		{
			std::cerr << "Erro ao abrir arquivo para escrita!" << std::endl;
			return 0;
		}

		file << std::fixed << std::setprecision(12); // formatação com 6 casas decimais

		file << "x y rho ux uy\n";

		// for (size_t y = 0; y < NY_FINE_GRID; y++)
		// {
		// 	for (size_t x = 0; x < NX_FINE_GRID; x++)
		// 	{
		// 		file << x_coord << " " << y_coord << " " << rho << " " << ux << " " << uy << std::endl;
		// 	}
		// }

		file.close();
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