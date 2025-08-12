
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

	dfloat sys_mass_in = 0.0;
	dfloat sys_mass_out = 0.0;
	dfloat coarse_mass_in = 0.0;
	dfloat fine_mass_in = 0.0;
	dfloat coarse_mass_out = 0.0;
	dfloat fine_mass_out = 0.0;

	bool calculate_mass = false;

	/* ---------------------------------------------------------------------- */
	/* ---------------------------- BEGIN LOOP ------------------------------ */
	/* ---------------------------------------------------------------------- */
	for (step = INI_STEP; step < N_STEPS; step++)
	{
		coarse_to_fine(coarse_nodes, fine_nodes);
		fine_to_coarse(fine_nodes, coarse_nodes);

		coarse_grid_solution(coarse_nodes);

		if (calculate_mass)
		{
			// Calculo das densidade de entrada da malha grossa
			sys_mass_in = 0.0;
			sys_mass_out = 0.0;

			coarse_mass_in = 0.0;
			coarse_mass_out = 0.0;

			for (size_t y = 0; y < NY_COARSE; ++y)
			{
				for (size_t x = 0; x < NX_COARSE; ++x)
				{
					size_t idx = coarse_idx(x, y);
					latticeNode node = coarse_nodes[idx];

					const arrayType<9> pop = regularization_mass(node);

					dfloat node_mass_in = 0.0;
					dfloat node_mass_out = 0.0;

					for (size_t i = 0; i < 9; ++i)
					{
						if (node.incomings[i])
						{
							node_mass_in += node.pop_in[i];
						}

						if (node.outgoings[i])
						{
							node_mass_out += pop.f[i];
						}
					}

					coarse_mass_in += node_mass_in;
					coarse_mass_out += node_mass_out;
				}
			}

			sys_mass_in += coarse_mass_in;
			sys_mass_out += coarse_mass_out;

			std::cout << "coarse_mass_in total: " << coarse_mass_in << std::endl;
			std::cout << "coarse_mass_in avg: " << coarse_mass_in / (NY_COARSE * NX_COARSE) << std::endl;

			std::cout << "coarse_mass_out total: " << coarse_mass_out << std::endl;
			std::cout << "coarse_mass_out avg: " << coarse_mass_out / (NY_COARSE * NX_COARSE) << std::endl;

			std::cout << "coarse_mass_net total: " << coarse_mass_out - coarse_mass_in << std::endl
					  << std::endl;

			// ---------------------------------------------
		}

		for (size_t fine_step = 0; fine_step < GRID_RATIO; ++fine_step)
		{

			fine_grid_solution(fine_nodes, fine_step == 0);

			if (calculate_mass)
			{
				// Calculo das densidade de entrada da malha fina
				fine_mass_in = 0.0;
				fine_mass_out = 0.0;

				for (size_t y = 0; y < NY_FINE; ++y)
				{
					for (size_t x = 0; x < NX_FINE; ++x)
					{
						size_t idx = fine_idx(x, y);
						latticeNode node = fine_nodes[idx];

						const arrayType<9> pop = regularization_mass(node);

						dfloat node_mass_in = 0.0;
						dfloat node_mass_out = 0.0;

						for (size_t i = 0; i < 9; ++i)
						{
							if (node.incomings[i])
							{
								node_mass_in += node.pop_in[i];
							}

							if (node.outgoings[i])
							{
								node_mass_out += pop.f[i];
							}
						}

						fine_mass_in += node_mass_in;
						fine_mass_out += node_mass_out;
					}
				}

				std::cout << "fine_mass_in total: " << fine_mass_in << std::endl;
				std::cout << "fine_mass_in avg: " << fine_mass_in / (NY_FINE * NX_FINE) << std::endl;

				std::cout << "fine_mass_out total: " << fine_mass_out << std::endl;
				std::cout << "fine_mass_out avg: " << fine_mass_out / (NY_FINE * NX_FINE) << std::endl;

				std::cout << "fine_mass_net total: " << fine_mass_out - fine_mass_in << std::endl
						  << std::endl;

				// ---------------------------------------------
			}
		}

		if (calculate_mass)
		{
			sys_mass_in += fine_mass_in;
			sys_mass_out += fine_mass_out;

			std::cout << "sys_mass_in total: " << sys_mass_in << std::endl;
			std::cout << "sys_mass_in avg: " << sys_mass_in / (NY_FINE * NX_FINE + NX_COARSE * NY_COARSE) << std::endl;

			std::cout << "sys_mass_out total: " << sys_mass_out << std::endl;
			std::cout << "sys_mass_out avg: " << sys_mass_out / (NY_FINE * NX_FINE + NX_COARSE * NY_COARSE) << std::endl;

			std::cout << "sys_mass_net total: " << sys_mass_out - sys_mass_in << std::endl
					  << std::endl;
		}

		if (step % MACR_SAVE == 0)
		{
			printf("\n--------------------------- Last Time Step %06zu ---------------------------\n", step);

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