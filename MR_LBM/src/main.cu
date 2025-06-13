
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
	dfloat *d_fMom;
	dfloat *d_fMom_old;
	ghostInterfaceData ghostInterface;
	cylinderProperties *d_cylinder_proporties;
	cylinderProperties *h_cylinder_proporties;

	unsigned int *dNodeType;
	unsigned int *hNodeType;

	dfloat D_Max;
	size_t countor_count;

	dfloat *h_fMom;

	dfloat *rho;

	dfloat *ux;
	dfloat *uy;

	dfloat *mxx;
	dfloat *mxy;
	dfloat *myy;

	/* ----------------- GRID AND THREADS DEFINITION FOR LBM ---------------- */
	dim3 threadBlock(BLOCK_NX, BLOCK_NY);
	dim3 gridBlock(NUM_BLOCK_X, NUM_BLOCK_Y);

	/* ------------------------- ALLOCATION FOR CPU ------------------------- */
	int step = 0;

	allocateHostMemory(&h_fMom, &rho, &ux, &uy, &mxx, &mxy, &myy);

	/* -------------- ALLOCATION FOR GPU ------------- */
	allocateDeviceMemory(&d_fMom, &dNodeType, &ghostInterface);

	cudaMalloc(&d_fMom_old, MEM_SIZE_MOM);

	// Setup Streams
	cudaStream_t streamsLBM[1];
	checkCudaErrors(cudaSetDevice(GPU_INDEX));
	checkCudaErrors(cudaStreamCreate(&streamsLBM[0]));
	checkCudaErrors(cudaDeviceSynchronize());

	initializeDomain(ghostInterface, d_fMom, h_fMom, hNodeType, dNodeType,
					 &step, gridBlock, threadBlock,
#ifdef CYLINDER
					 &D_Max, &h_cylinder_proporties,
					 d_cylinder_proporties, &countor_count
#endif
	);

	printf("count: %d, d_max:%f", countor_count, D_Max);

	const dfloat VISC = U_MAX * D_Max / RE;
	const dfloat TAU = 0.5 + 3.0 * VISC; // relaxation time

	const dfloat OMEGA = 1.0 / TAU;			   // (tau)^-1
	const dfloat OMEGAd2 = OMEGA / 2.0;		   // OMEGA/2
	const dfloat OMEGAd9 = OMEGA / 9.0;		   // OMEGA/9
	const dfloat T_OMEGA = 1.0 - OMEGA;		   // 1-OMEGA
	const dfloat TT_OMEGA = 1.0 - 0.5 * OMEGA; // 1.0 - OMEGA/2
	const dfloat OMEGA_P1 = 1.0 + OMEGA;	   // 1+ OMEGA
	const dfloat TT_OMEGA_T3 = TT_OMEGA * 3.0; // 3*(1-0.5*OMEGA)

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
		#ifdef CYLINDER
		streamingAndMom << <gridBlock, threadBlock >> > (d_fMom, OMEGA, dNodeType, ghostInterface, d_cylinder_proporties, step);
		checkCudaErrors(cudaDeviceSynchronize());
		checkCudaErrors(cudaMemcpy(d_fMom_old, d_fMom, sizeof(dfloat) * NUMBER_LBM_NODES * NUMBER_MOMENTS, cudaMemcpyDeviceToDevice));

		boundaryAndCollision << <gridBlock, threadBlock >> > (d_fMom, d_fMom_old, OMEGA, dNodeType, ghostInterface, d_cylinder_proporties, step);
		#else
		gpuMomCollisionStream<<<gridBlock, threadBlock>>>(d_fMom, dNodeType, ghostInterface, step);
		#endif

		// swap interface pointers
		swapGhostInterfaces(ghostInterface);

		if (MACR_SAVE != 0 && step % MACR_SAVE == 0)
		{
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