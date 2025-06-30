
#include "main.cuh"
#include <iostream>
#include <chrono>
#include "saveData.cuh"

#include <string>
#include <math.h>
#include <cuda.h>
#include <curand.h>
#include <cuda_runtime.h>
#include <builtin_types.h>
#include "var.h"

using namespace std;

__global__ void printDeviceProperties(cylinderProperties *cylinder_proporties, unsigned int *dNodeType)
{
	cylinderProperties property = cylinder_proporties[threadIdx.x];
	int xb = property.xb;
	int yb = property.yb;

	int nodeType = dNodeType[idxScalarBlock(xb % BLOCK_NX, yb % BLOCK_NY, xb / BLOCK_NX, yb / BLOCK_NY)];

	printf("nodeType: %d\n"
		   "Is = {%d, %d, %d, %d, %d, %d, %d, %d, %d}\n"
		   "Os = {%d, %d, %d, %d, %d, %d, %d, %d, %d}\n",
		   nodeType,
		   property.is[0], property.is[1], property.is[2], property.is[3], property.is[4], property.is[5], property.is[6], property.is[7], property.is[8],
		   property.os[0], property.os[1], property.os[2], property.os[3], property.os[4], property.os[5], property.os[6], property.os[7], property.os[8]);
}

int main()
{
	printf("BLOCK_NX: %d, BLOCK_NY: %d\n", BLOCK_NX, BLOCK_NY);

	folderSetup();

	// set cuda device
	checkCudaErrors(cudaSetDevice(GPU_INDEX));

	// variable declaration
	dfloat *d_fMom;
	ghostInterfaceData ghostInterface;
	cylinderProperties *d_cylinder_properties;
	cylinderProperties *h_cylinder_properties;

	unsigned int *dNodeType;
	unsigned int *hNodeType;

	dfloat D_Max;
	size_t countor_count; 

	dfloat *h_fMom;

	dfloat *rho;

	dfloat *ux;
	dfloat *uy;

	/* ----------------- GRID AND THREADS DEFINITION FOR LBM ---------------- */
	dim3 threadBlock(BLOCK_NX, BLOCK_NY);
	dim3 gridBlock(NUM_BLOCK_X, NUM_BLOCK_Y);

	/* ------------------------- ALLOCATION FOR CPU ------------------------- */
	int step = 0;

	allocateHostMemory(&h_fMom, &rho, &ux, &uy);

	/* -------------- ALLOCATION FOR GPU ------------- */
	allocateDeviceMemory(&d_fMom, &dNodeType, &ghostInterface);

	// Setup Streams
	cudaStream_t streamsLBM[1];
	checkCudaErrors(cudaSetDevice(GPU_INDEX));
	checkCudaErrors(cudaStreamCreate(&streamsLBM[0]));
	checkCudaErrors(cudaDeviceSynchronize());

	initializeDomain(ghostInterface, d_fMom, h_fMom, hNodeType, dNodeType,
					 &step, gridBlock, threadBlock,
#ifdef CYLINDER
					 &D_Max, &h_cylinder_properties,
					 d_cylinder_properties, &countor_count
#endif
	);

	printf("count: %zu, d_max:%f\n", countor_count, D_Max);

	const dfloat VISC = U_MAX * D_Max / RE;
	const dfloat TAU = 0.5 + 3.0 * VISC; // relaxation time
	const dfloat OMEGA = 1.0 / TAU;		 // (tau)^-1

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
		streamingAndMom<<<gridBlock, threadBlock>>>(d_fMom, OMEGA, countor_count, dNodeType, ghostInterface, d_cylinder_properties, step);
		checkCudaErrors(cudaDeviceSynchronize());
		updateInnerBoundaries<<<1, countor_count>>>(d_fMom, d_cylinder_properties, OMEGA, step);

		checkCudaErrors(cudaDeviceSynchronize());
		boundaryAndCollision<<<gridBlock, threadBlock>>>(d_fMom, countor_count, OMEGA, dNodeType, ghostInterface, d_cylinder_properties, step);
#else
		gpuMomCollisionStream<<<gridBlock, threadBlock>>>(d_fMom, dNodeType, ghostInterface, step);
#endif

		// swap interface pointers
		swapGhostInterfaces(ghostInterface);

#ifdef CYLINDER
		if (step >= STAT_BEGIN_TIME && step <= STAT_END_TIME)
		{
			checkCudaErrors(cudaDeviceSynchronize());
			checkCudaErrors(cudaMemcpy(h_cylinder_properties, d_cylinder_properties, sizeof(cylinderProperties) * countor_count, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(h_fMom, d_fMom, sizeof(dfloat) * NUMBER_LBM_NODES * NUMBER_MOMENTS, cudaMemcpyDeviceToHost));

			calculate_forces(h_cylinder_properties, countor_count, step);
			calculate_pressure(h_cylinder_properties, countor_count, step);
			calculate_inlet_density(h_fMom, step);
		}
#endif // CYLINDER



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
	interfaceFree(ghostInterface);
	return 0;
}