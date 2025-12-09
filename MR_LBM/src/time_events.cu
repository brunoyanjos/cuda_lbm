#include "time_events.cuh"
#include "errorDef.h"

void initializeCudaEvents(cudaEvent_t &start, cudaEvent_t &stop, cudaEvent_t &start_step, cudaEvent_t &stop_step)
{
    checkCudaErrors(cudaSetDevice(GPU_INDEX));
    checkCudaErrors(cudaEventCreate(&start));
    checkCudaErrors(cudaEventCreate(&stop));
    checkCudaErrors(cudaEventCreate(&start_step));
    checkCudaErrors(cudaEventCreate(&stop_step));

    checkCudaErrors(cudaEventRecord(start, 0));
    checkCudaErrors(cudaEventRecord(start_step, 0));
}

dfloat recordElapsedTime(cudaEvent_t &start_step, cudaEvent_t &stop_step, int step)
{
    checkCudaErrors(cudaEventRecord(stop_step, 0));
    checkCudaErrors(cudaEventSynchronize(stop_step));

    float elapsedTime;
    checkCudaErrors(cudaEventElapsedTime(&elapsedTime, start_step, stop_step));
    elapsedTime *= 0.001;

    size_t nodesUpdatedSync = step * NUMBER_LBM_NODES;
    dfloat MLUPS = (nodesUpdatedSync / 1e6) / elapsedTime;
    return MLUPS;
}