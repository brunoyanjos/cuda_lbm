#include "state.cuh"
#include "../errorDef.h"

LBMState init_state()
{
    LBMState state;

    // Defining Variable Sizes
    state.bytes_fields = NUMBER_LBM_NODES * sizeof(dfloat);
    state.bytes_types = NUMBER_LBM_NODES * sizeof(uint8_t);

    // Allocating HOST memory
    state.h_node_type = (uint8_t *)malloc(state.bytes_types);
    state.h_rho = (dfloat *)malloc(state.bytes_fields);

    state.h_ux = (dfloat *)malloc(state.bytes_fields);
    state.h_uy = (dfloat *)malloc(state.bytes_fields);

    state.h_mxx = (dfloat *)malloc(state.bytes_fields);
    state.h_mxy = (dfloat *)malloc(state.bytes_fields);
    state.h_myy = (dfloat *)malloc(state.bytes_fields);

    // Allocating DEVICE mesmory
    checkCudaErrors(cudaMalloc(&state.d_node_type, state.bytes_types));
    checkCudaErrors(cudaMalloc(&state.d_rho, state.bytes_fields));

    checkCudaErrors(cudaMalloc(&state.d_ux, state.bytes_fields));
    checkCudaErrors(cudaMalloc(&state.d_uy, state.bytes_fields));

    checkCudaErrors(cudaMalloc(&state.d_mxx, state.bytes_fields));
    checkCudaErrors(cudaMalloc(&state.d_mxy, state.bytes_fields));
    checkCudaErrors(cudaMalloc(&state.d_myy, state.bytes_fields));

    return state;
}

void upload_state_to_host(LBMState &state)
{
    checkCudaErrors(cudaMemcpy(state.h_rho, state.d_rho, state.bytes_fields, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(state.h_ux, state.d_ux, state.bytes_fields, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(state.h_uy, state.d_uy, state.bytes_fields, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(state.h_mxx, state.d_mxx, state.bytes_fields, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(state.h_mxy, state.d_mxy, state.bytes_fields, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(state.h_myy, state.d_myy, state.bytes_fields, cudaMemcpyDeviceToHost));
}

void free_state(LBMState &state)
{
    free(state.h_node_type);
    free(state.h_rho);
    free(state.h_ux);
    free(state.h_uy);
    free(state.h_mxx);
    free(state.h_mxy);
    free(state.h_myy);

    checkCudaErrors(cudaFree(state.d_node_type));
    checkCudaErrors(cudaFree(state.d_rho));
    checkCudaErrors(cudaFree(state.d_ux));
    checkCudaErrors(cudaFree(state.d_uy));
    checkCudaErrors(cudaFree(state.d_mxx));
    checkCudaErrors(cudaFree(state.d_mxy));
    checkCudaErrors(cudaFree(state.d_myy));
}