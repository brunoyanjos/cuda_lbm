#pragma once

#include "../var.h"
#include "../globalStructs.h"

__host__ __device__ inline void interfaceSwap(dfloat *&pt1, dfloat *&pt2)
{
    dfloat *temp = pt1;
    pt1 = pt2;
    pt2 = temp;
}

__host__ inline void interfaceFree(ghostInterfaceData &ghostInterface)
{
    cudaFree(ghostInterface.fGhost.X_0);
    cudaFree(ghostInterface.fGhost.X_1);
    cudaFree(ghostInterface.fGhost.Y_0);
    cudaFree(ghostInterface.fGhost.Y_1);

    cudaFree(ghostInterface.gGhost.X_0);
    cudaFree(ghostInterface.gGhost.X_1);
    cudaFree(ghostInterface.gGhost.Y_0);
    cudaFree(ghostInterface.gGhost.Y_1);
}

__host__ inline void swapGhostInterfaces(GhostInterfaceData &ghostInterface)
{

    // Swap interface pointers for fGhost and gGhost
    interfaceSwap(ghostInterface.fGhost.X_0, ghostInterface.gGhost.X_0);
    interfaceSwap(ghostInterface.fGhost.X_1, ghostInterface.gGhost.X_1);
    interfaceSwap(ghostInterface.fGhost.Y_0, ghostInterface.gGhost.Y_0);
    interfaceSwap(ghostInterface.fGhost.Y_1, ghostInterface.gGhost.Y_1);
}

__host__ inline void interfaceMalloc(ghostInterfaceData &ghostInterface)
{
    cudaMalloc((void **)&(ghostInterface.fGhost.X_0), sizeof(dfloat) * NUMBER_GHOST_FACE_X * QF);
    cudaMalloc((void **)&(ghostInterface.fGhost.X_1), sizeof(dfloat) * NUMBER_GHOST_FACE_X * QF);
    cudaMalloc((void **)&(ghostInterface.fGhost.Y_0), sizeof(dfloat) * NUMBER_GHOST_FACE_Y * QF);
    cudaMalloc((void **)&(ghostInterface.fGhost.Y_1), sizeof(dfloat) * NUMBER_GHOST_FACE_Y * QF);

    cudaMalloc((void **)&(ghostInterface.gGhost.X_0), sizeof(dfloat) * NUMBER_GHOST_FACE_X * QF);
    cudaMalloc((void **)&(ghostInterface.gGhost.X_1), sizeof(dfloat) * NUMBER_GHOST_FACE_X * QF);
    cudaMalloc((void **)&(ghostInterface.gGhost.Y_0), sizeof(dfloat) * NUMBER_GHOST_FACE_Y * QF);
    cudaMalloc((void **)&(ghostInterface.gGhost.Y_1), sizeof(dfloat) * NUMBER_GHOST_FACE_Y * QF);
}
