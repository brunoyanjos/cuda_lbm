#pragma once

#include "var.h"
#include "globalFunctions.h"

__device__ inline dfloat bilinear_velocity_interpolation(
    dfloat x, dfloat y,
    int x0, int y0,
    const dfloat *__restrict__ u_array)
{
    int x1 = x0 + 1;
    int y1 = y0 + 1;

    const dfloat u1 = u_array[idxBlockCoord(x0, y0)];
    const dfloat u2 = u_array[idxBlockCoord(x1, y0)];
    const dfloat u3 = u_array[idxBlockCoord(x0, y1)];
    const dfloat u4 = u_array[idxBlockCoord(x1, y1)];

    const dfloat xd = x - x0;
    const dfloat yd = y - y0;

    const dfloat u_y0 = fma(xd, u2 - u1, u1);
    const dfloat u_y1 = fma(xd, u4 - u3, u3);

    return fma(yd, u_y1 - u_y0, u_y0);
}