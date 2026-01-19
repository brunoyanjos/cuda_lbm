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

__device__ inline void bilinear_moment_interpolation(dfloat x, dfloat y,
                                                     int x0, int y0,
                                                     LBMState state, dfloat &mxx, dfloat &myy)
{
    int x1 = x0 + 1;
    int y1 = y0 + 1;

    const dfloat xd = x - dfloat(x0);
    const dfloat yd = y - dfloat(y0);

    const dfloat radius_1 = sqrt((dfloat(x0) - xc) * (dfloat(x0) - xc) + (dfloat(y0) - yc) * (dfloat(y0) - yc));
    const dfloat radius_2 = sqrt((dfloat(x1) - xc) * (dfloat(x1) - xc) + (dfloat(y0) - yc) * (dfloat(y0) - yc));
    const dfloat radius_3 = sqrt((dfloat(x0) - xc) * (dfloat(x0) - xc) + (dfloat(y1) - yc) * (dfloat(y1) - yc));
    const dfloat radius_4 = sqrt((dfloat(x1) - xc) * (dfloat(x1) - xc) + (dfloat(y1) - yc) * (dfloat(y1) - yc));

    const dfloat cos_theta_1 = (dfloat(x0) - xc) / radius_1;
    const dfloat cos_theta_2 = (dfloat(x1) - xc) / radius_2;
    const dfloat cos_theta_3 = (dfloat(x0) - xc) / radius_3;
    const dfloat cos_theta_4 = (dfloat(x1) - xc) / radius_4;

    const dfloat sen_theta_1 = (dfloat(y0) - yc) / radius_1;
    const dfloat sen_theta_2 = (dfloat(y0) - yc) / radius_2;
    const dfloat sen_theta_3 = (dfloat(y1) - yc) / radius_3;
    const dfloat sen_theta_4 = (dfloat(y1) - yc) / radius_4;

    const dfloat sen_two_theta_1 = dfloat(2.0) * cos_theta_1 * sen_theta_1;
    const dfloat sen_two_theta_2 = dfloat(2.0) * cos_theta_2 * sen_theta_2;
    const dfloat sen_two_theta_3 = dfloat(2.0) * cos_theta_3 * sen_theta_3;
    const dfloat sen_two_theta_4 = dfloat(2.0) * cos_theta_4 * sen_theta_4;

    const dfloat mxx1 = state.d_mxx[idxBlockCoord(x0, y0)];
    const dfloat mxx2 = state.d_mxx[idxBlockCoord(x1, y0)];
    const dfloat mxx3 = state.d_mxx[idxBlockCoord(x0, y1)];
    const dfloat mxx4 = state.d_mxx[idxBlockCoord(x1, y1)];

    const dfloat myy1 = state.d_myy[idxBlockCoord(x0, y0)];
    const dfloat myy2 = state.d_myy[idxBlockCoord(x1, y0)];
    const dfloat myy3 = state.d_myy[idxBlockCoord(x0, y1)];
    const dfloat myy4 = state.d_myy[idxBlockCoord(x1, y1)];

    const dfloat mxy1 = state.d_mxy[idxBlockCoord(x0, y0)];
    const dfloat mxy2 = state.d_mxy[idxBlockCoord(x1, y0)];
    const dfloat mxy3 = state.d_mxy[idxBlockCoord(x0, y1)];
    const dfloat mxy4 = state.d_mxy[idxBlockCoord(x1, y1)];

    const dfloat mxx_p1 = mxx1 * cos_theta_1 * cos_theta_1 + myy1 * sen_theta_1 * sen_theta_1 + mxy1 * sen_two_theta_1;
    const dfloat myy_p1 = mxx1 * sen_theta_1 * sen_theta_1 + myy1 * cos_theta_1 * cos_theta_1 - mxy1 * sen_two_theta_1;

    const dfloat mxx_p2 = mxx2 * cos_theta_2 * cos_theta_2 + myy2 * sen_theta_2 * sen_theta_2 + mxy2 * sen_two_theta_2;
    const dfloat myy_p2 = mxx2 * sen_theta_2 * sen_theta_2 + myy2 * cos_theta_2 * cos_theta_2 - mxy2 * sen_two_theta_2;

    const dfloat mxx_p3 = mxx3 * cos_theta_3 * cos_theta_3 + myy3 * sen_theta_3 * sen_theta_3 + mxy3 * sen_two_theta_3;
    const dfloat myy_p3 = mxx3 * sen_theta_3 * sen_theta_3 + myy3 * cos_theta_3 * cos_theta_3 - mxy3 * sen_two_theta_3;

    const dfloat mxx_p4 = mxx4 * cos_theta_4 * cos_theta_4 + myy4 * sen_theta_4 * sen_theta_4 + mxy4 * sen_two_theta_4;
    const dfloat myy_p4 = mxx4 * sen_theta_4 * sen_theta_4 + myy4 * cos_theta_4 * cos_theta_4 - mxy4 * sen_two_theta_4;

    const dfloat mxx_temp = (1.0 - xd) * mxx_p1 + xd * mxx_p2;
    const dfloat mxx_temp2 = (1.0 - xd) * mxx_p3 + xd * mxx_p4;

    const dfloat myy_temp = (1.0 - xd) * myy_p1 + xd * myy_p2;
    const dfloat myy_temp2 = (1.0 - xd) * myy_p3 + xd * myy_p4;

    mxx = (1.0 - yd) * mxx_temp + yd * mxx_temp2;
    myy = (1.0 - yd) * myy_temp + yd * myy_temp2;
}