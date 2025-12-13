#pragma once

#include "var.h"
#include "globalFunctions.h"
#include "init/state.cuh"

__device__ inline void bilinear_density_interpolation(dfloat x, dfloat y, int x0, int y0, int x1, int y1, dfloat *moms, int mom_index, dfloat *r_f)
{
    const dfloat xd = x - x0;
    const dfloat yd = y - y0;

    const dfloat r1 = RHO_0 + moms[idxMom(x0 % BLOCK_NX, y0 % BLOCK_NY, mom_index, x0 / BLOCK_NX, y0 / BLOCK_NY)];
    const dfloat r2 = RHO_0 + moms[idxMom(x1 % BLOCK_NX, y0 % BLOCK_NY, mom_index, x1 / BLOCK_NX, y0 / BLOCK_NY)];
    const dfloat r3 = RHO_0 + moms[idxMom(x0 % BLOCK_NX, y1 % BLOCK_NY, mom_index, x0 / BLOCK_NX, y1 / BLOCK_NY)];
    const dfloat r4 = RHO_0 + moms[idxMom(x1 % BLOCK_NX, y1 % BLOCK_NY, mom_index, x1 / BLOCK_NX, y1 / BLOCK_NY)];

    const dfloat r_y0 = (static_cast<dfloat>(1.0) - xd) * r1 + xd * r2;
    const dfloat r_y1 = (static_cast<dfloat>(1.0) - xd) * r3 + xd * r4;

    *r_f = (static_cast<dfloat>(1.0) - yd) * r_y0 + yd * r_y1;
}

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

__device__ inline void bilinear_moment_interpolation(dfloat x, dfloat y, int x0, int y0, LBMState state, dfloat &mxx_f, dfloat &myy_f)
{
    int x1 = x0 + 1;
    int y1 = y0 + 1;

    const dfloat xd = x - dfloat(x0);
    const dfloat yd = y - dfloat(y0);

    const dfloat x0f = dfloat(x0);
    const dfloat x1f = dfloat(x1);
    const dfloat y0f = dfloat(y0);
    const dfloat y1f = dfloat(y1);

    const dfloat dx0 = x0f - xc;
    const dfloat dx1 = x1f - xc;
    const dfloat dy0 = y0f - yc;
    const dfloat dy1 = y1f - yc;

    const dfloat r1_2 = dx0 * dx0 + dy0 * dy0;
    const dfloat r2_2 = dx1 * dx1 + dy0 * dy0;
    const dfloat r3_2 = dx0 * dx0 + dy1 * dy1;
    const dfloat r4_2 = dx1 * dx1 + dy1 * dy1;

    const dfloat inv_r1 = dfloat(1.0) / dsqrt(r1_2);
    const dfloat inv_r2 = dfloat(1.0) / dsqrt(r2_2);
    const dfloat inv_r3 = dfloat(1.0) / dsqrt(r3_2);
    const dfloat inv_r4 = dfloat(1.0) / dsqrt(r4_2);

    const dfloat cos_theta_1 = dx0 * inv_r1;
    const dfloat cos_theta_2 = dx1 * inv_r2;
    const dfloat cos_theta_3 = dx0 * inv_r3;
    const dfloat cos_theta_4 = dx1 * inv_r4;

    const dfloat sen_theta_1 = dy0 * inv_r1;
    const dfloat sen_theta_2 = dy0 * inv_r2;
    const dfloat sen_theta_3 = dy1 * inv_r3;
    const dfloat sen_theta_4 = dy1 * inv_r4;

    const dfloat sen_two_theta_1 = static_cast<dfloat>(2.0) * cos_theta_1 * sen_theta_1;
    const dfloat sen_two_theta_2 = static_cast<dfloat>(2.0) * cos_theta_2 * sen_theta_2;
    const dfloat sen_two_theta_3 = static_cast<dfloat>(2.0) * cos_theta_3 * sen_theta_3;
    const dfloat sen_two_theta_4 = static_cast<dfloat>(2.0) * cos_theta_4 * sen_theta_4;

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

    const dfloat mxx_y0 = fma(xd, mxx_p2 - mxx_p1, mxx_p1);
    const dfloat mxx_y1 = fma(xd, mxx_p4 - mxx_p3, mxx_p3);

    const dfloat myy_y0 = fma(xd, myy_p2 - myy_p1, myy_p1);
    const dfloat myy_y1 = fma(xd, myy_p4 - myy_p3, myy_p3);

    mxx_f = fma(yd, mxx_y1 - mxx_y0, mxx_y0);
    myy_f = fma(yd, myy_y1 - myy_y0, myy_y0);
}

__device__ inline void pressure_extrapolation_old(dfloat xw, dfloat yw, dfloat x1, dfloat y1, dfloat x2, dfloat y2, dfloat x3, dfloat y3, dfloat rho1, dfloat rho2, dfloat rho3, dfloat *pressure)
{

    // pressure interpolation
    dfloat xw_diff = xw - xc;
    dfloat yw_diff = yw - yc;

    dfloat x1_diff = x1 - xc;
    dfloat y1_diff = y1 - yc;

    dfloat x2_diff = x2 - xc;
    dfloat y2_diff = y2 - yc;

    dfloat x3_diff = x3 - xc;
    dfloat y3_diff = y3 - yc;

    dfloat rw2 = xw_diff * xw_diff + yw_diff * yw_diff;
    dfloat r12 = x1_diff * x1_diff + y1_diff * y1_diff;
    dfloat r22 = x2_diff * x2_diff + y2_diff * y2_diff;
    dfloat r32 = x3_diff * x3_diff + y3_diff * y3_diff;

    dfloat rw = dsqrt(rw2);
    dfloat r1 = dsqrt(r12);
    dfloat r2 = dsqrt(r22);
    dfloat r3 = dsqrt(r32);

    dfloat denom = (r1 - r2) * (r1 - r3) * (r2 - r3);

    dfloat p1 = rho1 * cs2;
    dfloat p2 = rho2 * cs2;
    dfloat p3 = rho3 * cs2;

    dfloat a0 = (r1 * r3 * p2 * (r3 - r1) + (r2 * r2) * (r3 * p1 - r1 * p3) + r2 * ((r1 * r1) * p3 - (r3 * r3) * p1)) / denom;
    dfloat a1 = ((r3 * r3) * (p1 - p2) + (r1 * r1) * (p2 - p3) + (r2 * r2) * (p3 - p1)) / denom;
    dfloat a2 = (r3 * (p2 - p1) + r2 * (p1 - p3) + r1 * (p3 - p2)) / denom;

    *pressure = a0 + a1 * rw + a2 * (rw * rw);
}

__device__ inline void pressure_extrapolation(dfloat rhob, dfloat rho1, dfloat rho2, dfloat *pressure, dfloat delta)
{
    dfloat pb = rhob * cs2;

    dfloat p1 = rho1 * cs2;
    dfloat p2 = rho2 * cs2;

    const dfloat deltax = dsqrt(static_cast<dfloat>(2.0));

    const dfloat deltax2 = deltax * deltax;
    const dfloat inv_deltax2 = static_cast<dfloat>(1) / deltax2;

    const dfloat delta2 = delta * delta;

    const dfloat first_term = (static_cast<dfloat>(2) * deltax2 - delta2 + static_cast<dfloat>(3) * delta * deltax) * static_cast<dfloat>(0.5) * inv_deltax2;
    const dfloat second_term = delta * (delta - static_cast<dfloat>(2.0) * deltax) * inv_deltax2;
    const dfloat third_term = delta * (delta - deltax) * static_cast<dfloat>(0.5) * inv_deltax2;

    const dfloat inv_first = static_cast<dfloat>(1) / first_term;

    *pressure = (pb - second_term * p1 - third_term * p2) * inv_first;
    ;
}

__device__ inline dfloat extrapolation(dfloat delta, dfloat value1, dfloat value2)
{
    const dfloat deltax = dsqrt(static_cast<dfloat>(2.0));

    return (delta * (delta - static_cast<dfloat>(2.0) * deltax) / (deltax * deltax)) * value1 - (delta * (delta - deltax) / (static_cast<dfloat>(2.0) * deltax * deltax)) * value2;
}