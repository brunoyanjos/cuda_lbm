#ifndef INTERPOLATION_UTILITIES_CUH
#define INTERPOLATION_UTILITIES_CUH

#include "../../var.h"
#include "../../globalFunctions.h"

__device__
inline void bilinear_velocity_interpolation(dfloat x, dfloat y, int x0, int y0, int x1, int y1, dfloat* moms, int mom_index, dfloat* u_f) {
	const dfloat xd = x - x0;
	const dfloat yd = y - y0;

	const dfloat ux1 = moms[idxMom(x0 % BLOCK_NX, y0 % BLOCK_NY, mom_index, x0 / BLOCK_NX, y0 / BLOCK_NY)];
	const dfloat ux2 = moms[idxMom(x1 % BLOCK_NX, y0 % BLOCK_NY, mom_index, x1 / BLOCK_NX, y0 / BLOCK_NY)];
	const dfloat ux3 = moms[idxMom(x0 % BLOCK_NX, y1 % BLOCK_NY, mom_index, x0 / BLOCK_NX, y1 / BLOCK_NY)];
	const dfloat ux4 = moms[idxMom(x1 % BLOCK_NX, y1 % BLOCK_NY, mom_index, x1 / BLOCK_NX, y1 / BLOCK_NY)];

	const dfloat ux_y0 = (1.0 - xd) * ux1 + xd * ux2;
	const dfloat ux_y1 = (1.0 - xd) * ux3 + xd * ux4;

	*u_f = (1.0 - yd) * ux_y0 + yd * ux_y1;
}

__device__
inline void bilinear_moment_interpolation(dfloat x, dfloat y, int x0, int y0, int x1, int y1, dfloat* moms, dfloat* mxx_f, dfloat* myy_f) {

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

	const dfloat sen_two_theta_1 = 2.0 * cos_theta_1 * sen_theta_1;
	const dfloat sen_two_theta_2 = 2.0 * cos_theta_2 * sen_theta_2;
	const dfloat sen_two_theta_3 = 2.0 * cos_theta_3 * sen_theta_3;
	const dfloat sen_two_theta_4 = 2.0 * cos_theta_4 * sen_theta_4;

	const dfloat mxx1 = moms[idxMom(x0 % BLOCK_NX, y0 % BLOCK_NY, M_MXX_INDEX, x0 / BLOCK_NX, y0 / BLOCK_NY)];
	const dfloat mxx2 = moms[idxMom(x1 % BLOCK_NX, y0 % BLOCK_NY, M_MXX_INDEX, x1 / BLOCK_NX, y0 / BLOCK_NY)];
	const dfloat mxx3 = moms[idxMom(x0 % BLOCK_NX, y1 % BLOCK_NY, M_MXX_INDEX, x0 / BLOCK_NX, y1 / BLOCK_NY)];
	const dfloat mxx4 = moms[idxMom(x1 % BLOCK_NX, y1 % BLOCK_NY, M_MXX_INDEX, x1 / BLOCK_NX, y1 / BLOCK_NY)];

	const dfloat myy1 = moms[idxMom(x0 % BLOCK_NX, y0 % BLOCK_NY, M_MYY_INDEX, x0 / BLOCK_NX, y0 / BLOCK_NY)];
	const dfloat myy2 = moms[idxMom(x1 % BLOCK_NX, y0 % BLOCK_NY, M_MYY_INDEX, x1 / BLOCK_NX, y0 / BLOCK_NY)];
	const dfloat myy3 = moms[idxMom(x0 % BLOCK_NX, y1 % BLOCK_NY, M_MYY_INDEX, x0 / BLOCK_NX, y1 / BLOCK_NY)];
	const dfloat myy4 = moms[idxMom(x1 % BLOCK_NX, y1 % BLOCK_NY, M_MYY_INDEX, x1 / BLOCK_NX, y1 / BLOCK_NY)];

	const dfloat mxy1 = moms[idxMom(x0 % BLOCK_NX, y0 % BLOCK_NY, M_MXY_INDEX, x0 / BLOCK_NX, y0 / BLOCK_NY)];
	const dfloat mxy2 = moms[idxMom(x1 % BLOCK_NX, y0 % BLOCK_NY, M_MXY_INDEX, x1 / BLOCK_NX, y0 / BLOCK_NY)];
	const dfloat mxy3 = moms[idxMom(x0 % BLOCK_NX, y1 % BLOCK_NY, M_MXY_INDEX, x0 / BLOCK_NX, y1 / BLOCK_NY)];
	const dfloat mxy4 = moms[idxMom(x1 % BLOCK_NX, y1 % BLOCK_NY, M_MXY_INDEX, x1 / BLOCK_NX, y1 / BLOCK_NY)];

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

	*mxx_f = (1.0 - yd) * mxx_temp + yd * mxx_temp2;
	*myy_f = (1.0 - yd) * myy_temp + yd * myy_temp2;
}

__device__
inline dfloat extrapolation(dfloat delta, dfloat value1, dfloat value2) {
	const dfloat deltax = sqrt(2.0);

	return (delta * (delta - 2.0 * deltax) / (deltax * deltax)) * value1 - (delta * (delta - deltax) / (2.0 * deltax * deltax)) * value2;
}

#endif