#ifndef COLLISION_AND_RECONSTRUCTION
#define COLLISION_AND_RECONSTRUCTION

// CUDA INCLUDE
#include <cuda.h>
#include <cuda_runtime.h>
#include "device_launch_parameters.h"

#include "../var.h"

namespace second
{
	__device__ inline void pop_reconstruction(dfloat rhoVar, dfloat ux, dfloat uy, dfloat mxx, dfloat myy, dfloat mxy, dfloat *pop)
	{
		dfloat pics2 = dfloat(1.0) - cs2 * (mxx + myy);

		dfloat multiplyTerm = W0 * rhoVar;
		pop[0] = multiplyTerm * pics2 - W0;

		multiplyTerm = W1 * rhoVar;
		pop[1] = multiplyTerm * (pics2 + ux + mxx) - W1;
		pop[2] = multiplyTerm * (pics2 + uy + myy) - W1;
		pop[3] = multiplyTerm * (pics2 - ux + mxx) - W1;
		pop[4] = multiplyTerm * (pics2 - uy + myy) - W1;

		multiplyTerm = W2 * rhoVar;
		pop[5] = multiplyTerm * (pics2 + ux + uy + mxx + myy + mxy) - W2;
		pop[6] = multiplyTerm * (pics2 - ux + uy + mxx + myy - mxy) - W2;
		pop[7] = multiplyTerm * (pics2 - ux - uy + mxx + myy + mxy) - W2;
		pop[8] = multiplyTerm * (pics2 + ux - uy + mxx + myy - mxy) - W2;
	}
}

namespace third
{
	__device__ inline void pop_reconstruction(dfloat rhoVar, dfloat ux, dfloat uy, dfloat mxx, dfloat myy, dfloat mxy, dfloat *pop)
	{
		dfloat mxxy = ux * mxy + uy * mxx;
		dfloat mxyy = uy * mxy + ux * myy;

		dfloat one_cs2 = dfloat(1.0) - cs2;
		dfloat cs2_mxxy = cs2 * mxxy;
		dfloat cs2_mxyy = cs2 * mxyy;

		dfloat pics2 = dfloat(1.0) - cs2 * (mxx + myy);

		dfloat multiplyTerm = W0 * rhoVar;
		pop[0] = multiplyTerm * (pics2);

		multiplyTerm = W1 * rhoVar;
		pop[1] = multiplyTerm * (pics2 + ux + mxx - cs2_mxyy);
		pop[2] = multiplyTerm * (pics2 + uy + myy - cs2_mxxy);
		pop[3] = multiplyTerm * (pics2 - ux + mxx + cs2_mxyy);
		pop[4] = multiplyTerm * (pics2 - uy + myy + cs2_mxxy);

		multiplyTerm = W2 * rhoVar;
		pop[5] = multiplyTerm * (pics2 + ux + uy + mxx + myy + mxy + one_cs2 * mxxy + one_cs2 * mxyy);
		pop[6] = multiplyTerm * (pics2 - ux + uy + mxx + myy - mxy + one_cs2 * mxxy - one_cs2 * mxyy);
		pop[7] = multiplyTerm * (pics2 - ux - uy + mxx + myy + mxy - one_cs2 * mxxy - one_cs2 * mxyy);
		pop[8] = multiplyTerm * (pics2 + ux - uy + mxx + myy - mxy - one_cs2 * mxxy + one_cs2 * mxyy);
	}
}

namespace fourth
{
	__device__ inline void pop_reconstruction(dfloat rhoVar, dfloat ux, dfloat uy, dfloat mxx, dfloat myy, dfloat mxy, dfloat *pop)
	{
		dfloat mxxy = ux * mxy + uy * mxx;
		dfloat mxyy = uy * mxy + ux * myy;

		dfloat half = dfloat(0.5);

		dfloat mxxyy = half * (uy * mxxy + ux * ux * myy + ux * uy * mxy);

		dfloat one_cs2 = dfloat(1.0) - cs2;
		dfloat cs2_mxxy = cs2 * mxxy;
		dfloat cs2_mxyy = cs2 * mxyy;

		dfloat one_cs2_onde_cs2 = one_cs2 * one_cs2;

		dfloat cs2_mxxyy = cs2 * mxxyy;

		dfloat pics2 = dfloat(1.0) - cs2 * (mxx + myy);

		dfloat multiplyTerm = W0 * rhoVar;
		pop[0] = multiplyTerm * (pics2);

		multiplyTerm = W1 * rhoVar;
		pop[1] = multiplyTerm * (pics2 + ux + mxx - cs2_mxyy - one_cs2 * cs2_mxxyy);
		pop[2] = multiplyTerm * (pics2 + uy + myy - cs2_mxxy - one_cs2 * cs2_mxxyy);
		pop[3] = multiplyTerm * (pics2 - ux + mxx + cs2_mxyy - one_cs2 * cs2_mxxyy);
		pop[4] = multiplyTerm * (pics2 - uy + myy + cs2_mxxy - one_cs2 * cs2_mxxyy);

		multiplyTerm = W2 * rhoVar;
		pop[5] = multiplyTerm * (pics2 + ux + uy + mxx + myy + mxy + one_cs2 * mxxy + one_cs2 * mxyy + one_cs2_onde_cs2 * mxxyy);
		pop[6] = multiplyTerm * (pics2 - ux + uy + mxx + myy - mxy + one_cs2 * mxxy - one_cs2 * mxyy + one_cs2_onde_cs2 * mxxyy);
		pop[7] = multiplyTerm * (pics2 - ux - uy + mxx + myy + mxy - one_cs2 * mxxy - one_cs2 * mxyy + one_cs2_onde_cs2 * mxxyy);
		pop[8] = multiplyTerm * (pics2 + ux - uy + mxx + myy - mxy - one_cs2 * mxxy + one_cs2 * mxyy + one_cs2_onde_cs2 * mxxyy);
	}
}

__device__ inline void moment_collision(dfloat ux, dfloat uy, dfloat &mxx, dfloat &mxy, dfloat &myy, dfloat OMEGA)
{
	const dfloat omegaVar = OMEGA;
	const dfloat t_omegaVar = static_cast<dfloat>(1.0) - omegaVar;
	const dfloat omegaVar_d2 = omegaVar * static_cast<dfloat>(0.5);

	mxx = (t_omegaVar * mxx + omegaVar_d2 * ux * ux);
	myy = (t_omegaVar * myy + omegaVar_d2 * uy * uy);

	mxy = (t_omegaVar * mxy + omegaVar * ux * uy);
}

#endif // COLLISION_AND_RECONSTRUCTION
