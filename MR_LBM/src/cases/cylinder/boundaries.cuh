#ifndef BOUNDARIES_CUH
#define BOUNDARIES_CUH

// CUDA INCLUDE
#include <cuda.h>
#include <cuda_runtime.h>
#include "device_launch_parameters.h"

#include "../../var.h"
#include "../../globalFunctions.h"
#include "../../nodeTypeMap.h"
#include "numerical_solutions.cuh"
#include "aux_functions.cuh"
#include "interpolation_utilities.cuh"

__host__ __device__ inline void boundary_definition(unsigned int *nodeType, unsigned int x, unsigned int y)
{
	if (x == 0)
	{
		*nodeType = WEST;
	}
	else if (x == (NX - 1))
	{
		*nodeType = EAST;
	}
	else
	{
		*nodeType = BULK;
	}
}

__device__ inline void boundary_calculation(unsigned int nodeType, dfloat &rhoVar,
											dfloat &ux, dfloat &uy,
											dfloat &mxx, dfloat &myy, dfloat &mxy,
											dfloat *pop, dfloat *fMom, int x, int y, dfloat OMEGA)
{
	const dfloat pop_0 = pop[0] + W0;

	const dfloat pop_1 = pop[1] + W1;
	const dfloat pop_2 = pop[2] + W1;
	const dfloat pop_3 = pop[3] + W1;
	const dfloat pop_4 = pop[4] + W1;

	const dfloat pop_5 = pop[5] + W2;
	const dfloat pop_6 = pop[6] + W2;
	const dfloat pop_7 = pop[7] + W2;
	const dfloat pop_8 = pop[8] + W2;

	switch (nodeType)
	{
	case WEST:
	{
		const dfloat rhoIn = pop_0 + pop_2 + pop_3 + pop_4 + pop_6 + pop_7;
		const dfloat inv_rhoIn = static_cast<float>(1) / rhoIn;

		const dfloat mxxIn = (pop_3 + pop_6 + pop_7) * inv_rhoIn - cs2;
		const dfloat mxyIn = (pop_7 - pop_6) * inv_rhoIn;
		const dfloat myyIn = (pop_2 + pop_4 + pop_6 + pop_7) * inv_rhoIn - cs2;

		uy = static_cast<float>(0);
		ux = U_MAX;

		const dfloat rho = (static_cast<float>(4) * rhoIn + static_cast<float>(3) * rhoIn * mxxIn) /
						   (static_cast<float>(3) - static_cast<float>(3) * (ux));
		mxx = (rho + static_cast<float>(9) * rhoIn * mxxIn + static_cast<float>(3) * rho * (ux)) / (static_cast<float>(6) * rho);
		mxy = static_cast<float>(2) * rhoIn * mxyIn / rho;
		myy = static_cast<float>(6) * rhoIn * myyIn / (static_cast<float>(5) * rho);

		rhoVar = rho;

		break;
	}
	case EAST:
	{
		const dfloat rhoIn = pop_0 + pop_1 + pop_2 + pop_4 + pop_5 + pop_8;
		const dfloat inv_rhoIn = static_cast<float>(1) / rhoIn;

		const dfloat mxxIn = (pop_1 + pop_5 + pop_8) * inv_rhoIn - cs2;
		const dfloat mxyIn = (pop_5 - pop_8) * inv_rhoIn;
		const dfloat myyIn = (pop_2 + pop_4 + pop_5 + pop_8) * inv_rhoIn - cs2;

		const dfloat rho = RHO_0 + fMom[idxMom(threadIdx.x - 1, threadIdx.y, M_RHO_INDEX, blockIdx.x, blockIdx.y)];
		ux = fMom[idxMom(threadIdx.x - 1, threadIdx.y, M_UX_INDEX, blockIdx.x, blockIdx.y)] / F_M_I_SCALE;
		uy = fMom[idxMom(threadIdx.x - 1, threadIdx.y, M_UY_INDEX, blockIdx.x, blockIdx.y)] / F_M_I_SCALE;

		rhoVar = rho;

		mxx = (rho + static_cast<float>(9) * rhoIn * mxxIn - static_cast<float>(3) * rho * (ux)) / (static_cast<float>(6) * rho);
		mxy = (static_cast<float>(6) * rhoIn * mxyIn - rho * (uy)) / (static_cast<float>(3) * rho);
		myy = static_cast<float>(6) * rhoIn * myyIn / (static_cast<float>(5) * rho);

		break;
	}
	default:
		break;
	}
}

__device__ inline void fluid_boundary_evaluation(unsigned int nodeType,
												 cylinderProperties *cylinder_properties, size_t cylinder_counter,
												 const dfloat (&pop)[9], dfloat OMEGA,
												 dfloat &rhoVar, dfloat &ux_t30, dfloat &uy_t30,
												 dfloat &m_xx_t45, dfloat &m_yy_t45, dfloat &m_xy_t90,
												 size_t x, size_t y)
{
	dfloat rho_I = static_cast<dfloat>(0);

	dfloat ux_I = static_cast<dfloat>(0);
	dfloat uy_I = static_cast<dfloat>(0);

	dfloat m_xx_I = static_cast<dfloat>(0);
	dfloat m_yy_I = static_cast<dfloat>(0);
	dfloat m_xy_I = static_cast<dfloat>(0);

	cylinderProperties *bc_property = findCylindeProperty(cylinder_properties, cylinder_counter, x, y);

	for (size_t i = 0; i < 9; i++)
	{
		if ((*bc_property).is[i] == 1)
		{

			const dfloat Hxx = (cx[i] * cx[i]) - cs2;
			const dfloat Hyy = (cy[i] * cy[i]) - cs2;
			const dfloat Hxy = cx[i] * cy[i];

			rho_I += (pop[i] + w[i]);

			ux_I += (pop[i] + w[i]) * cx[i];
			uy_I += (pop[i] + w[i]) * cy[i];

			m_xx_I += (pop[i] + w[i]) * Hxx;
			m_yy_I += (pop[i] + w[i]) * Hyy;
			m_xy_I += (pop[i] + w[i]) * Hxy;
		}
	}

	if (nodeType == 201)
	{
		double linear_part = (5202.0 + 612.0 * m_xx_I - 1836.0 * m_xy_I + 612.0 * m_yy_I -
							  720.0 * m_xx_I * OMEGA + 2160.0 * m_xy_I * OMEGA - 720.0 * m_yy_I * OMEGA +
							  255.0 * ux_I + 159.0 * OMEGA * ux_I - 255.0 * uy_I - 159.0 * OMEGA * uy_I) *
							 rho_I;

		double inner_expr = (1734.0 + 204.0 * m_yy_I + m_xx_I * (204.0 - 240.0 * OMEGA) -
							 240.0 * m_yy_I * OMEGA + 36.0 * m_xy_I * (-17.0 + 20.0 * OMEGA) +
							 85.0 * ux_I + 53.0 * OMEGA * ux_I - 85.0 * uy_I - 53.0 * OMEGA * uy_I);

		double squared_term = 3.0 * inner_expr * inner_expr;

		double quadratic_terms = (45.0 * m_xx_I * m_xx_I + 405.0 * m_xy_I * m_xy_I + 45.0 * m_yy_I * m_yy_I -
								  345.0 * m_yy_I * ux_I + 589.0 * ux_I * ux_I -
								  15.0 * m_xx_I * (18.0 * m_xy_I - 6.0 * m_yy_I + 23.0 * ux_I - 23.0 * uy_I) +
								  345.0 * m_yy_I * uy_I - 1467.0 * ux_I * uy_I + 589.0 * uy_I * uy_I -
								  45.0 * m_xy_I * (6.0 * m_yy_I - 23.0 * ux_I + 23.0 * uy_I));

		double OMEGA_factor = 2.0 * OMEGA * (4998.0 + 103.0 * OMEGA);

		double sqrt_expr = std::sqrt(3.0 * rho_I * rho_I * (squared_term + OMEGA_factor * quadratic_terms));

		rhoVar = (linear_part + sqrt_expr) / (9996.0 + 206.0 * OMEGA);

		ux_t30 = -(3.0 * m_xx_I * rho_I - 9.0 * m_xy_I * rho_I + 3.0 * m_yy_I * rho_I -
				   20.0 * rho_I * ux_I + 3.0 * rho_I * uy_I + rhoVar) /
				 (17.0 * rhoVar);

		uy_t30 = -(-3.0 * m_xx_I * rho_I + 9.0 * m_xy_I * rho_I - 3.0 * m_yy_I * rho_I +
				   3.0 * rho_I * ux_I - 20.0 * rho_I * uy_I - rhoVar) /
				 (17.0 * rhoVar);

		m_xx_t45 = -(-57.0 * m_xx_I * rho_I + 18.0 * m_xy_I * rho_I - 6.0 * m_yy_I * rho_I +
					 6.0 * rho_I * ux_I - 6.0 * rho_I * uy_I - 2.0 * rhoVar) /
				   (51.0 * rhoVar);

		m_yy_t45 = (6.0 * m_xx_I * rho_I - 18.0 * m_xy_I * rho_I + 57.0 * m_yy_I * rho_I -
					6.0 * rho_I * ux_I + 6.0 * rho_I * uy_I + 2.0 * rhoVar) /
				   (51.0 * rhoVar);

		m_xy_t90 = -(3.0 * m_xx_I * rho_I - 26.0 * m_xy_I * rho_I + 3.0 * m_yy_I * rho_I -
					 3.0 * rho_I * ux_I + 3.0 * rho_I * uy_I + rhoVar) /
				   (17.0 * rhoVar);
	}
	else if (nodeType == 202)
	{
		// BCFLUID-TYPE2
		double linear_part = (5202.0 * rho_I + 612.0 * m_xx_I * rho_I + 1836.0 * m_xy_I * rho_I +
							  612.0 * m_yy_I * rho_I - 720.0 * m_xx_I * OMEGA * rho_I -
							  2160.0 * m_xy_I * OMEGA * rho_I - 720.0 * m_yy_I * OMEGA * rho_I -
							  255.0 * rho_I * ux_I - 159.0 * OMEGA * rho_I * ux_I -
							  255.0 * rho_I * uy_I - 159.0 * OMEGA * rho_I * uy_I);

		double inner_expr = (-1734.0 - 204.0 * m_yy_I + 240.0 * m_yy_I * OMEGA +
							 12.0 * m_xx_I * (-17.0 + 20.0 * OMEGA) +
							 36.0 * m_xy_I * (-17.0 + 20.0 * OMEGA) +
							 85.0 * ux_I + 53.0 * OMEGA * ux_I +
							 85.0 * uy_I + 53.0 * OMEGA * uy_I);

		double squared_term = 3.0 * inner_expr * inner_expr;

		double quadratic_terms = (45.0 * m_xx_I * m_xx_I + 405.0 * m_xy_I * m_xy_I + 45.0 * m_yy_I * m_yy_I +
								  345.0 * m_yy_I * ux_I + 589.0 * ux_I * ux_I +
								  345.0 * m_yy_I * uy_I + 1467.0 * ux_I * uy_I +
								  589.0 * uy_I * uy_I +
								  45.0 * m_xy_I * (6.0 * m_yy_I + 23.0 * (ux_I + uy_I)) +
								  15.0 * m_xx_I * (18.0 * m_xy_I + 6.0 * m_yy_I + 23.0 * (ux_I + uy_I)));

		double OMEGA_factor = 2.0 * OMEGA * (4998.0 + 103.0 * OMEGA);

		double sqrt_expr = std::sqrt(3.0 * rho_I * rho_I * (squared_term + OMEGA_factor * quadratic_terms));

		rhoVar = (linear_part + sqrt_expr) / (9996.0 + 206.0 * OMEGA);

		ux_t30 = -(-3.0 * m_xx_I * rho_I - 9.0 * m_xy_I * rho_I - 3.0 * m_yy_I * rho_I -
				   20.0 * rho_I * ux_I - 3.0 * rho_I * uy_I - rhoVar) /
				 (17.0 * rhoVar);

		uy_t30 = -(-3.0 * m_xx_I * rho_I - 9.0 * m_xy_I * rho_I - 3.0 * m_yy_I * rho_I -
				   3.0 * rho_I * ux_I - 20.0 * rho_I * uy_I - rhoVar) /
				 (17.0 * rhoVar);

		m_xx_t45 = -(-57.0 * m_xx_I * rho_I - 18.0 * m_xy_I * rho_I - 6.0 * m_yy_I * rho_I -
					 6.0 * rho_I * ux_I - 6.0 * rho_I * uy_I - 2.0 * rhoVar) /
				   (51.0 * rhoVar);

		m_yy_t45 = (6.0 * m_xx_I * rho_I + 18.0 * m_xy_I * rho_I + 57.0 * m_yy_I * rho_I +
					6.0 * rho_I * ux_I + 6.0 * rho_I * uy_I + 2.0 * rhoVar) /
				   (51.0 * rhoVar);

		m_xy_t90 = (3.0 * m_xx_I * rho_I + 26.0 * m_xy_I * rho_I + 3.0 * m_yy_I * rho_I +
					3.0 * rho_I * ux_I + 3.0 * rho_I * uy_I + rhoVar) /
				   (17.0 * rhoVar);
	}
	else if (nodeType == 203)
	{
		// BCFLUID-TYPE3
		double linear_part = (5202.0 * rho_I + 612.0 * m_xx_I * rho_I + 1836.0 * m_xy_I * rho_I +
							  612.0 * m_yy_I * rho_I - 720.0 * m_xx_I * OMEGA * rho_I -
							  2160.0 * m_xy_I * OMEGA * rho_I - 720.0 * m_yy_I * OMEGA * rho_I +
							  255.0 * rho_I * ux_I + 159.0 * OMEGA * rho_I * ux_I +
							  255.0 * rho_I * uy_I + 159.0 * OMEGA * rho_I * uy_I);

		double inner_expr = (1734.0 + 204.0 * m_xx_I + 612.0 * m_xy_I + 204.0 * m_yy_I -
							 240.0 * m_xx_I * OMEGA - 720.0 * m_xy_I * OMEGA - 240.0 * m_yy_I * OMEGA +
							 85.0 * ux_I + 53.0 * OMEGA * ux_I + 85.0 * uy_I + 53.0 * OMEGA * uy_I);

		double squared_term = 3.0 * inner_expr * inner_expr;

		double quadratic_terms = (45.0 * m_xx_I * m_xx_I + 405.0 * m_xy_I * m_xy_I + 45.0 * m_yy_I * m_yy_I -
								  345.0 * m_yy_I * ux_I + 589.0 * ux_I * ux_I -
								  345.0 * m_yy_I * uy_I + 1467.0 * ux_I * uy_I +
								  589.0 * uy_I * uy_I +
								  45.0 * m_xy_I * (6.0 * m_yy_I - 23.0 * (ux_I + uy_I)) +
								  15.0 * m_xx_I * (18.0 * m_xy_I + 6.0 * m_yy_I - 23.0 * (ux_I + uy_I)));

		double OMEGA_factor = 2.0 * OMEGA * (4998.0 + 103.0 * OMEGA);

		double sqrt_expr = std::sqrt(3.0 * rho_I * rho_I * (squared_term + OMEGA_factor * quadratic_terms));

		rhoVar = (linear_part + sqrt_expr) / (9996.0 + 206.0 * OMEGA);

		ux_t30 = -(3.0 * m_xx_I * rho_I + 9.0 * m_xy_I * rho_I + 3.0 * m_yy_I * rho_I -
				   20.0 * rho_I * ux_I - 3.0 * rho_I * uy_I + rhoVar) /
				 (17.0 * rhoVar);

		uy_t30 = -(3.0 * m_xx_I * rho_I + 9.0 * m_xy_I * rho_I + 3.0 * m_yy_I * rho_I -
				   3.0 * rho_I * ux_I - 20.0 * rho_I * uy_I + rhoVar) /
				 (17.0 * rhoVar);

		m_xx_t45 = -(-57.0 * m_xx_I * rho_I - 18.0 * m_xy_I * rho_I - 6.0 * m_yy_I * rho_I +
					 6.0 * rho_I * ux_I + 6.0 * rho_I * uy_I - 2.0 * rhoVar) /
				   (51.0 * rhoVar);

		m_yy_t45 = (6.0 * m_xx_I * rho_I + 18.0 * m_xy_I * rho_I + 57.0 * m_yy_I * rho_I -
					6.0 * rho_I * ux_I - 6.0 * rho_I * uy_I + 2.0 * rhoVar) /
				   (51.0 * rhoVar);

		m_xy_t90 = (3.0 * m_xx_I * rho_I + 26.0 * m_xy_I * rho_I + 3.0 * m_yy_I * rho_I -
					3.0 * rho_I * ux_I - 3.0 * rho_I * uy_I + rhoVar) /
				   (17.0 * rhoVar);
	}
	else
	{
		// BCFLUID-TYPE4
		double linear_part = (5202.0 * rho_I + 612.0 * m_xx_I * rho_I - 1836.0 * m_xy_I * rho_I +
							  612.0 * m_yy_I * rho_I - 720.0 * m_xx_I * OMEGA * rho_I +
							  2160.0 * m_xy_I * OMEGA * rho_I - 720.0 * m_yy_I * OMEGA * rho_I -
							  255.0 * rho_I * ux_I - 159.0 * OMEGA * rho_I * ux_I +
							  255.0 * rho_I * uy_I + 159.0 * OMEGA * rho_I * uy_I);

		double inner_expr = (1734.0 + 204.0 * m_yy_I + m_xx_I * (204.0 - 240.0 * OMEGA) -
							 240.0 * m_yy_I * OMEGA + 36.0 * m_xy_I * (-17.0 + 20.0 * OMEGA) -
							 85.0 * ux_I - 53.0 * OMEGA * ux_I + 85.0 * uy_I + 53.0 * OMEGA * uy_I);

		double squared_term = 3.0 * inner_expr * inner_expr;

		double quadratic_terms = (45.0 * m_xx_I * m_xx_I + 405.0 * m_xy_I * m_xy_I + 45.0 * m_yy_I * m_yy_I +
								  345.0 * m_yy_I * ux_I + 589.0 * ux_I * ux_I -
								  45.0 * m_xy_I * (6.0 * m_yy_I + 23.0 * ux_I - 23.0 * uy_I) -
								  345.0 * m_yy_I * uy_I - 1467.0 * ux_I * uy_I + 589.0 * uy_I * uy_I -
								  15.0 * m_xx_I * (18.0 * m_xy_I - 6.0 * m_yy_I - 23.0 * ux_I + 23.0 * uy_I));

		double OMEGA_factor = 2.0 * OMEGA * (4998.0 + 103.0 * OMEGA);

		double sqrt_expr = std::sqrt(3.0 * rho_I * rho_I * (squared_term + OMEGA_factor * quadratic_terms));

		rhoVar = (linear_part + sqrt_expr) / (9996.0 + 206.0 * OMEGA);

		ux_t30 = -(-3.0 * m_xx_I * rho_I + 9.0 * m_xy_I * rho_I - 3.0 * m_yy_I * rho_I -
				   20.0 * rho_I * ux_I + 3.0 * rho_I * uy_I - rhoVar) /
				 (17.0 * rhoVar);

		uy_t30 = -(3.0 * m_xx_I * rho_I - 9.0 * m_xy_I * rho_I + 3.0 * m_yy_I * rho_I +
				   3.0 * rho_I * ux_I - 20.0 * rho_I * uy_I + rhoVar) /
				 (17.0 * rhoVar);

		m_xx_t45 = -(-57.0 * m_xx_I * rho_I + 18.0 * m_xy_I * rho_I - 6.0 * m_yy_I * rho_I -
					 6.0 * rho_I * ux_I + 6.0 * rho_I * uy_I - 2.0 * rhoVar) /
				   (51.0 * rhoVar);

		m_yy_t45 = (6.0 * m_xx_I * rho_I - 18.0 * m_xy_I * rho_I + 57.0 * m_yy_I * rho_I +
					6.0 * rho_I * ux_I - 6.0 * rho_I * uy_I + 2.0 * rhoVar) /
				   (51.0 * rhoVar);

		m_xy_t90 = -(3.0 * m_xx_I * rho_I - 26.0 * m_xy_I * rho_I + 3.0 * m_yy_I * rho_I +
					 3.0 * rho_I * ux_I - 3.0 * rho_I * uy_I + rhoVar) /
				   (17.0 * rhoVar);
	}
}

#endif // BOUNDARY_FUNCTIONS_CUH
