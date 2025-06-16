#ifndef AUX_FUNCTIONS_CUH
#define AUX_FUNCTIONS_CUH

#include "../../globalStructs.h"

__device__ inline cylinderProperties findCylindeProperty(
    cylinderProperties *cylinderProperties, size_t counter,
    size_t x, size_t y)
{
    for(size_t i = 0; i < counter; i++) {
        size_t xb = cylinderProperties[i].xb;
        size_t yb = cylinderProperties[i].yb;
    
        if(x == xb && y == yb) {
            return cylinderProperties[i];
        }
    }
}

__device__
inline void immersedBoundaryLoop(
	const int(&incomings)[9],
	const dfloat(&pop)[9],
	dfloat* rhoVar,
	dfloat* m_xx_t45,
	dfloat* m_yy_t45,
	dfloat* m_xy_t90,
	int x,
	int y
)
{
	dfloat rho_I = 0;
	dfloat m_xx_I = 0;
	dfloat m_yy_I = 0;
	dfloat m_xy_I = 0;

	const dfloat radius = sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc));

	const dfloat cos_theta = (x - xc) / radius;
	const dfloat sen_theta = (y - yc) / radius;
	const dfloat sen_two_theta = 2.0 * cos_theta * sen_theta;
	const dfloat cos_two_theta = cos_theta * cos_theta - sen_theta * sen_theta;

	for (std::size_t i = 0; i < 9; i++)
	{
		if (incomings[i] == 1)
		{
			if (ROTATIONAL_COORDINATES)
			{
				const dfloat Hxx = (cx[i] * cx[i]) - cs2;
				const dfloat Hyy = (cy[i] * cy[i]) - cs2;
				const dfloat Hxy = cx[i] * cy[i];

				const dfloat Hxx_p = Hxx * cos_theta * cos_theta + Hyy * sen_theta * sen_theta + Hxy * sen_two_theta;
				const dfloat Hyy_p = Hxx * sen_theta * sen_theta + Hyy * cos_theta * cos_theta - Hxy * sen_two_theta;
				const dfloat Hxy_p = (Hyy - Hxx) * 0.5 * sen_two_theta + Hxy * cos_two_theta;

				rho_I += pop[i];
				m_xx_I += pop[i] * Hxx_p;
				m_yy_I += pop[i] * Hyy_p;
				m_xy_I += pop[i] * Hxy_p;
			}
			else
			{
				const dfloat Hxx = (cx[i] * cx[i]) - cs2;
				const dfloat Hyy = (cy[i] * cy[i]) - cs2;
				const dfloat Hxy = cx[i] * cy[i];
                
				rho_I += pop[i];
				m_xx_I += pop[i] * Hxx;
				m_yy_I += pop[i] * Hyy;
				m_xy_I += pop[i] * Hxy;
			}
		}
	}

	const dfloat inv_rho_I = 1.0 / rho_I;

	*rhoVar = rho_I;
	*m_xx_t45 = m_xx_I * inv_rho_I;
	*m_yy_t45 = m_yy_I * inv_rho_I;
	*m_xy_t90 = m_xy_I * inv_rho_I;
}

#endif