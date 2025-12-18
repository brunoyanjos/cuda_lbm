#pragma once

#include "globalStructs.h"
#include "nodeTypeMap.h"
#include "globalFunctions.h"
#include <fstream>
#include <iostream>

__host__ inline void write_tke(LBMState state)
{
    std::string path = std::string(PATH_FILES) + "/" + ID_SIM + "/kinetic_energy.bin";

    std::ofstream file(path, std::ios::binary | std::ios::app | std::ios::out);

    dfloat tke = 0;
    int total_points = 0;

    for (int y = 0; y < NY; ++y)
    {
        for (int x = 0; x < NX; ++x)
        {
            int idx = idxBlockCoord(x, y);

            uint8_t node_type = state.h_node_type[idx];

            if (node_type == SOLID_NODE)
                continue;

            dfloat ux = state.h_ux[idx] / F_M_I_SCALE;
            dfloat uy = state.h_uy[idx] / F_M_I_SCALE;

            dfloat ux2 = ux * ux;
            dfloat uy2 = uy * uy;

            dfloat u2 = ux2 + uy2;

            tke += dfloat(0.5) * u2;
            total_points++;
        }
    }

    const dfloat denom = U_MAX * U_MAX * total_points;

    tke /= denom;

    file.write(reinterpret_cast<char *>(&tke), sizeof(dfloat));
    file.close();
}