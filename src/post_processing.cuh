#pragma once

#include "globalStructs.h"
#include "nodeTypeMap.h"
#include "globalFunctions.h"
#include <fstream>
#include <iostream>
#include <vector>

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

__host__ inline void write_velocity_profile(LBMState state)
{
    std::string path_uy_theta_0 = std::string(PATH_FILES) + "/" + ID_SIM + "/uy_theta_0.bin";
    std::string path_ux_theta_90 = std::string(PATH_FILES) + "/" + ID_SIM + "/ux_theta_90.bin";
    std::string path_uy_theta_180 = std::string(PATH_FILES) + "/" + ID_SIM + "/uy_theta_180.bin";
    std::string path_ux_theta_270 = std::string(PATH_FILES) + "/" + ID_SIM + "/ux_theta_270.bin";

    std::ofstream file_theta_0(path_uy_theta_0, std::ios::binary | std::ios::app | std::ios::out);
    std::ofstream file_theta_90(path_ux_theta_90, std::ios::binary | std::ios::app | std::ios::out);
    std::ofstream file_theta_180(path_uy_theta_180, std::ios::binary | std::ios::app | std::ios::out);
    std::ofstream file_theta_270(path_ux_theta_270, std::ios::binary | std::ios::app | std::ios::out);

    std::vector<dfloat> uy_theta_0;
    std::vector<dfloat> ux_theta_90;
    std::vector<dfloat> uy_theta_180;
    std::vector<dfloat> ux_theta_270;

    int y0 = NY / 2 - 1;
    int y1 = NY / 2;

    for (int x = 0; x < NX; ++x)
    {
        int idx0 = idxBlockCoord(x, y0);
        int idx1 = idxBlockCoord(x, y1);

        uint8_t node_type = state.h_node_type[idx0];

        if (node_type == SOLID_NODE)
            continue;

        const dfloat uy0 = state.h_uy[idx0] / F_M_I_SCALE;
        const dfloat uy1 = state.h_uy[idx1] / F_M_I_SCALE;

        const dfloat ux0 = state.h_ux[idx0] / F_M_I_SCALE;
        const dfloat ux1 = state.h_ux[idx1] / F_M_I_SCALE;

        const dfloat dx = (x - xc);
        const dfloat dy0 = (y0 - yc);
        const dfloat dy1 = (y1 - yc);

        const dfloat dx2 = dx * dx;
        const dfloat dy02 = dy0 * dy0;
        const dfloat dy12 = dy1 * dy1;

        const dfloat radii0 = dsqrt(dx2 + dy02);
        const dfloat radii1 = dsqrt(dx2 + dy12);

        const dfloat u0 = dx / radii0 * uy0 - dy0 / radii0 * ux0;
        const dfloat u1 = dx / radii1 * uy1 - dy1 / radii1 * ux1;

        dfloat u = (u0 + u1) * 0.5f;

        if (x < NX / 2)
        {
            uy_theta_0.push_back(u);
        }
        else
        {
            uy_theta_180.push_back(u);
        }
    }

    int x0 = NX / 2 - 1;
    int x1 = NX / 2;

    for (int y = 0; y < NY; ++y)
    {
        int idx0 = idxBlockCoord(x0, y);
        int idx1 = idxBlockCoord(x1, y);

        uint8_t node_type = state.h_node_type[idx0];

        if (node_type == SOLID_NODE)
            continue;

        dfloat ux0 = state.h_ux[idx0] / F_M_I_SCALE;
        dfloat ux1 = state.h_ux[idx1] / F_M_I_SCALE;

        dfloat uy0 = state.h_uy[idx0] / F_M_I_SCALE;
        dfloat uy1 = state.h_uy[idx1] / F_M_I_SCALE;

        const dfloat dx0 = (x0 - xc);
        const dfloat dx1 = (x1 - xc);
        const dfloat dy = (y - yc);

        const dfloat dx02 = dx0 * dx0;
        const dfloat dx12 = dx1 * dx1;
        const dfloat dy2 = dy * dy;

        const dfloat radii0 = dsqrt(dx02 + dy2);
        const dfloat radii1 = dsqrt(dx12 + dy2);

        const dfloat u0 = dx0 / radii0 * uy0 - dy / radii0 * ux0;
        const dfloat u1 = dx1 / radii1 * uy1 - dy / radii1 * ux1;

        dfloat u = (u0 + u1) * 0.5f;

        if (y < NY / 2)
        {
            ux_theta_270.push_back(u);
        }
        else
        {
            ux_theta_90.push_back(u);
        }
    }

    file_theta_0.write(reinterpret_cast<char *>(uy_theta_0.data()), uy_theta_0.size() * sizeof(dfloat));
    file_theta_90.write(reinterpret_cast<char *>(ux_theta_90.data()), ux_theta_90.size() * sizeof(dfloat));
    file_theta_180.write(reinterpret_cast<char *>(uy_theta_180.data()), uy_theta_180.size() * sizeof(dfloat));
    file_theta_270.write(reinterpret_cast<char *>(ux_theta_270.data()), ux_theta_270.size() * sizeof(dfloat));

    file_theta_0.close();
    file_theta_90.close();
    file_theta_180.close();
    file_theta_270.close();
}

__host__ inline void write_pressure_profile(LBMState state)
{
    std::string path_pressure_theta_0 = std::string(PATH_FILES) + "/" + ID_SIM + "/uy_theta_0.bin";
    std::string path_pressure_theta_90 = std::string(PATH_FILES) + "/" + ID_SIM + "/ux_theta_90.bin";
    std::string path_pressure_theta_180 = std::string(PATH_FILES) + "/" + ID_SIM + "/uy_theta_180.bin";
    std::string path_pressure_theta_270 = std::string(PATH_FILES) + "/" + ID_SIM + "/ux_theta_270.bin";

    std::ofstream file_theta_0(path_pressure_theta_0, std::ios::binary | std::ios::app | std::ios::out);
    std::ofstream file_theta_90(path_pressure_theta_90, std::ios::binary | std::ios::app | std::ios::out);
    std::ofstream file_theta_180(path_pressure_theta_180, std::ios::binary | std::ios::app | std::ios::out);
    std::ofstream file_theta_270(path_pressure_theta_270, std::ios::binary | std::ios::app | std::ios::out);

    std::vector<dfloat> pressure_theta_0;
    std::vector<dfloat> pressure_theta_90;
    std::vector<dfloat> pressure_theta_180;
    std::vector<dfloat> pressure_theta_270;

    int y0 = NY / 2 - 1;
    int y1 = NY / 2;

    for (int x = 0; x < NX; ++x)
    {
        int idx0 = idxBlockCoord(x, y0);
        int idx1 = idxBlockCoord(x, y1);

        uint8_t node_type = state.h_node_type[idx0];

        if (node_type == SOLID_NODE)
            continue;

        const dfloat uy0 = state.h_uy[idx0] / F_M_I_SCALE;
        const dfloat uy1 = state.h_uy[idx1] / F_M_I_SCALE;

        const dfloat ux0 = state.h_ux[idx0] / F_M_I_SCALE;
        const dfloat ux1 = state.h_ux[idx1] / F_M_I_SCALE;

        const dfloat dx = (x - xc);
        const dfloat dy0 = (y0 - yc);
        const dfloat dy1 = (y1 - yc);

        const dfloat dx2 = dx * dx;
        const dfloat dy02 = dy0 * dy0;
        const dfloat dy12 = dy1 * dy1;

        const dfloat radii0 = dsqrt(dx2 + dy02);
        const dfloat radii1 = dsqrt(dx2 + dy12);

        const dfloat u0 = dx / radii0 * uy0 - dy0 / radii0 * ux0;
        const dfloat u1 = dx / radii1 * uy1 - dy1 / radii1 * ux1;

        dfloat u = (u0 + u1) * 0.5f;

        if (x < NX / 2)
        {
            pressure_theta_0.push_back(u);
        }
        else
        {
            pressure_theta_180.push_back(u);
        }
    }

    int x0 = NX / 2 - 1;
    int x1 = NX / 2;

    for (int y = 0; y < NY; ++y)
    {
        int idx0 = idxBlockCoord(x0, y);
        int idx1 = idxBlockCoord(x1, y);

        uint8_t node_type = state.h_node_type[idx0];

        if (node_type == SOLID_NODE)
            continue;

        dfloat ux0 = state.h_ux[idx0] / F_M_I_SCALE;
        dfloat ux1 = state.h_ux[idx1] / F_M_I_SCALE;

        dfloat uy0 = state.h_uy[idx0] / F_M_I_SCALE;
        dfloat uy1 = state.h_uy[idx1] / F_M_I_SCALE;

        const dfloat dx0 = (x0 - xc);
        const dfloat dx1 = (x1 - xc);
        const dfloat dy = (y - yc);

        const dfloat dx02 = dx0 * dx0;
        const dfloat dx12 = dx1 * dx1;
        const dfloat dy2 = dy * dy;

        const dfloat radii0 = dsqrt(dx02 + dy2);
        const dfloat radii1 = dsqrt(dx12 + dy2);

        const dfloat u0 = dx0 / radii0 * uy0 - dy / radii0 * ux0;
        const dfloat u1 = dx1 / radii1 * uy1 - dy / radii1 * ux1;

        dfloat u = (u0 + u1) * 0.5f;

        if (y < NY / 2)
        {
            pressure_theta_270.push_back(u);
        }
        else
        {
            pressure_theta_90.push_back(u);
        }
    }

    file_theta_0.write(reinterpret_cast<char *>(pressure_theta_0.data()), pressure_theta_0.size() * sizeof(dfloat));
    file_theta_90.write(reinterpret_cast<char *>(pressure_theta_90.data()), pressure_theta_90.size() * sizeof(dfloat));
    file_theta_180.write(reinterpret_cast<char *>(pressure_theta_180.data()), pressure_theta_180.size() * sizeof(dfloat));
    file_theta_270.write(reinterpret_cast<char *>(pressure_theta_270.data()), pressure_theta_270.size() * sizeof(dfloat));

    file_theta_0.close();
    file_theta_90.close();
    file_theta_180.close();
    file_theta_270.close();
}