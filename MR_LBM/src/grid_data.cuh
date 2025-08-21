#ifndef GRID_DATA_CUH
#define GRID_DATA_CUH

#include <string>
#include <fstream>
#include <sstream>
#include <iostream> // std::cout, std::fixed
#include <iomanip>  // std::setprecision

#include "var.h"
#include "globalFunctions.h"

__host__ inline void coarse_velocity_profile(dfloat *moments)
{
    // 1. Defining variables to store the path
    std::ostringstream source_path;
    std::ostringstream rho_path;
    std::ostringstream ux_path;
    std::ostringstream uy_path;
    std::ostringstream mxx_path;
    std::ostringstream mxy_path;
    std::ostringstream myy_path;

    // 2. Defining the source path for all files
    source_path << PATH_FILES << "/" << ID_SIM << "/";

    // 3. Definign velocities files name
    rho_path << source_path.str() << "coarse_rho.bin";
    ux_path << source_path.str() << "coarse_ux.bin";
    uy_path << source_path.str() << "coarse_uy.bin";
    mxx_path << source_path.str() << "coarse_mxx.bin";
    mxy_path << source_path.str() << "coarse_mxy.bin";
    myy_path << source_path.str() << "coarse_myy.bin";

    // 4. Now we open the files in binary mode
    std::ofstream rho_file(rho_path.str(), std::ios::binary);
    std::ofstream ux_file(ux_path.str(), std::ios::binary);
    std::ofstream uy_file(uy_path.str(), std::ios::binary);
    std::ofstream mxx_file(mxx_path.str(), std::ios::binary);
    std::ofstream mxy_file(mxy_path.str(), std::ios::binary);
    std::ofstream myy_file(myy_path.str(), std::ios::binary);

    const size_t x_coord_0 = NX_COARSE / 2 - 1;
    const size_t x_coord_1 = NX_COARSE / 2;

    for (size_t y = 0; y < NY_COARSE + N_OVERLAP_LAYER; ++y)
    {
        const dfloat rho_0 = moments[coarse_moment_idx(x_coord_0, y, M_RHO_INDEX)];
        const dfloat ux_0 = moments[coarse_moment_idx(x_coord_0, y, M_UX_INDEX)] / F_M_I_SCALE;
        const dfloat uy_0 = moments[coarse_moment_idx(x_coord_0, y, M_UY_INDEX)] / F_M_I_SCALE;
        const dfloat mxx_0 = moments[coarse_moment_idx(x_coord_0, y, M_MXX_INDEX)] / F_M_II_SCALE;
        const dfloat mxy_0 = moments[coarse_moment_idx(x_coord_0, y, M_MXY_INDEX)] / F_M_IJ_SCALE;
        const dfloat myy_0 = moments[coarse_moment_idx(x_coord_0, y, M_MYY_INDEX)] / F_M_II_SCALE;

        const dfloat rho_1 = moments[coarse_moment_idx(x_coord_1, y, M_RHO_INDEX)];
        const dfloat ux_1 = moments[coarse_moment_idx(x_coord_1, y, M_UX_INDEX)] / F_M_I_SCALE;
        const dfloat uy_1 = moments[coarse_moment_idx(x_coord_1, y, M_UY_INDEX)] / F_M_I_SCALE;
        const dfloat mxx_1 = moments[coarse_moment_idx(x_coord_1, y, M_MXX_INDEX)] / F_M_II_SCALE;
        const dfloat mxy_1 = moments[coarse_moment_idx(x_coord_1, y, M_MXY_INDEX)] / F_M_IJ_SCALE;
        const dfloat myy_1 = moments[coarse_moment_idx(x_coord_1, y, M_MYY_INDEX)] / F_M_II_SCALE;

        const dfloat rho = (rho_0 + rho_1) * static_cast<dfloat>(0.5);
        const dfloat ux = (ux_0 + ux_1) * static_cast<dfloat>(0.5);
        const dfloat uy = (uy_0 + uy_1) * static_cast<dfloat>(0.5);
        const dfloat mxx = (mxx_0 + mxx_1) * static_cast<dfloat>(0.5);
        const dfloat mxy = (mxy_0 + mxy_1) * static_cast<dfloat>(0.5);
        const dfloat myy = (myy_0 + myy_1) * static_cast<dfloat>(0.5);

        rho_file.write(reinterpret_cast<const char *>(&rho), sizeof(dfloat));
        ux_file.write(reinterpret_cast<const char *>(&ux), sizeof(dfloat));
        uy_file.write(reinterpret_cast<const char *>(&uy), sizeof(dfloat));
        mxx_file.write(reinterpret_cast<const char *>(&mxx), sizeof(dfloat));
        mxy_file.write(reinterpret_cast<const char *>(&mxy), sizeof(dfloat));
        myy_file.write(reinterpret_cast<const char *>(&myy), sizeof(dfloat));
    }
}

__host__ inline void fine_velocity_profile(dfloat *moments)
{
    // 1. Defining variables to store the path
    std::ostringstream source_path;
    std::ostringstream rho_path;
    std::ostringstream ux_path;
    std::ostringstream uy_path;
    std::ostringstream mxx_path;
    std::ostringstream mxy_path;
    std::ostringstream myy_path;

    // 2. Defining the source path for all files
    source_path << PATH_FILES << "/" << ID_SIM << "/";

    // 3. Definign velocities files name
    rho_path << source_path.str() << "fine_rho.bin";
    ux_path << source_path.str() << "fine_ux.bin";
    uy_path << source_path.str() << "fine_uy.bin";
    mxx_path << source_path.str() << "fine_mxx.bin";
    mxy_path << source_path.str() << "fine_mxy.bin";
    myy_path << source_path.str() << "fine_myy.bin";

    // 4. Now we open the files in binary mode
    std::ofstream rho_file(rho_path.str(), std::ios::binary);
    std::ofstream ux_file(ux_path.str(), std::ios::binary);
    std::ofstream uy_file(uy_path.str(), std::ios::binary);
    std::ofstream mxx_file(mxx_path.str(), std::ios::binary);
    std::ofstream mxy_file(mxy_path.str(), std::ios::binary);
    std::ofstream myy_file(myy_path.str(), std::ios::binary);

    const size_t x_coord_0 = NX_FINE / 2 - 1;
    const size_t x_coord_1 = NX_FINE / 2;

    for (size_t y = 0; y < NY_FINE; ++y)
    {
        const dfloat rho_0 = moments[fine_moment_idx(x_coord_0, y, M_RHO_INDEX)];
        const dfloat ux_0 = moments[fine_moment_idx(x_coord_0, y, M_UX_INDEX)] / F_M_I_SCALE;
        const dfloat uy_0 = moments[fine_moment_idx(x_coord_0, y, M_UY_INDEX)] / F_M_I_SCALE;
        const dfloat mxx_0 = moments[fine_moment_idx(x_coord_0, y, M_MXX_INDEX)] / F_M_II_SCALE;
        const dfloat mxy_0 = moments[fine_moment_idx(x_coord_0, y, M_MXY_INDEX)] / F_M_IJ_SCALE;
        const dfloat myy_0 = moments[fine_moment_idx(x_coord_0, y, M_MYY_INDEX)] / F_M_II_SCALE;

        const dfloat rho_1 = moments[fine_moment_idx(x_coord_1, y, M_RHO_INDEX)];
        const dfloat ux_1 = moments[fine_moment_idx(x_coord_1, y, M_UX_INDEX)] / F_M_I_SCALE;
        const dfloat uy_1 = moments[fine_moment_idx(x_coord_1, y, M_UY_INDEX)] / F_M_I_SCALE;
        const dfloat mxx_1 = moments[fine_moment_idx(x_coord_1, y, M_MXX_INDEX)] / F_M_II_SCALE;
        const dfloat mxy_1 = moments[fine_moment_idx(x_coord_1, y, M_MXY_INDEX)] / F_M_IJ_SCALE;
        const dfloat myy_1 = moments[fine_moment_idx(x_coord_1, y, M_MYY_INDEX)] / F_M_II_SCALE;

        const dfloat rho = (rho_0 + rho_1) * static_cast<dfloat>(0.5);
        const dfloat ux = (ux_0 + ux_1) * static_cast<dfloat>(0.5);
        const dfloat uy = (uy_0 + uy_1) * static_cast<dfloat>(0.5);
        const dfloat mxx = (mxx_0 + mxx_1) * static_cast<dfloat>(0.5);
        const dfloat mxy = (mxy_0 + mxy_1) * static_cast<dfloat>(0.5);
        const dfloat myy = (myy_0 + myy_1) * static_cast<dfloat>(0.5);

        rho_file.write(reinterpret_cast<const char *>(&rho), sizeof(dfloat));
        ux_file.write(reinterpret_cast<const char *>(&ux), sizeof(dfloat));
        uy_file.write(reinterpret_cast<const char *>(&uy), sizeof(dfloat));
        mxx_file.write(reinterpret_cast<const char *>(&mxx), sizeof(dfloat));
        mxy_file.write(reinterpret_cast<const char *>(&mxy), sizeof(dfloat));
        myy_file.write(reinterpret_cast<const char *>(&myy), sizeof(dfloat));
    }
}

#endif