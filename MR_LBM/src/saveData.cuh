#ifndef __SAVE_DATA_H
#define __SAVE_DATA_H

#include <string>
#include <math.h>
#include <cuda.h>
#include <curand.h>
#include <cuda_runtime.h>
#include <builtin_types.h>

#include <vector>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream> // std::cout, std::fixed
#include <iomanip>  // std::setprecision

#include "errorDef.h"
#include "globalStructs.h"
#include "globalFunctions.h"

std::string getSimInfoString(int step, dfloat MLUPS);

void saveSimInfo(int step, dfloat MLUPS);

void saveVarBin(
    std::string strFile,
    dfloat *var,
    size_t memSize);

void folderSetup();

std::string getVarFilename(const std::string varName, unsigned int step, const std::string ext);

__host__ inline void write_average_velocity_profile(dfloat *ux)
{
    std::string ux_file_path = std::string(PATH_FILES) + "/" + ID_SIM + "/ux_dy_average.bin";

    std::ofstream ux_file(ux_file_path, std::ios::binary | std::ios::out);

    // buffer temporário para armazenar em ordem global
    std::vector<dfloat> ux_vec(NY);

    int x0 = NX / 2 - 1;
    int x1 = NX / 2;

    // reconstruir indices globais (x,y,z)
    for (int y = 0; y < NY; y++)
    {
        int tx0 = x0 % BLOCK_NX;
        int tx1 = x1 % BLOCK_NX;

        int bx0 = x0 / BLOCK_NX;
        int bx1 = x1 / BLOCK_NX;

        int ty = y % BLOCK_NY;
        int by = y / BLOCK_NY;

        int bid0 = idxScalarBlock(tx0, ty, bx0, by);
        int bid1 = idxScalarBlock(tx0, ty, bx0, by);
        int bid2 = idxScalarBlock(tx1, ty, bx1, by);
        int bid3 = idxScalarBlock(tx1, ty, bx1, by);

        dfloat ux0 = ux[bid0];
        dfloat ux1 = ux[bid1];
        dfloat ux2 = ux[bid2];
        dfloat ux3 = ux[bid3];

        ux_vec[y] = (ux0 + ux1 + ux2 + ux3) * static_cast<dfloat>(0.25);
    }

    ux_file.write(reinterpret_cast<char *>(ux_vec.data()), NY * sizeof(dfloat));

    ux_file.close();
}

#endif