#include "post_processing.cuh"
#include "globalFunctions.h"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>

void create_vtk(const LBMState &state, const unsigned int &step)
{
    std::ostringstream filename;
    filename << "output_" << std::setw(6) << std::setfill('0') << step << ".vtk";

    std::ofstream file(filename.str());
    if (!file.is_open())
    {
        std::cerr << "Could not open VTK file for writing: "
                  << filename.str() << std::endl;
        return;
    }

    // --- HEADER ---
    file << "# vtk DataFile Version 3.0\n";
    file << "LBM output\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_POINTS\n";
    file << "DIMENSIONS " << NX << " " << NY << " 1\n";
    file << "ORIGIN 0 0 0\n";
    file << "SPACING 1 1 1\n";
    file << "POINT_DATA " << NX * NY << "\n";

    // --- DENSITY ---
    file << "SCALARS rho float 1\n";
    file << "LOOKUP_TABLE default\n";

    for (int j = 0; j < NY; j++)
    {
        for (int i = 0; i < NX; i++)
        {
            int idx = idxScalarBlock(i % BLOCK_NX, j % BLOCK_NY, i / BLOCK_NX, j / BLOCK_NY);
            float rho = state.h_rho[idx] + RHO_0;
            file << rho << "\n";
        }
    }

    // --- VELOCITY ---
    file << "VECTORS velocity float\n";
    for (int j = 0; j < NY; j++)
    {
        for (int i = 0; i < NX; i++)
        {
            int idx = idxScalarBlock(i % BLOCK_NX, j % BLOCK_NY, i / BLOCK_NX, j / BLOCK_NY);
            file << state.h_ux[idx] << " "
                 << state.h_uy[idx] << " "
                 << 0.0f << "\n";
        }
    }

    file.close();
}