#ifndef __SAVE_DATA_H
#define __SAVE_DATA_H

#include <string>
#include <math.h>
#include <cuda.h>
#include <curand.h>
#include <cuda_runtime.h>
#include <builtin_types.h>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream> // std::cout, std::fixed
#include <iomanip>  // std::setprecision

#include "globalFunctions.h"
#include "errorDef.h"
#include "globalStructs.h"


/*
 *   Get string with simulation information
 *   @param step: simulation's step
 *   @return string with simulation info
 */
std::string getSimInfoString(int step, dfloat MLUPS);

/*
 *   Save simulation's information
 *   @param info: simulation's informations
 */
void saveSimInfo(int step, dfloat MLUPS);

/*
 *   @brief Save array content to binary file
 *   @param strFile: filename to save
 *   @param var: float variable to save
 *   @param memSize: sizeof var

 */
void saveVarBin(
    std::string strFile,
    dfloat* var,
    size_t memSize
);


__host__ void velocity_profiles(dfloat* fMom, unsigned int step);

__host__ void kinetic_energy(dfloat *fMom, unsigned int step);

void folderSetup();



__host__ void saveMacr(dfloat* h_fMom, dfloat* rho, dfloat* ux, dfloat* uy, unsigned int nSteps);

std::string getVarFilename(const std::string varName, unsigned int step, const std::string ext);

#endif