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

void folderSetup();

__host__ void saveMacr_coarse(
    dfloat *moments, unsigned int nSteps);

__host__ void saveMacr_fine(
    dfloat *moments, unsigned int nSteps);

std::string getVarFilename(const std::string varName, unsigned int step, const std::string ext);

#endif