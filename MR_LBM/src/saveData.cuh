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

#include "errorDef.h"
#include "globalStructs.h"

std::string getSimInfoString(int step, dfloat MLUPS);

void saveSimInfo(int step, dfloat MLUPS);

void saveVarBin(
    std::string strFile,
    dfloat *var,
    size_t memSize);

void folderSetup();

std::string getVarFilename(const std::string varName, unsigned int step, const std::string ext);

#endif