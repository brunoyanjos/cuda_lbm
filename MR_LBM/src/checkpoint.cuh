#ifndef CHECKPOINT_CUH
#define CHECKPOINT

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

#include "var.h"

__host__ void save_checkpoint(int current_step, dfloat *moments);

__host__ bool load_checkpoint(int *current_step, dfloat *&moments);

#endif