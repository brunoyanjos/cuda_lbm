#pragma once

#include "var.h"

__host__ void initializeCudaEvents(cudaEvent_t &start, cudaEvent_t &stop, cudaEvent_t &start_step, cudaEvent_t &stop_step);

[[nodiscard]] __host__ dfloat recordElapsedTime(cudaEvent_t &start_step, cudaEvent_t &stop_step, int step);