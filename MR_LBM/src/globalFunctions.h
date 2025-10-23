
#ifndef __GLOBAL_FUNCTIONS_H
#define __GLOBAL_FUNCTIONS_H

#include <builtin_types.h> // for device variables
#include "var.h"

__host__
    size_t __forceinline__
    idx_grid(size_t x, size_t y, size_t nx)
{
    return x + nx * y;
}

__host__
    size_t __forceinline__
    idx_mom(size_t x, size_t y, size_t mom_idx, size_t nx)
{
    return (x + nx * y) * NUMBER_MOMENTS + mom_idx;
}

__host__
    size_t __forceinline__
    idx_pop(size_t x, size_t y, size_t pop_idx, size_t nx)
{
    return (x + nx * y) * Q + pop_idx;
}

#endif // !__GLOBAL_FUNCTIONS_H
