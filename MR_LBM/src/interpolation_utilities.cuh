#ifndef INTERPOLATION_UTILITIES_CUH
#define INTERPOLATION_UTILITIES_CUH

#include "var.h"

__host__ inline dfloat first_order_interpolation(const dfloat g_minus_h, const dfloat g_plus_h)
{
    const dfloat inv_two = static_cast<dfloat>(1) / static_cast<dfloat>(2);

    return (g_minus_h + g_plus_h) * inv_two;
}

__host__ inline dfloat second_order_interpolation(const dfloat g_minus_h, const dfloat g_plus_h, const dfloat d_plus_3h)
{

    const dfloat three_by_eight = static_cast<dfloat>(3) / static_cast<dfloat>(8);
    const dfloat three_by_four = static_cast<dfloat>(3) / static_cast<dfloat>(4);
    const dfloat one_by_eight = static_cast<dfloat>(1) / static_cast<dfloat>(8);

    return three_by_eight * g_minus_h + three_by_four * g_plus_h - one_by_eight * d_plus_3h;
}

__host__ inline dfloat third_order_interpolation(const dfloat g_minus_3h, const dfloat g_minus_h, const dfloat g_plus_h, const dfloat d_plus_3h)
{
    const dfloat gph_plus_gmh = g_minus_h + g_plus_h;
    const dfloat gp3h_plus_gm3h = g_minus_3h + d_plus_3h;

    const dfloat inv_sixteen = static_cast<dfloat>(1) / static_cast<dfloat>(16);

    return static_cast<dfloat>(9) * inv_sixteen * gph_plus_gmh - inv_sixteen * gp3h_plus_gm3h;
}

#endif