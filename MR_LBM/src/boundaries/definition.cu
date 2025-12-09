#include "definition.cuh"
#include "node_type.h"
#include "../var.h"

namespace boundary
{
    __device__ uint8_t definition(int x, int y)
    {
        const dfloat xc_local = xc;
        const dfloat yc_local = yc;

        // ---------------------------
        // Distance from the cell to the center
        // ---------------------------
        const dfloat dx = dfloat(x) - xc_local;
        const dfloat dy = dfloat(y) - yc_local;
        const dfloat dist = dsqrt(dx * dx + dy * dy);

        // ---------------------------
        // Radii definitions
        // ---------------------------
        const dfloat inner_radius = dfloat(D) / 2.0;
        const dfloat outer_radius = dfloat(NX - 1) / 2.0;

        // Smoothing thickness (LBM-friendly transition band)
        const dfloat smooth = 0.5;

        // ---------------------------
        // INNER BOUNDARY (smooth band)
        // ---------------------------
        const dfloat d_inner = dabs(dist - inner_radius);

        // Inside the inner solid region OR within the smoothing band → solid
        if (dist <= inner_radius || d_inner <= smooth)
        {
            return SOLID_NODE;
        }

        // ---------------------------
        // OUTER BOUNDARY (smooth band)
        // ---------------------------
        const dfloat d_outer = dabs(dist - outer_radius);

        // Outside the outer radius → solid region
        // Within the smoothing band around the outer boundary → solid
        if (d_outer <= smooth || dist >= outer_radius)
        {
            return SOLID_NODE;
        }

        // ---------------------------
        // FLUID REGION (between inner and outer radii)
        // ---------------------------
        return BULK;
    }
}