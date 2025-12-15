#ifndef __GLOBAL_STRUCTS_H
#define __GLOBAL_STRUCTS_H

#include "var.h"
#include "errorDef.h"

typedef struct ghostData
{
    dfloat *X_0;
    dfloat *X_1;
    dfloat *Y_0;
    dfloat *Y_1;
} GhostData;

typedef struct ghostInterfaceData
{
    ghostData fGhost;
    ghostData gGhost;
    ghostData h_fGhost;
} GhostInterfaceData;

struct LBMState
{
    // DEVICE
    uint8_t *d_node_type;

    dfloat *d_rho;

    dfloat *d_ux;
    dfloat *d_uy;

    dfloat *d_mxx;
    dfloat *d_mxy;
    dfloat *d_myy;

    // HOST
    uint8_t *h_node_type;
    dfloat *h_rho;

    dfloat *h_ux;
    dfloat *h_uy;

    dfloat *h_mxx;
    dfloat *h_mxy;
    dfloat *h_myy;

    // SIZES
    size_t bytes_fields;
    size_t bytes_types;
};

#endif //__GLOBAL_STRUCTS_H
