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

#endif //__GLOBAL_STRUCTS_H
