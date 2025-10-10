#ifndef __VAR_H
#define __VAR_H

#include <builtin_types.h> // for devices variables
#include <stdint.h>		   // for uint32_t
#include <map>

#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include <limits>
#include <cstring>

/* ----------------------------- PROBLEM DEFINE ---------------------------- */
typedef float dfloat;

#define GPU_INDEX 0
/* --------------------------  SIMULATION DEFINES -------------------------- */

#define STR_IMPL(A) #A
#define STR(A) STR_IMPL(A)

#include "cases/outputs_and_model.h"
#include "definitions.h"
#endif //__VAR_H