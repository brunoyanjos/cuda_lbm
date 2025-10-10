#ifndef OUTPUTS_AND_MODEL_H
#define OUTPUTS_AND_MODEL_H

#define PATH_FILES "GRID"

#ifndef ID_SIM
#define ID_SIM "000"
#endif

#define BC_PROBLEM channel_flow
#define CASE_DIRECTORY cases
#define REG_ORDER 2nd_order

// clang-format off

#define COLREC STR(colrec/REG_ORDER/collision_and_reconstruction.cuh)
#define CASE_CONSTANTS STR(BC_PROBLEM/constants.h)
#define CASE_BC STR(CASE_DIRECTORY/BC_PROBLEM/boundaries.cuh)

// clang-format on

#include CASE_CONSTANTS

#endif // !OUTPUTS_AND_MODEL_H
