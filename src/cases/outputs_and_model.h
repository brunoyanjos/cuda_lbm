#ifndef OUTPUTS_AND_MODEL_H
#define OUTPUTS_AND_MODEL_H

#define PATH_FILES "ANNUL"

#ifndef ID_SIM
#define ID_SIM "001"
#endif

#define BC_PROBLEM ldc
#define CASE_DIRECTORY cases
#define REG_ORDER 2nd_order

// clang-format off

#define COLREC STR(colrec/REG_ORDER/collision_and_reconstruction.cuh)
#define CASE_CONSTANTS STR(constants.h)
#define CASE_BC STR(CASE_DIRECTORY/boundaries.cuh)

// clang-format on

#include CASE_CONSTANTS

constexpr size_t CHECKPOINT_STEP = 1000000;
constexpr bool LOAD_CHECKPOINT = false;

#endif // !OUTPUTS_AND_MODEL_H
