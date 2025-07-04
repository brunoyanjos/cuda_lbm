#ifndef OUTPUTS_AND_MODEL_H
#define OUTPUTS_AND_MODEL_H

#define PATH_FILES "JET_FLOW"
#define ID_SIM "001"

#define BC_PROBLEM jet_flow
#define CASE_DIRECTORY cases
#define REG_ORDER 2nd_order

#define COLREC STR(colrec/REG_ORDER/collision_and_reconstruction.cuh)
#define CASE_CONSTANTS STR(BC_PROBLEM/constants.h)
#define CASE_BC STR(CASE_DIRECTORY/BC_PROBLEM/boundaries.cuh)

#include CASE_CONSTANTS

#endif // !OUTPUTS_AND_MODEL_H
