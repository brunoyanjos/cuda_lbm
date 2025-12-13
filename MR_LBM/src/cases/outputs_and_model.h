#ifndef OUTPUTS_AND_MODEL_H
#define OUTPUTS_AND_MODEL_H

#define PATH_FILES "ANNUL"

#ifndef ID_SIM
#define ID_SIM "000"
#endif

#define BC_PROBLEM cylinder
#define CASE_DIRECTORY cases
#define REG_ORDER 2nd_order

#include "constants.h"

constexpr bool CALCULATE_PRESSURE = true;
constexpr bool CALCULATE_FORCES = true;
constexpr bool CALCULATE_RHO = true;
constexpr bool SAVE_DOMAIN = true;

constexpr unsigned int STAT_BEGIN_TIME = (tstar * D / U_MAX);
constexpr unsigned int STAT_END_TIME = N_STEPS;

constexpr int INI_MEAN_STEP = 0;

#endif // !OUTPUTS_AND_MODEL_H
