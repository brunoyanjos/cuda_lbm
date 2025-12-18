#ifndef TIME_ELAPSING_CUH
#define TIME_ELAPSING_CUH

#include "var.h"

__host__ inline void time_elapsing_count(timestep &step_end, timestep &step_start, size_t step)
{
    step_end = std::chrono::high_resolution_clock::now();
    double step_time = std::chrono::duration<double>(step_end - step_start).count();

    // Calculate MLUPS for the current step
    dfloat MLUPS = (NUMBER_LBM_NODES * MACR_SAVE / 1e6) / step_time;

    std::cout << "Elapsed time: " << step_time << " seconds" << std::endl;
    std::cout << "MLUPS: " << MLUPS << std::endl;

    // Calculate remaining time
    size_t steps_remaining = N_STEPS - step;
    double total_seconds_remaining = steps_remaining * NUMBER_LBM_NODES / 1e6 / MLUPS;

    // Convert to hours, minutes, seconds
    size_t hours = static_cast<size_t>(total_seconds_remaining) / 3600;
    size_t minutes = static_cast<size_t>((total_seconds_remaining - hours * 3600) / 60);
    size_t seconds = static_cast<size_t>(total_seconds_remaining) % 60;

    std::cout << "Estimated time left: "
              << hours << "h "
              << minutes << "min "
              << seconds << "s" << std::endl;

    step_start = std::chrono::high_resolution_clock::now();
}

#endif