#include <iostream>
#include <iomanip>

#include "post_stats.h"
#include "file_utils.h"

int main()
{
    std::string PATH_FILES = "./../CYLINDER";
    std::string ID_SIM;

    std::cout << "Enter your simulation ID:";
    std::cin >> ID_SIM;

    // -----------------Reading input file---------------------------------------------------
    std::string filename = PATH_FILES + "/" + ID_SIM + "/" + "input_stat_" + ID_SIM + ".dat";
    InputParameters params = read_input(filename);
    const double D_cy = params.Dcy;
    const double uo = params.uo;
    const double rho_infty = params.rho_infty;
    const double p_infty = rho_infty / 3.0;
    const double f_norm = 0.50 * rho_infty * (uo * uo) * D_cy;
    const double p_norm = 0.50 * rho_infty * (uo * uo);
    const double t_norm = D_cy / uo;

    // -------------------Getting number of time stamps in a file------------------------------
    filename = PATH_FILES + "/" + ID_SIM + "/" + "forces_" + ID_SIM + ".dat";
    const int nsteps = count_lines_in_file(filename);
    std::cout << "The number of time steps available: " << nsteps << std::endl;

    // ----------------------------------------------------------------------------------------------
    // Reading theta file
    filename = PATH_FILES + "/" + ID_SIM + "/" + "theta_" + ID_SIM + ".dat";
    const int npoints = get_npoints_on_cylinder(filename);

    std::cout << "nsteps = " << nsteps << "\n";
    std::cout << "npoints = " << npoints << "\n";
    std::cout << "ncy = " << npoints << "\n";
    std::cout << "Dcy = " << params.Dcy << "\n";
    std::cout << "uo = " << params.uo << "\n";
    std::cout << "rho_infty = " << params.rho_infty << "\n";

    CFDdata data = allocate_vectors(nsteps, npoints);

    filename = PATH_FILES + "/" + ID_SIM + "/" + "theta_" + ID_SIM + ".dat";
    int status = read_theta(filename, npoints, data.theta);

    filename = PATH_FILES + "/" + ID_SIM + "/" + "forces_" + ID_SIM + ".dat";
    status = read_forces(filename, nsteps, data.time, data.F_drag, data.F_lift);

    filename = PATH_FILES + "/" + ID_SIM + "/" + "pressure_" + ID_SIM + ".dat";
    status = read_pressure(filename, nsteps, npoints, data.ps);
    // status = read_pressure_binary(filename, nsteps, npoints, data.ps);

    std::cout << data.ps[nsteps-1][npoints-1] << std::endl;

    // Calculate the lift and drag coefficients
    for (int i = 0; i < nsteps; i++)
    {
        data.ntime[i] = data.time[i] / t_norm;
        data.Cd[i] = data.F_drag[i] / f_norm;
        data.Cl[i] = data.F_lift[i] / f_norm;
    }
    filename = PATH_FILES + "/" + ID_SIM + "/" + "coeff_" + ID_SIM + ".dat";
    std::ofstream coeff_file(filename);
    if (!coeff_file.is_open())
    {
        std::cout << "Error creating coeff_file!!" << std::endl;
        return 1;
    }
    for (int i = 0; i < nsteps; i++)
    {
        coeff_file << std::setprecision(12) << std::fixed << data.ntime[i] << " " << data.Cd[i] << " " << data.Cl[i] << std::endl;
    }
    coeff_file.close();

    // gaussian_smoothing(data.ntime,data.Cl, data.Cl_smooth,0.1); // sigma value 0.01
    moving_average_smoothing(data.ntime, data.Cl, data.Cl_smooth, 21); // odd window numer 21
    filename = PATH_FILES + "/" + ID_SIM + "/" + "coeff_smooth_" + ID_SIM + ".dat";
    std::ofstream coeff_file2(filename);
    if (!coeff_file2.is_open())
    {
        std::cout << "Error creating cl_smooth_file!!" << std::endl;
        return 1;
    }
    for (int i = 0; i < nsteps; i++)
    {
        coeff_file2 << std::setprecision(12) << std::fixed << data.ntime[i] << " " << data.Cl[i] << " " << data.Cl_smooth[i] << std::endl;
    }
    coeff_file2.close();

    const cycleAndPeriod cycle_and_period = find_peaks(data.Cl_smooth);

    const int cycle_start = cycle_and_period.start;
    const int cycle_period = cycle_and_period.period;
    const int n_cycle = (nsteps - cycle_start) / cycle_period;
    const int cycle_end = cycle_start + (n_cycle * cycle_period);

    std::cout << "start: " << cycle_start << std::endl;
    std::cout << "period: " << cycle_period << std::endl;
    std::cout << "start time: " << data.time[cycle_start] << std::endl;
    std::cout << "end time (one cycle): " << data.time[cycle_start + cycle_period] << std::endl;
    std::cout << "number of cycles: " << n_cycle << std::endl;

    // Strouhal number
    double Str = D_cy / (cycle_period * uo);

    for (int i = cycle_start; i <= cycle_end; i++)
    {
        int mean_counter = i - cycle_start;

        data.Cd_avg = ((double)mean_counter * data.Cd_avg + data.Cd[i]) / ((double)(mean_counter + 1));
        data.Cl_avg = ((double)mean_counter * data.Cl_avg + (data.Cl[i] * data.Cl[i])) / ((double)(mean_counter + 1));

        for (int j = 0; j < npoints; j++)
        {
            data.ps_avg[j] = ((double)mean_counter * data.ps_avg[j] + data.ps[i][j]) / ((double)(mean_counter + 1));
            data.Cp[j] = (data.ps_avg[j] - p_infty) / p_norm;
        }
    }
    data.Cl_avg = sqrt(data.Cl_avg);

    std::cout << "---------------------------------------------------------" << std::endl;
    std::cout << "Averaging over\t" << n_cycle << "\tcycles" << std::endl;
    std::cout << "---------------------------------------------------------" << std::endl;
    std::cout << "Strouhal Number:" << Str << std::endl;
    std::cout << "Drag Coefficient:" << data.Cd_avg << std::endl;
    std::cout << "Lift Coefficient:" << data.Cl_avg << std::endl;
    std::cout << "---------------------------------------------------------" << std::endl;

    filename = PATH_FILES + "/" + ID_SIM + "/" + "cl_time_" + ID_SIM + ".dat";
    std::string filename2 = PATH_FILES + "/" + ID_SIM + "/" + "cd_time_" + ID_SIM + ".dat";
    std::ofstream data_file(filename);
    std::ofstream data_file2(filename2);
    if (!data_file.is_open())
    {
        std::cout << "Error creating cl_file!!" << std::endl;
        return 1;
    }
    if (!data_file2.is_open())
    {
        std::cout << "Error creating cd_file!!" << std::endl;
        return 1;
    }

    for (int i = cycle_start; i <= cycle_end; i++)
    {
        double tstar = (data.time[i] - data.time[cycle_start]) / (double)cycle_period;
        data_file << tstar << " " << data.Cl[i] << std::endl;
        data_file2 << tstar << " " << data.Cd[i] << std::endl;
    }
    data_file.close();
    data_file2.close();

    for (int i = 0; i < npoints; i++)
    {
        data.idx[i] = i;
    }

    for (int i = 0; i < npoints - 1; i++)
    {
        for (int j = i + 1; j < npoints; j++)
        {
            int index_i = data.idx[i];
            int index_j = data.idx[j];
            if (data.theta[index_i] > data.theta[index_j])
            {
                int temp = index_i;
                data.idx[i] = index_j;
                data.idx[j] = temp;
            }
        }
    }

    for (int i = 0; i < npoints; i++)
    {
        int index = data.idx[i];
        data.theta_sort[i] = data.theta[index];
        data.Cp_sort[i] = data.Cp[index];
    }

    // gaussian_smoothing(data.theta_sort,data.Cp_sort, data.Cp_smooth,0.1);
    moving_average_smoothing(data.theta_sort, data.Cp_sort, data.Cp_smooth, 11);

    filename = PATH_FILES + "/" + ID_SIM + "/" + "cp_" + ID_SIM + ".dat";
    std::ofstream cp_file(filename);
    if (!cp_file.is_open())
    {
        std::cout << "Error creating Cp_file!!" << std::endl;
        return 1;
    }

    for (int i = 0; i < npoints; i++)
    {
        if (data.theta_sort[i] < PI)
        {
            double theta = (data.theta_sort[i] * 180.0 / PI) - 180.0;
            cp_file << std::abs(theta) << " " << data.Cp_sort[i] << " " << data.Cp_smooth[i] << std::endl;
        }
    }

    return 0;
}