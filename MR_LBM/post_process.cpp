#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<sstream>
#include<vector>

struct cycleAndPeriod
{
    const int start;
    const int period;
};

// void find_peaks(const std::vector<double>& signal, int& cycle_start, int& period)
// {
//     int peak1 = -1;
//     int peak2 = -1;

//     for (int i = 1; i < signal.size() - 1; ++i) {
//         if (signal[i] > signal[i - 1] && signal[i] > signal[i + 1]) {
//             if (peak1 == -1) {
//                 peak1 = i;
//             } else if (peak2 == -1) {
//                 peak2 = i;
//                 break;  // exit loop after finding second peak
//             }
//         }
//     }

//     if (peak1 == -1 || peak2 == -1) {
//         throw std::runtime_error("Error: Could not detect two peaks.");
//     }

//     cycle_start = peak1;
//     period = peak2 - peak1;
// }

[[nodiscard]] cycleAndPeriod find_peaks(const std::vector<double>& signal)
{
    int peak1 = -1;
    int peak2 = -1;

    for (int i = 1; i < signal.size() - 1; ++i) {
        if (signal[i] > signal[i - 1] && signal[i] > signal[i + 1]) {
            if (peak1 == -1) {
                peak1 = i;
            } else if (peak2 == -1) {
                peak2 = i;
                break;  // exit loop after finding second peak
            }
        }
    }

    if (peak1 == -1 || peak2 == -1) {
        throw std::runtime_error("Error: Could not detect two peaks.");
    }

   
    return {peak1, peak2 - peak1};
}

int main()
{
    std::string PATH_FILES = "CYLINDER";
    std::string ID_SIM;

    std::cout << "Enter your simulation ID:";
    std::cin >> ID_SIM;

    std::ostringstream filename_stream;
    filename_stream << PATH_FILES << "/" << ID_SIM << "/" << "forces_" << ID_SIM << ".dat";
    std::string filename = filename_stream.str();


    // Finding number of lines in the file
    std::ifstream count_file(filename);
    if(!count_file.is_open()){
        std::cerr << "Failed to open the file" << std::endl;
        return 1;
    }

    int line_count = 0;
    std::string dummy;
    while(std::getline(count_file, dummy)){
        line_count += 1;
    }
    count_file.close();

    std::cout << "The number of time steps available: " << line_count <<std::endl;
    const int nstep = line_count;


    // Allocating the variables:
    std::vector<float> F_drag(line_count), F_lift(line_count);
    std::vector<float> Cd(line_count), Cl(line_count);
    std::vector<int> time(line_count);

    // Reading force file:
    std::ifstream data_file(filename);
    if(!data_file.is_open()){
        std::cout << "Failed to open the force file!" << std::endl;
        return 1;
    }

    int step;
    double fx,fy;

    for(int i=0; i<=line_count; i++){
        data_file >> step >> fx >> fy;
        time[i] = step;
        F_drag[i] = fx;
        F_lift[i] = fy;
    }
    data_file.close();

    // for(int i=0; i<=line_count; i++){
    //     std::cout << time[i] << "\t" << F_drag[i] << "\t" << F_lift[i] << std::endl;
    // }

    // int cycle_start=0;
    // int period = 0;
    // find_peaks(F_lift, cycle_start,period);

    const cycleAndPeriod cycle_and_period = find_peaks(F_lift);

    

    std::cout << cycle_and_period.start  << std::endl;
    std::cout << cycle_and_period.period << std::endl;
    std::cout << time[cycle_and_period.start] << std::endl;

    const int cycle_start = cycle_and_period.start;
    const int cycle_period = cycle_and_period.period;
    const int n_cycle = (nstep - cycle_start)/cycle_period;
    const int cycle_end = cycle_start + (n_cycle * cycle_period);
    
    std:: cout << "Averaging over %d cycles" << n_cycle <<std::endl;

    // Strouhal number
    const float D_cy = 32.9848440;
    const float uo = 0.1;
    const float rho_infty = 1.169;
    const float f_norm = 0.50 * rho_infty * (uo * uo) * D_cy;
    float Str = D_cy/(cycle_period * uo);

    std:: cout << Str << std::endl;


    for(int i=cycle_start; i<=cycle_end; i++){
        int mean_counter = i-cycle_start;

        Cd(i) = fd_signal(i)/f_norm
        Cl(i) = fl_signal(i)/f_norm


    }













    return 0;
}