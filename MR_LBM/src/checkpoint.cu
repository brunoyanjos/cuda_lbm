#include "checkpoint.cuh"

// Saves simulation state to a binary file
__host__ void save_checkpoint(int current_step, dfloat *moments)
{
    std::ostringstream file_path;

    file_path << PATH_FILES << "/" << ID_SIM << "/checkpoint.bin";

    std::ofstream file(file_path.str(), std::ios::binary);

    if (!file)
    {
        std::cerr << "Error opening checkpoint file" << " - Reason: " << std::strerror(errno) << std::endl;
        return;
    }

    // 1. Write current simulation step
    file.write(reinterpret_cast<char *>(&current_step), sizeof(int));

    // 3. Write moments
    file.write(reinterpret_cast<char *>(moments), NUMBER_LBM_NODES * NUMBER_MOMENTS * sizeof(dfloat));

    // Cleanup
    file.close();
}

__host__ bool load_checkpoint(int *current_step, dfloat *&moments)
{
    std::ostringstream file_path;
    file_path << PATH_FILES << "/" << ID_SIM << "/checkpoint.bin";
    std::ifstream file(file_path.str(), std::ios::binary);

    if (!file)
    {
        std::cerr << "Error opening checkpoint file - Reason: "
                  << std::strerror(errno) << std::endl;
        return false;
    }

    // 1. Read current simulation step
    file.read(reinterpret_cast<char *>(current_step), sizeof(int));

    // 2. Create host buffer
    const size_t total_elements = NUMBER_LBM_NODES * NUMBER_MOMENTS;

    // 3. Read data directly into host buffer
    file.read(reinterpret_cast<char *>(moments),
              total_elements * sizeof(dfloat));

    // 4. Check if read was successful
    if (!file)
    {
        std::cerr << "Error reading checkpoint data - Reason: "
                  << std::strerror(errno) << std::endl;
        return false;
    }

    return true;
}

__host__ void truncate_tke_file(unsigned int last_valid_step)
{
    std::ostringstream tke_path;
    tke_path << PATH_FILES << "/" << ID_SIM << "/total_kinetic_energy.bin";
    const std::string filename = tke_path.str();

    // 1. Read all valid records
    std::ifstream in_file(filename, std::ios::binary);
    std::vector<std::tuple<dfloat, dfloat>> records;

    dfloat t_star, tke_val;
    dfloat valid_t_star =  last_valid_step * U_MAX / NX;
    
    while (in_file.read(reinterpret_cast<char *>(&t_star), sizeof(dfloat)) &&
           in_file.read(reinterpret_cast<char *>(&tke_val), sizeof(dfloat)))
    {
        if (t_star <= valid_t_star)
        {
            records.emplace_back(t_star, tke_val);
        }
    }

    in_file.close();

    // 2. Rewrite file with only valid records
    std::ofstream out_file(filename, std::ios::binary | std::ios::trunc);
    for (const auto &record : records)
    {
        t_star = std::get<0>(record);
        tke_val = std::get<1>(record);

        out_file.write(reinterpret_cast<const char *>(&t_star), sizeof(dfloat));
        out_file.write(reinterpret_cast<const char *>(&tke_val), sizeof(dfloat));
    }
}