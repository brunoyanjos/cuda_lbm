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