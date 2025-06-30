#include "saveData.cuh"

__host__ void velocity_profiles(dfloat *fMom, unsigned int step)
{
    // 1. Defining variables to store the path
    std::ostringstream source_path;
    std::ostringstream ux_path;
    std::ostringstream uy_path;

    // 2. Defining the source path for all files
    source_path << PATH_FILES << "/" << ID_SIM << "/";

    // 3. Definign velocities files name
    ux_path << source_path.str() << "velocity_x_" << std::setw(8) << std::setfill('0') << step << ".bin";
    uy_path << source_path.str() << "velocity_y_" << std::setw(8) << std::setfill('0') << step << ".bin";

    // 4. Now we open the files in binary mode
    std::ofstream ux_file(ux_path.str(), std::ios::binary);
    std::ofstream uy_file(uy_path.str(), std::ios::binary);

    // 5. Make sure the files were create properly
    // 5.1. First we make sure that ux_file is fine
    if (!ux_file)
    {
        std::cerr << "Error opening ux_file" << " - Reason: " << std::strerror(errno) << std::endl;
        return;
    }

    // 5.2. Then we do the same to the uy_file
    if (!uy_file)
    {
        std::cerr << "Error opening uy_file" << " - Reason: " << std::strerror(errno) << std::endl;
        return;
    }

    // 6. Now we are able to star the loop to calculate the velocities on he center line,
    // and the write then on the file

    const std::size_t x_coord = NX / 2;
    const std::size_t y_coord = NY / 2;

    for (std::size_t y = 0; y < NY; ++y)
    {
        const std::size_t x_thread_right = x_coord % BLOCK_NX;
        const std::size_t x_block_right = x_coord / BLOCK_NX;

        const std::size_t x_thread_left = (x_coord - 1) % BLOCK_NX;
        const std::size_t x_block_left = (x_coord - 1) / BLOCK_NX;

        const std::size_t y_thread = y % BLOCK_NY;
        const std::size_t y_block = y / BLOCK_NY;

        const dfloat ux_left = fMom[idxMom(x_thread_left, y_thread, M_UX_INDEX, x_block_left, y_block)] / F_M_I_SCALE;
        const dfloat ux_right = fMom[idxMom(x_thread_right, y_thread, M_UX_INDEX, x_block_right, y_block)] / F_M_I_SCALE;

        const dfloat ux = (ux_left + ux_right) * 0.5;

        ux_file.write(reinterpret_cast<const char *>(&ux), sizeof(dfloat));
    }

    for (std::size_t x = 0; x < NX; ++x)
    {
        const std::size_t x_thread = x % BLOCK_NX;
        const std::size_t x_block = x / BLOCK_NX;

        const std::size_t y_thread_top = y_coord % BLOCK_NY;
        const std::size_t y_block_top = y_coord / BLOCK_NY;

        const std::size_t y_thread_bottom = (y_coord - 1) % BLOCK_NY;
        const std::size_t y_block_bottom = (y_coord - 1) / BLOCK_NY;

        const dfloat uy_top = fMom[idxMom(x_thread, y_thread_top, M_UY_INDEX, x_block, y_block_top)] / F_M_I_SCALE;
        const dfloat uy_bottom = fMom[idxMom(x_thread, y_thread_bottom, M_UY_INDEX, x_block, y_block_bottom)] / F_M_I_SCALE;

        const dfloat uy = (uy_top + uy_bottom) * 0.5;

        uy_file.write(reinterpret_cast<const char *>(&uy), sizeof(dfloat));
    }
}

__host__ void kinetic_energy(dfloat *fMom, unsigned int step)
{
    // 1. Defining a variable that will store the totl kinetic energy (TKE) path
    std::ostringstream tke_path;

    // 2. Defining the path to TKE
    tke_path << PATH_FILES << "/" << ID_SIM << "/" << "total_kinetic_energy.bin";

    // 3. Open the file as a binary
    std::ofstream tke_file(tke_path.str(), std::ios::binary | std::ios::app);

    // 4. Now we make sure that the file open properly
    if (!tke_file)
    {
        std::cerr << "Error opening tke_file" << " - Reason: " << std::strerror(errno) << std::endl;
        return;
    }

    // 5. Finally we will sum over all the domain the kinetic energy, to have a total kinetic energy

    dfloat tke_sum = 0.0;

    for (std::size_t y = 0; y < NY; ++y)
    {
        for (std::size_t x = 0; x < NX; ++x)
        {
            const std::size_t x_thread = x % BLOCK_NX;
            const std::size_t y_thread = y % BLOCK_NY;

            const std::size_t x_block = x / BLOCK_NX;
            const std::size_t y_block = y / BLOCK_NY;

            const dfloat ux = fMom[idxMom(x_thread, y_thread, M_UX_INDEX, x_block, y_block)] / F_M_I_SCALE;
            const dfloat uy = fMom[idxMom(x_thread, y_thread, M_UY_INDEX, x_block, y_block)] / F_M_I_SCALE;

            const dfloat ux2 = ux * ux;
            const dfloat uy2 = uy * uy;

            const dfloat u2 = (ux2 + uy2) * 0.5f;

            tke_sum += u2;
        }
    }

    tke_sum /= (U_MAX * U_MAX * NX * NY);
    const dfloat t_star = step * U_MAX / NX;

    tke_file.write(reinterpret_cast<const char*> (&t_star), sizeof(dfloat));
    tke_file.write(reinterpret_cast<const char*> (&tke_sum), sizeof(dfloat));
}