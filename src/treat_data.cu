#include "treat_data.cuh"

__global__ void velocity_average(dfloat *fMom, dfloat *ux_mean, dfloat *uy_mean, unsigned int step)
{
    const std::size_t idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx >= NX)
        return;

    const std::size_t x_coord = NX / 2;
    const std::size_t y_coord = NY / 2;

    const std::size_t time_counter = (step) / MACR_SAVE;
    const dfloat inv_count = 1.0f / (1.0f + time_counter);

    const std::size_t x_thread_right = x_coord % BLOCK_NX;
    const std::size_t x_block_right = x_coord / BLOCK_NX;

    const std::size_t x_thread_left = (x_coord - 1) % BLOCK_NX;
    const std::size_t x_block_left = (x_coord - 1) / BLOCK_NX;

    const std::size_t y_thread = idx % BLOCK_NY;
    const std::size_t y_block = idx / BLOCK_NY;

    const dfloat ux_left = fMom[idxMom(x_thread_left, y_thread, M_UX_INDEX, x_block_left, y_block)] / F_M_I_SCALE;
    const dfloat ux_right = fMom[idxMom(x_thread_right, y_thread, M_UX_INDEX, x_block_right, y_block)] / F_M_I_SCALE;

    const dfloat ux = (ux_left + ux_right) * 0.5;

    ux_mean[idx] = (ux_mean[idx] * time_counter + ux) * inv_count;

    const std::size_t x_thread = idx % BLOCK_NX;
    const std::size_t x_block = idx / BLOCK_NX;

    const std::size_t y_thread_top = y_coord % BLOCK_NY;
    const std::size_t y_block_top = y_coord / BLOCK_NY;

    const std::size_t y_thread_bottom = (y_coord - 1) % BLOCK_NY;
    const std::size_t y_block_bottom = (y_coord - 1) / BLOCK_NY;

    const dfloat uy_top = fMom[idxMom(x_thread, y_thread_top, M_UY_INDEX, x_block, y_block_top)] / F_M_I_SCALE;
    const dfloat uy_bottom = fMom[idxMom(x_thread, y_thread_bottom, M_UY_INDEX, x_block, y_block_bottom)] / F_M_I_SCALE;

    const dfloat uy = (uy_top + uy_bottom) * 0.5;

    uy_mean[idx] = (uy_mean[idx] * time_counter + uy) * inv_count;
}