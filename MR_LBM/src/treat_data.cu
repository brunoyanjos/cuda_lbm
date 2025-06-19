#include "treat_data.cuh"

__host__
void calculate_pressure(cylinderProperties* h_cylinder_properties, unsigned int count, unsigned int step) {
	if (step == STAT_BEGIN_TIME) {
		std::ostringstream filename_stream;
		filename_stream << PATH_FILES << "/" << ID_SIM << "/" << "theta_" << ID_SIM << ".dat";
		std::string filename = filename_stream.str();

		std::ofstream data_file(filename.c_str(), std::ios::app);

		data_file << count;
		for (int i = 0; i < count; i++) {
			cylinderProperties property = h_cylinder_properties[i];

			data_file << " " << property.theta;
		}

		data_file << std::endl;
		data_file.close();
	}

	std::ostringstream filename_stream;
	filename_stream << PATH_FILES << "/" << ID_SIM << "/" << "pressure_" << ID_SIM << ".dat";
	std::string filename = filename_stream.str();

	std::ofstream data_file(filename.c_str(), std::ios::app);

	data_file << step;

	for (int i = 0; i < count; i++) {
		cylinderProperties property = h_cylinder_properties[i];

		data_file << " " << property.ps;
	}

	data_file << std::endl;

	data_file.close();
}

__host__
void calculate_forces(cylinderProperties* h_cylinder_properties, unsigned int count, unsigned int step) {
	dfloat f_x_net = 0.0;
	dfloat f_y_net = 0.0;

	for (int i = 0; i < count; i++) {
		f_x_net += h_cylinder_properties[i].Fx;
		f_y_net += h_cylinder_properties[i].Fy;
	}

	// write to a file
	std::ostringstream filename_stream;
	filename_stream << PATH_FILES << "/" << ID_SIM << "/" << "forces_" << ID_SIM << ".dat";
	std::string filename = filename_stream.str();

	std::ofstream data_file(filename.c_str(), std::ios::app);

	data_file << step << " " << f_x_net << " " << f_y_net << std::endl;
	data_file.close();
}

__host__
void calculate_inlet_density(dfloat* h_fMom, unsigned int step) {
	dfloat rho_inlet = 0;

	for (int y = 0; y < NY; y++) {
		rho_inlet += RHO_0 + h_fMom[idxMom(0 % BLOCK_NX, y % BLOCK_NY, M_RHO_INDEX, 0 / BLOCK_NX, y / BLOCK_NY)];
	}

	rho_inlet /= NY;

	std::ostringstream filename_rho;
	filename_rho << PATH_FILES << "/" << ID_SIM << "/" << "rho_inlet_" << ID_SIM << ".dat";
	std::string filename_r = filename_rho.str();

	std::ofstream rho_file(filename_r.c_str(), std::ios::app);

	rho_file << step << " " << rho_inlet << std::endl;
	rho_file.close();

}

__global__
void domain_avg(dfloat* fMom, dfloat* ux_mean, dfloat* uy_mean, unsigned int step) {
	unsigned int x = threadIdx.x + blockDim.x * blockIdx.x;
	unsigned int y = threadIdx.y + blockDim.y * blockIdx.y;

	unsigned int global_index = x + y * NX;

	int meanCounter = step - INI_MEAN_STEP;
	dfloat invCount = 1.0 / (meanCounter + 1.0);

	dfloat ux = fMom[idxMom(x % BLOCK_NX, y % BLOCK_NY, M_UX_INDEX, x / BLOCK_NX, y / BLOCK_NY)] / F_M_I_SCALE;
	dfloat uy = fMom[idxMom(x % BLOCK_NX, y % BLOCK_NY, M_UY_INDEX, x / BLOCK_NX, y / BLOCK_NY)] / F_M_I_SCALE;

	ux_mean[global_index] = (ux_mean[global_index] * meanCounter + ux) * invCount;
	uy_mean[global_index] = (uy_mean[global_index] * meanCounter + uy) * invCount;
}