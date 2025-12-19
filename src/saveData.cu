#include "saveData.cuh"

__host__ void saveMacr(
	dfloat *h_fMom, dfloat *rho, dfloat *ux, dfloat *uy, unsigned int nSteps)
{
	int x, y;
	size_t indexMacr;

	// ---------------------------------------------------------------------
	// 1) Reconstrói campos macroscópicos (pode zerar se quiser testar)
	// ---------------------------------------------------------------------
	for (y = 0; y < NY; y++)
	{
		for (x = 0; x < NX; x++)
		{
			indexMacr = idxScalarGlobal(x, y);

			rho[indexMacr] =
				RHO_0 +
				h_fMom[idxMom(x % BLOCK_NX, y % BLOCK_NY,
							  M_RHO_INDEX, x / BLOCK_NX, y / BLOCK_NY)];

			ux[indexMacr] =
				h_fMom[idxMom(x % BLOCK_NX, y % BLOCK_NY,
							  M_UX_INDEX, x / BLOCK_NX, y / BLOCK_NY)] /
				F_M_I_SCALE;

			uy[indexMacr] =
				h_fMom[idxMom(x % BLOCK_NX, y % BLOCK_NY,
							  M_UY_INDEX, x / BLOCK_NX, y / BLOCK_NY)] /
				F_M_I_SCALE;
		}
	}

	const int nprocs = 1;
	const int nz = 1; // 2D em PLOT3D
	const int nf = 3; // rho, ux, uy
	const std::string prefix = "data";
	const std::string suffix = ".f";

	// Caminho base dos arquivos
	std::string basePath = std::string(PATH_FILES) + "/" + ID_SIM + "/";

	// ---------------------------------------------------------------------
	// 2) master.p3d – só precisa ser criado uma vez (nSteps == 0)
	//    aqui desligamos o auto-detect e descrevemos o formato explicitamente
	// ---------------------------------------------------------------------
	if (nSteps == 0)
	{
		std::string metaFile = basePath + "master.p3d";
		std::ofstream out(metaFile.c_str());
		if (!out)
		{
			std::cerr << "Erro abrindo master.p3d\n";
		}
		else
		{
			out << "{\n";
			out << "  \"auto-detect-format\" : false,\n";
			out << "  \"format\"            : \"binary\",\n";
			out << "  \"byte-order\"        : \"little\",\n";
			out << "  \"precision\"         : 32,\n";
			out << "  \"multi-grid\"        : true,\n";
			out << "  \"language\"          : \"C\",\n";
			out << "  \"blanking\"          : false,\n";
			out << "  \"2D\"                : true,\n";
			out << "  \"filenames\"         : [\n";

			bool first = true;
			for (int iter = 1; iter < N_STEPS; ++iter)
			{
				if (iter % MACR_SAVE == 0)
				{
					if (!first)
						out << ",\n";
					first = false;

					out << "    { \"time\" : "
						<< (iter / MACR_SAVE)
						<< ", \"xyz\" : \"grid.x\", \"function\" : \""
						<< prefix << (10000000 + iter) << suffix << "\" }";
				}
			}
			out << "\n  ]\n";
			out << "}\n";
		}
	}

	// ---------------------------------------------------------------------
	// 3) grid.x – geometria PLOT3D (XYZ)
	//    também só precisa uma vez, no início da simulação
	// ---------------------------------------------------------------------
	if (nSteps == 0)
	{
		std::string gridFile = basePath + "grid.x";
		std::ofstream grid(gridFile.c_str(), std::ios::binary);
		if (!grid)
		{
			std::cerr << "Erro abrindo grid.x\n";
			return;
		}

		// número de blocos
		grid.write(reinterpret_cast<const char *>(&nprocs), sizeof(int));

		// (NI, NJ, NK) de cada bloco
		for (int p = 0; p < nprocs; ++p)
		{
			grid.write(reinterpret_cast<const char *>(&NX), sizeof(int));
			grid.write(reinterpret_cast<const char *>(&NY), sizeof(int));
			grid.write(reinterpret_cast<const char *>(&nz), sizeof(int));
		}

		// X(i,j,k)
		for (int k = 0; k < nz; ++k)
			for (int j = 0; j < NY; ++j)
				for (int i = 0; i < NX; ++i)
				{
					float v = static_cast<float>(i);
					grid.write(reinterpret_cast<const char *>(&v), sizeof(float));
				}

		// Y(i,j,k)
		for (int k = 0; k < nz; ++k)
			for (int j = 0; j < NY; ++j)
				for (int i = 0; i < NX; ++i)
				{
					float v = static_cast<float>(j);
					grid.write(reinterpret_cast<const char *>(&v), sizeof(float));
				}

		// Z(i,j,k) = 0
		for (int k = 0; k < nz; ++k)
			for (int j = 0; j < NY; ++j)
				for (int i = 0; i < NX; ++i)
				{
					float v = 0.0f;
					grid.write(reinterpret_cast<const char *>(&v), sizeof(float));
				}

		grid.close();
	}

	// ---------------------------------------------------------------------
	// 4) Arquivo de dados PLOT3D (functions) – um por nSteps salvo
	// ---------------------------------------------------------------------
	{
		std::ostringstream fname;
		fname << basePath << prefix << (10000000 + nSteps) << suffix;
		std::string filename = fname.str();

		std::ofstream datafile(filename.c_str(), std::ios::binary);
		if (!datafile)
		{
			std::cerr << "Erro abrindo " << filename << "\n";
			return;
		}

		// nprocs
		datafile.write(reinterpret_cast<const char *>(&nprocs), sizeof(int));

		// (NI, NJ, NK, NF) por bloco
		for (int p = 0; p < nprocs; ++p)
		{
			datafile.write(reinterpret_cast<const char *>(&NX), sizeof(int));
			datafile.write(reinterpret_cast<const char *>(&NY), sizeof(int));
			datafile.write(reinterpret_cast<const char *>(&nz), sizeof(int));
			datafile.write(reinterpret_cast<const char *>(&nf), sizeof(int));
		}

		// Campos na ordem: rho, ux, uy (3 * NX*NY*nz floats)
		for (int k = 0; k < nz; ++k)
		{
			// rho
			for (int j = 0; j < NY; ++j)
				for (int i = 0; i < NX; ++i)
				{
					indexMacr = idxScalarGlobal(i, j);
					float val = static_cast<float>(rho[indexMacr]);
					datafile.write(reinterpret_cast<const char *>(&val), sizeof(float));
				}

			// ux
			for (int j = 0; j < NY; ++j)
				for (int i = 0; i < NX; ++i)
				{
					indexMacr = idxScalarGlobal(i, j);
					float val = static_cast<float>(ux[indexMacr]);
					datafile.write(reinterpret_cast<const char *>(&val), sizeof(float));
				}

			// uy
			for (int j = 0; j < NY; ++j)
				for (int i = 0; i < NX; ++i)
				{
					indexMacr = idxScalarGlobal(i, j);
					float val = static_cast<float>(uy[indexMacr]);
					datafile.write(reinterpret_cast<const char *>(&val), sizeof(float));
				}
		}

		datafile.close();
	}

	// Saídas auxiliares da tua infra antiga (binários separados)
	std::string strFileRho = getVarFilename("rho", nSteps, ".bin");
	std::string strFileUx = getVarFilename("ux", nSteps, ".bin");
	std::string strFileUy = getVarFilename("uy", nSteps, ".bin");
}

std::string getVarFilename(
	const std::string varName,
	unsigned int step,
	const std::string ext)
{
	unsigned int n_zeros = 0, pot_10 = 10;
	unsigned int aux1 = 1000000; // 6 numbers on step
	// calculate number of zeros
	if (step != 0)
		for (n_zeros = 0; step * pot_10 < aux1; pot_10 *= 10)
			n_zeros++;
	else
		n_zeros = 5;

	// generates the file name as "PATH_FILES/id/id_varName000000.bin"
	std::string strFile = PATH_FILES;
	strFile += "/";
	strFile += ID_SIM;
	strFile += "/";
	strFile += ID_SIM;
	strFile += "_";
	strFile += varName;
	for (unsigned int i = 0; i < n_zeros; i++)
		strFile += "0";
	strFile += std::to_string(step);
	strFile += ext;

	return strFile;
}

void saveVarBin(
	std::string strFile,
	dfloat *var,
	size_t memSize)
{
	FILE *outFile = nullptr;

	outFile = fopen(strFile.c_str(), "wb");

	if (outFile != nullptr)
	{
		fwrite(var, memSize, 1, outFile);
		fclose(outFile);
	}
	else
	{
		printf("Error saving \"%s\" \nProbably wrong path!\n", strFile.c_str());
	}
}

std::string getSimInfoString(int step, dfloat MLUPS)
{
#define BOOL_TO_STRING(b) ((b) ? "true" : "false")
	std::ostringstream strSimInfo("");

	strSimInfo << std::scientific;
	strSimInfo << std::setprecision(6);

	strSimInfo << "---------------------------- SIMULATION INFORMATION ----------------------------\n";
	strSimInfo << "      Simulation ID: " << ID_SIM << "\n";
	strSimInfo << "       Velocity set: D2Q9\n";
	strSimInfo << "          Reg Order: " << STR(REG_ORDER) << "\n";
	strSimInfo << "                 Re: " << RE << "\n";
	strSimInfo << "          Precision: float\n";
	strSimInfo << "                 NX: " << NX << "\n";
	strSimInfo << "                 NY: " << NY << "\n";
	strSimInfo << std::scientific << std::setprecision(6);
	/*strSimInfo << "                Tau: " << TAU << "\n";*/
	strSimInfo << "               Umax: " << U_MAX << "\n";
	strSimInfo << "             Macr_save: " << MACR_SAVE << "\n";
	strSimInfo << "             Nsteps: " << step << "\n";
	strSimInfo << "              MLUPS: " << MLUPS << "\n";
	strSimInfo << std::scientific << std::setprecision(0);
	strSimInfo << "                 BX: " << BLOCK_NX << "\n";
	strSimInfo << "                 BY: " << BLOCK_NY << "\n";
	strSimInfo << "--------------------------------------------------------------------------------\n";

	return strSimInfo.str();
}

void folderSetup()
{
	// Windows
#if defined(_WIN32)
	std::string strPath;
	strPath = PATH_FILES;
	strPath += "\\\\"; // adds "\\"
	strPath += ID_SIM;
	std::string cmd = "md ";
	cmd += strPath;
	system(cmd.c_str());
	return;
#endif // !_WIN32

	// Unix
#if defined(__APPLE__) || defined(__MACH__) || defined(__linux__)
	std::string strPath;
	strPath = PATH_FILES;
	strPath += "/";
	strPath += ID_SIM;
	std::string cmd = "mkdir -p ";
	cmd += strPath;
	const int i = system(cmd.c_str());
	static_cast<void>(i);
	return;
#endif // !Unix
	printf("I don't know how to setup folders for your operational system :(\n");
	return;
}

void saveSimInfo(int step, dfloat MLUPS)
{
	std::string strInf = PATH_FILES;
	strInf += "/";
	strInf += ID_SIM;
	strInf += "/";
	strInf += "info.txt"; // generate file name (with path)
	FILE *outFile = nullptr;

	outFile = fopen(strInf.c_str(), "w");
	if (outFile != nullptr)
	{
		std::string strSimInfo = getSimInfoString(step, MLUPS);
		fprintf(outFile, "%s\n", strSimInfo.c_str());
		fclose(outFile);
	}
	else
	{
		printf("Error saving \"%s\" \nProbably wrong path!\n", strInf.c_str());
	}
}
/**/

void create_vtk(const LBMState &state, const unsigned int &step)
{
	std::ostringstream filename;
	filename << PATH_FILES << "/" << ID_SIM << "/vtk/" << "output_" << std::setw(6) << std::setfill('0') << step << ".vtk";

	std::ofstream file(filename.str());
	if (!file.is_open())
	{
		std::cerr << "Could not open VTK file for writing: "
				  << filename.str() << std::endl;
		return;
	}

	// --- HEADER ---
	file << "# vtk DataFile Version 3.0\n";
	file << "LBM output\n";
	file << "ASCII\n";
	file << "DATASET STRUCTURED_POINTS\n";
	file << "DIMENSIONS " << NX << " " << NY << " 1\n";
	file << "ORIGIN 0 0 0\n";
	file << "SPACING 1 1 1\n";
	file << "POINT_DATA " << NX * NY << "\n";

	// --- DENSITY ---
	file << "SCALARS rho float 1\n";
	file << "LOOKUP_TABLE default\n";

	for (int y = 0; y < NY; y++)
	{
		for (int x = 0; x < NX; x++)
		{
			int idx = idxScalarBlock(x % BLOCK_NX, y % BLOCK_NY, x / BLOCK_NX, y / BLOCK_NY);
			float rho = state.h_rho[idx] + RHO_0;
			file << rho << "\n";
		}
	}

	// --- VELOCITY UX ---
	file << "SCALARS ux float 1\n";
	file << "LOOKUP_TABLE default\n";

	for (int y = 0; y < NY; y++)
	{
		for (int x = 0; x < NX; x++)
		{
			int idx = idxScalarBlock(x % BLOCK_NX, y % BLOCK_NY, x / BLOCK_NX, y / BLOCK_NY);
			file << state.h_ux[idx] / F_M_I_SCALE << "\n";
		}
	}

	// --- VELOCITY UY ---
	file << "SCALARS uy float 1\n";
	file << "LOOKUP_TABLE default\n";

	for (int y = 0; y < NY; y++)
	{
		for (int x = 0; x < NX; x++)
		{
			int idx = idxScalarBlock(x % BLOCK_NX, y % BLOCK_NY, x / BLOCK_NX, y / BLOCK_NY);
			file << state.h_uy[idx] / F_M_I_SCALE << "\n";
		}
	}

	// --- VELOCITY ---
	file << "VECTORS velocity float\n";
	for (int y = 0; y < NY; y++)
	{
		for (int x = 0; x < NX; x++)
		{
			int idx = idxScalarBlock(x % BLOCK_NX, y % BLOCK_NY, x / BLOCK_NX, y / BLOCK_NY);
			file << state.h_ux[idx] / F_M_I_SCALE << " "
				 << state.h_uy[idx] / F_M_I_SCALE << " "
				 << 0.0f << "\n";
		}
	}

	file.close();
}