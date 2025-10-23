#include "saveData.cuh"

__host__ void saveMacr_coarse(dfloat *moments, unsigned int nSteps)
{
	const std::string basePath = std::string(PATH_FILES) + "/" + ID_SIM + "/COARSE/";
	const int nprocs = 1;

	auto writeInts = [](std::ofstream &f, const std::initializer_list<int> &vals)
	{
		for (auto v : vals)
			f.write(reinterpret_cast<const char *>(&v), sizeof(int));
	};

	auto writeField = [&](std::ofstream &f, int fieldIdx, float scale = 1.0f)
	{
		for (int j = 0; j < NY_COARSE; ++j)
			for (int i = 0; i < NX_COARSE; ++i)
			{
				float val = static_cast<float>(moments[idx_mom(i, j, fieldIdx, NX_COARSE)]) / scale;
				f.write(reinterpret_cast<const char *>(&val), sizeof(float));
			}
	};

	// =======================================================================
	// MASTER FILE (.p3d)
	// =======================================================================
	{
		std::ofstream out(basePath + "master.p3d");
		if (!out)
		{
			std::cerr << "Failed to open master.p3d\n";
			return;
		}

		out << "{\n\n \"auto-detect-format\": true,\n\n \"filenames\": [\n\n";
		for (int iter = 1; iter < N_STEPS; ++iter)
			if (iter % MACR_SAVE == 0)
				out << "{ \"time\": " << iter / MACR_SAVE
					<< ", \"xyz\": \"grid.x\", \"function\": \"data" << (10000000 + iter)
					<< ".f\" },\n";
		out << "\n]\n}\n";
	}

	// =======================================================================
	// GRID FILE (grid.x)
	// =======================================================================
	if (nSteps == 0)
	{
		std::ofstream grid(basePath + "grid.x", std::ios::binary);
		if (!grid)
		{
			std::cerr << "Error opening grid file\n";
			return;
		}

		writeInts(grid, {nprocs});
		for (int m = 0; m < nprocs; ++m)
			writeInts(grid, {NX_COARSE, NY_COARSE});

		// write x, y coordinates
		for (int m = 0; m < nprocs; ++m)
		{
			for (int j = 0; j < NY_COARSE; ++j)
				for (int i = 0; i < NX_COARSE; ++i)
				{
					float val = static_cast<float>(i);
					grid.write(reinterpret_cast<const char *>(&val), sizeof(float));
				}
			for (int j = 0; j < NY_COARSE; ++j)
				for (int i = 0; i < NX_COARSE; ++i)
				{
					float val = static_cast<float>(j);
					grid.write(reinterpret_cast<const char *>(&val), sizeof(float));
				}
		}
	}

	// =======================================================================
	// DATA FILE (dataXXXXXXX.f)
	// =======================================================================
	const std::string dataFileName = basePath + "data" + std::to_string(10000000 + nSteps) + ".f";
	std::ofstream data(dataFileName, std::ios::binary);
	if (!data)
	{
		std::cerr << "Error opening data file\n";
		return;
	}

	writeInts(data, {nprocs});
	for (int l = 0; l < nprocs; ++l)
		writeInts(data, {NX_COARSE, NY_COARSE, 3}); // 3 fields: rho, ux, uy

	writeField(data, M_RHO_INDEX);
	writeField(data, M_UX_INDEX, F_M_I_SCALE);
	writeField(data, M_UY_INDEX, F_M_I_SCALE);
}

__host__ void saveMacr_fine(dfloat *moments, unsigned int nSteps)
{
	const std::string basePath = std::string(PATH_FILES) + "/" + ID_SIM + "/FINE/";
	const int nprocs = 1;

	auto writeInts = [](std::ofstream &f, const std::initializer_list<int> &vals)
	{
		for (auto v : vals)
			f.write(reinterpret_cast<const char *>(&v), sizeof(int));
	};

	auto writeField = [&](std::ofstream &f, int fieldIdx, float scale = 1.0f)
	{
		for (int j = 0; j < NY_FINE; ++j)
			for (int i = 0; i < NX_FINE; ++i)
			{
				float val = static_cast<float>(moments[idx_mom(i, j, fieldIdx, NX_FINE)]) / scale;
				f.write(reinterpret_cast<const char *>(&val), sizeof(float));
			}
	};

	// =======================================================================
	// MASTER FILE (.p3d)
	// =======================================================================
	{
		std::ofstream out(basePath + "master.p3d");
		if (!out)
		{
			std::cerr << "Failed to open master.p3d\n";
			return;
		}

		out << "{\n\n \"auto-detect-format\": true,\n\n \"filenames\": [\n\n";
		for (int iter = 1; iter < N_STEPS; ++iter)
			if (iter % MACR_SAVE == 0)
				out << "{ \"time\": " << iter / MACR_SAVE
					<< ", \"xyz\": \"grid.x\", \"function\": \"data" << (10000000 + iter)
					<< ".f\" },\n";
		out << "\n]\n}\n";
	}

	// =======================================================================
	// GRID FILE (grid.x)
	// =======================================================================
	if (nSteps == 0)
	{
		std::ofstream grid(basePath + "grid.x", std::ios::binary);
		if (!grid)
		{
			std::cerr << "Error opening grid file\n";
			return;
		}

		writeInts(grid, {nprocs});
		for (int m = 0; m < nprocs; ++m)
			writeInts(grid, {NX_FINE, NY_FINE});

		// write x, y coordinates
		for (int m = 0; m < nprocs; ++m)
		{
			for (int j = 0; j < NY_FINE; ++j)
				for (int i = 0; i < NX_FINE; ++i)
				{
					float val = static_cast<float>(i) * 0.5 + NX_COARSE - 2;
					grid.write(reinterpret_cast<const char *>(&val), sizeof(float));
				}
			for (int j = 0; j < NY_FINE; ++j)
				for (int i = 0; i < NX_FINE; ++i)
				{
					float val = static_cast<float>(j) * 0.5;
					grid.write(reinterpret_cast<const char *>(&val), sizeof(float));
				}
		}
	}

	// =======================================================================
	// DATA FILE (dataXXXXXXX.f)
	// =======================================================================
	const std::string dataFileName = basePath + "data" + std::to_string(10000000 + nSteps) + ".f";
	std::ofstream data(dataFileName, std::ios::binary);
	if (!data)
	{
		std::cerr << "Error opening data file\n";
		return;
	}

	writeInts(data, {nprocs});
	for (int l = 0; l < nprocs; ++l)
		writeInts(data, {NX_FINE, NY_FINE, 3}); // 3 fields: rho, ux, uy

	writeField(data, M_RHO_INDEX);
	writeField(data, M_UX_INDEX, F_M_I_SCALE);
	writeField(data, M_UY_INDEX, F_M_I_SCALE);
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
	strSimInfo << "                 Re: " << RE << "\n";
	strSimInfo << "          Precision: float\n";
	strSimInfo << "                 NX: " << NX << "\n";
	strSimInfo << "                 NY: " << NY << "\n";
	strSimInfo << std::scientific << std::setprecision(6);
	strSimInfo << "               Umax: " << U_MAX << "\n";
	strSimInfo << "             Macr_save: " << MACR_SAVE << "\n";
	strSimInfo << "             Nsteps: " << step << "\n";
	strSimInfo << "              MLUPS: " << MLUPS << "\n";
	strSimInfo << std::scientific << std::setprecision(0);
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