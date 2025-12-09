#include "saveData.cuh"
#include "globalFunctions.h"
#include <iostream>

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

	strSimInfo << std::fixed;
	strSimInfo << std::setprecision(6);

	strSimInfo << "---------------------------- SIMULATION INFORMATION ----------------------------\n";
	strSimInfo << "      Simulation ID =  " << ID_SIM << "\n";
	strSimInfo << "       Velocity set = D2Q9\n";
	strSimInfo << "                 Re = " << RE << "\n";
	strSimInfo << "          Precision = float\n";
	strSimInfo << "                 NX = " << NX << "\n";
	strSimInfo << "                 NY = " << NY << "\n";
	strSimInfo << std::fixed << std::setprecision(6);
	strSimInfo << "                 uo = " << U_MAX << "\n";
	strSimInfo << "          Macr_save = " << MACR_SAVE << "\n";
	strSimInfo << "             Nsteps = " << step << "\n";
	strSimInfo << "              MLUPS = " << MLUPS << "\n";
	strSimInfo << std::fixed << std::setprecision(0);
	strSimInfo << "                 BX = " << BLOCK_NX << "\n";
	strSimInfo << "                 BY = " << BLOCK_NY << "\n";
	strSimInfo << "                  D = " << D << "\n";
	strSimInfo << std::fixed << std::setprecision(5);
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
	strPath += "/";
	strPath += "plots";
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
	strInf += ID_SIM;
	strInf += "_info.txt"; // generate file name (with path)
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