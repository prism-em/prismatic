// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// Prismatic is distributed under the GNU General Public License (GPL)
// If you use Prismatic, we kindly ask that you cite the following papers:

// 1. Ophus, C.: A fast image simulation algorithm for scanning
//    transmission electron microscopy. Advanced Structural and
//    Chemical Imaging 3(1), 13 (2017)

// 2. Pryor, Jr., A., Ophus, C., and Miao, J.: A Streaming Multi-GPU
//    Implementation of Image Simulation Algorithms for Scanning
//	  Transmission Electron Microscopy. arXiv:1706.08563 (2017)

#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include "PRISM_entry.h"
#include "configure.h"
#include "ArrayND.h"
#include "PRISM01_calcPotential.h"
#include "PRISM02_calcSMatrix.h"
#include "PRISM03_calcOutput.h"
#include "params.h"
#include "H5Cpp.h"
#include "utility.h"
#include "fileIO.h"

namespace Prismatic
{
using namespace std;
Parameters<PRISMATIC_FLOAT_PRECISION> PRISM_entry(Metadata<PRISMATIC_FLOAT_PRECISION> &meta)
{
	Parameters<PRISMATIC_FLOAT_PRECISION> prismatic_pars;
	try
	{ // read atomic coordinates
		prismatic_pars = Parameters<PRISMATIC_FLOAT_PRECISION>(meta);
	}
	catch (...)
	{
		std::cout << "Terminating" << std::endl;
		exit(1);
	}
	prismatic_pars.meta.toString();

	//        to_xyz(prismatic_pars.atoms, "/Users/ajpryor/Documents/MATLAB/multislice/PRISM/build/test.XYZ", "comment", 5.43,5.43,5.43);

	prismatic_pars.outputFile = H5::H5File(prismatic_pars.meta.filenameOutput.c_str(), H5F_ACC_TRUNC);
	setupOutputFile(prismatic_pars);
	// compute projected potentials
	prismatic_pars.fpFlag = 0;

	if(prismatic_pars.meta.importSMatrix)
	{
		std::cout << "Skipping PRISM01. Using precalculated scattering matrix from: "  << prismatic_pars.meta.importFile << std::endl;
	}
	else if(prismatic_pars.meta.importPotential)
	{
		std::cout << "Using precalculated potential from " << prismatic_pars.meta.importFile << std::endl;
		std::string groupName = "4DSTEM_simulation/data/realslices/";
		std::string baseName = "ppotential_fp";
		H5::H5File inFile = H5::H5File(prismatic_pars.meta.importFile.c_str(), H5F_ACC_RDONLY);
		H5::Group realslices = inFile.openGroup(groupName.c_str());
		int configurations = countDataGroups(realslices, baseName);

		std::cout << configurations << " frozen phonon configurations available in " << prismatic_pars.meta.importFile << std::endl;
		int tmp_fp = prismatic_pars.meta.numFP;

		//if user requests more than the available configurations, only run number available
		if(prismatic_pars.meta.numFP > 1)
		{
			prismatic_pars.meta.numFP = (tmp_fp > configurations) ? configurations : tmp_fp;
			std::cout << "User requested " << tmp_fp  << " frozen phonons." << std::endl;
			std::cout << "Running " << prismatic_pars.meta.numFP << " frozen phonons out of " << configurations << " available configurations." << std::endl;
		}
		else
		{
			if( not prismatic_pars.meta.userSpecifiedNumFP)
			{
				//if user specifically specifies to run a single frozen phonon, this is skipped and only the first configuration will run
				prismatic_pars.meta.numFP = configurations;
			}
		}
	}

	// compute compact S-matrix
	if(prismatic_pars.meta.importSMatrix)
	{
		std::string groupName = "4DSTEM_simulation/data/realslices/";
		std::string baseName = "smatrix_fp";
		H5::H5File inFile = H5::H5File(prismatic_pars.meta.importFile.c_str(), H5F_ACC_RDONLY);
		H5::Group realslices = inFile.openGroup(groupName.c_str());
		int configurations = countDataGroups(realslices, baseName);

		std::cout << configurations << " frozen phonon configurations available in " << prismatic_pars.meta.importFile << std::endl;
		int tmp_fp = prismatic_pars.meta.numFP;

		//if user requests more than the available configurations, only run number available
		if(prismatic_pars.meta.numFP > 1)
		{
			prismatic_pars.meta.numFP = (tmp_fp > configurations) ? configurations : tmp_fp;
			std::cout << "User requested " << tmp_fp  << " frozen phonons." << std::endl;
			std::cout << "Running " << prismatic_pars.meta.numFP << " frozen phonons out of " << configurations << " available configurations." << std::endl;
		}
		else
		{
			if( not prismatic_pars.meta.userSpecifiedNumFP)
			{
				//if user specifically specifies to run a single frozen phonon, this is skipped and only the first configuration will run
				prismatic_pars.meta.numFP = configurations;
			}
		}
	}

	for(auto i = 0; i < prismatic_pars.meta.numFP; i++)
	{
		PRISM_runFP(prismatic_pars, i);
	}

	//average data by fp
	for (auto &i : prismatic_pars.net_output)
		i /= prismatic_pars.meta.numFP;
	for (auto &i : prismatic_pars.net_output_c)
		i /= prismatic_pars.meta.numFP;

	if (prismatic_pars.meta.saveDPC_CoM)
	{
		for (auto &j : prismatic_pars.net_DPC_CoM)
			j /= prismatic_pars.meta.numFP; //since squared intensities are used to calculate DPC_CoM, this is incoherent averaging
	}

	saveSTEM(prismatic_pars);

#ifdef PRISMATIC_ENABLE_GPU
	cout << "peak GPU memory usage = " << prismatic_pars.maxGPUMem << '\n';
#endif //PRISMATIC_ENABLE_GPU
	std::cout << "PRISM Calculation complete.\n"
			  << std::endl;
	return prismatic_pars;
};

void PRISM_runFP(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, size_t fpNum)
{
	pars.meta.randomSeed = rand() % 100000;
	pars.meta.fpNum = fpNum;
	cout << "Frozen Phonon #" << fpNum << endl;
	pars.meta.toString();

	pars.outputFile = H5::H5File(pars.meta.filenameOutput.c_str(), H5F_ACC_RDWR);
	pars.fpFlag = fpNum;

	if(pars.meta.importSMatrix)
	{
		std::cout << "Skipping PRISM01. Using precalculated scattering matrix from: "  << pars.meta.importFile << std::endl;
	}
	else if(pars.meta.importPotential)
	{
		std::cout << "Using precalculated potential from " << pars.meta.importFile << std::endl;
		PRISM01_importPotential(pars);
	}
	else
	{
		PRISM01_calcPotential(pars);
	}

	// compute compact S-matrix
	if(pars.meta.importSMatrix)
	{
		PRISM02_importSMatrix(pars);
	}
	else
	{
		PRISM02_calcSMatrix(pars);
	}

	PRISM03_calcOutput(pars);
	pars.outputFile.close();

	if(fpNum >= 1)
	{
		pars.net_output += pars.output;
		pars.net_output_c += pars.output_c;
		if (pars.meta.saveDPC_CoM) pars.net_DPC_CoM += pars.DPC_CoM;
	}
	else
	{
		pars.net_output = pars.output;
		pars.net_output_c = pars.output_c;
		if (pars.meta.saveDPC_CoM) pars.net_DPC_CoM = pars.DPC_CoM;
	}
	
};

} // namespace Prismatic
