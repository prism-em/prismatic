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

	prismatic_pars.outputFile = H5::H5File(prismatic_pars.meta.filenameOutput.c_str(), H5F_ACC_TRUNC);
	setupOutputFile(prismatic_pars);

	if(prismatic_pars.meta.importPotential or prismatic_pars.meta.importSMatrix) configureImportFP(prismatic_pars);

	if(prismatic_pars.meta.simSeries)
	{
		for(auto i = 0; i < prismatic_pars.meta.numFP; i++)
		{
			PRISM_series_runFP(prismatic_pars, i);
		}

		for(auto i = 0; i < prismatic_pars.meta.seriesTags.size(); i++)
		{
			std::cout << "writing output for series iter " << i << std::endl;
			std::string currentName = prismatic_pars.meta.seriesTags[i];
			prismatic_pars.currentTag = currentName;
			readRealDataSet_inOrder(prismatic_pars.net_output, "prismatic_scratch.h5", "scratch/"+currentName);
			if(prismatic_pars.meta.saveDPC_CoM)
				readRealDataSet_inOrder(prismatic_pars.net_DPC_CoM, "prismatic_scratch.h5", "scratch/"+currentName+"_DPC");
			//average data by fp
			for (auto &i : prismatic_pars.net_output)
				i /= prismatic_pars.meta.numFP;

			if (prismatic_pars.meta.saveDPC_CoM)
			{
				for (auto &j : prismatic_pars.net_DPC_CoM)
					j /= prismatic_pars.meta.numFP; //since squared intensities are used to calculate DPC_CoM, this is incoherent averaging
			}

			saveSTEM(prismatic_pars);
		}
	}
	else
	{
		for(auto i = 0; i < prismatic_pars.meta.numFP; i++)
		{
			PRISM_runFP(prismatic_pars, i);
		}

		std::cout << "All frozen phonon configurations complete. Writing data to output file." << std::endl;
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
	}
	
	prismatic_pars.outputFile = H5::H5File(prismatic_pars.meta.filenameOutput.c_str(), H5F_ACC_RDWR);
	
	//perhaps have this check against the keys
	if(prismatic_pars.meta.simSeries) CCseriesSG(prismatic_pars.outputFile);

	writeMetadata(prismatic_pars);
	prismatic_pars.outputFile.close();
	removeScratchFile(prismatic_pars);

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

void PRISM_series_runFP(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, size_t fpNum)
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

	for(auto i = 0; i < pars.meta.seriesVals[0].size(); i++)
	{
		std::cout << "------------------- Series iter " << i << " -------------------" << std::endl;
		updateSeriesParams(pars, i);
		PRISM03_calcOutput(pars);

		if(i == 0 and fpNum == 0) createScratchFile(pars);
		updateScratchData(pars);
	}
	pars.outputFile.close();

};

} // namespace Prismatic
