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
#include "aberration.h"

namespace Prismatic
{
using namespace std;
void PRISM_entry(Metadata<PRISMATIC_FLOAT_PRECISION> &meta)
{
	Parameters<PRISMATIC_FLOAT_PRECISION> pars;
	try
	{ // read atomic coordinates
		pars = Parameters<PRISMATIC_FLOAT_PRECISION>(meta);
	}
	catch (...)
	{
		std::cout << "Terminating" << std::endl;
		exit(1);
	}

    PRISM_entry_pars(pars);
}

void PRISM_entry_pars(Parameters<PRISMATIC_FLOAT_PRECISION> &pars)
{
	pars.outputFile = H5::H5File(pars.meta.filenameOutput.c_str(), H5F_ACC_TRUNC);
	setupOutputFile(pars);

	if(pars.meta.importPotential or pars.meta.importSMatrix) configureImportFP(pars);

	if(pars.meta.simSeries)
	{
		for(auto i = 0; i < pars.meta.numFP; i++)
		{
			PRISM_series_runFP(pars, i);
		}

		for(auto i = 0; i < pars.meta.seriesTags.size(); i++)
		{
			std::string currentName = pars.meta.seriesTags[i];
			pars.currentTag = currentName;
			pars.meta.probeDefocus = pars.meta.seriesVals[0][i]; //TODO: later, if expanding sim series past defocus, need to pull current val more generally

			readRealDataSet_inOrder(pars.net_output, "prismatic_scratch.h5", "scratch/"+currentName);
			if(pars.meta.saveDPC_CoM)
				readRealDataSet_inOrder(pars.net_DPC_CoM, "prismatic_scratch.h5", "scratch/"+currentName+"_DPC");
			//average data by fp
			for (auto &i : pars.net_output)
				i /= pars.meta.numFP;

			if (pars.meta.saveDPC_CoM)
			{
				for (auto &j : pars.net_DPC_CoM)
					j /= pars.meta.numFP; //since squared intensities are used to calculate DPC_CoM, this is incoherent averaging
			}

			saveSTEM(pars);
		}
	}
	else
	{
		pars.meta.aberrations = updateAberrations(pars.meta.aberrations, pars.meta.probeDefocus, pars.meta.C3, pars.meta.C5, pars.lambda);
		for(auto i = 0; i < pars.meta.numFP; i++)
		{
			PRISM_runFP(pars, i);
		}

		std::cout << "All frozen phonon configurations complete. Writing data to output file." << std::endl;
		//average data by fp
		for (auto &i : pars.net_output)
			i /= pars.meta.numFP;

		if (pars.meta.saveDPC_CoM)
		{
			for (auto &j : pars.net_DPC_CoM)
				j /= pars.meta.numFP; //since squared intensities are used to calculate DPC_CoM, this is incoherent averaging
		}

		saveSTEM(pars);
	}
	
	pars.outputFile = H5::H5File(pars.meta.filenameOutput.c_str(), H5F_ACC_RDWR);
	
	//perhaps have this check against the keys
	if(pars.meta.simSeries) CCseriesSG(pars.outputFile);

	writeMetadata(pars);
	pars.outputFile.close();

	if (pars.meta.simSeries) removeScratchFile(pars);

#ifdef PRISMATIC_ENABLE_GPU
	cout << "peak GPU memory usage = " << pars.maxGPUMem << '\n';
#endif //PRISMATIC_ENABLE_GPU
	std::cout << "PRISM Calculation complete.\n"
			  << std::endl;
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
	else if(!pars.potentialReady){
        if(pars.meta.importPotential)
        {
            std::cout << "Using precalculated potential from " << pars.meta.importFile << std::endl;
            PRISM01_importPotential(pars);
        }
        else
        {
            PRISM01_calcPotential(pars);
        }
    }
    else{
        //reset flag if more than one FP
		if(pars.meta.numFP > 1) pars.potentialReady = false;
    }

#ifdef PRISMATIC_BUILDING_GUI
    pars.parent_thread->passPotentialToParent(pars.pot);
#endif

	// compute compact S-matrix
	if(pars.meta.importSMatrix)
	{
		PRISM02_importSMatrix(pars);
	}
	else
	{
		PRISM02_calcSMatrix(pars);
	}

	if(pars.meta.matrixRefocus)
	{
		refocus(pars);
	}


	PRISM03_calcOutput(pars);
	pars.outputFile.close();

	if(fpNum >= 1)
	{
		pars.net_output += pars.output;
		if (pars.meta.saveDPC_CoM) pars.net_DPC_CoM += pars.DPC_CoM;
	}
	else
	{
		pars.net_output = pars.output;
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
	else if(!pars.potentialReady){
        if(pars.meta.importPotential)
        {
            std::cout << "Using precalculated potential from " << pars.meta.importFile << std::endl;
            PRISM01_importPotential(pars);
        }
        else
        {
            PRISM01_calcPotential(pars);
        }
    }
    else{
        //reset flag if more than one FP
		if(pars.meta.numFP > 1) pars.potentialReady = false;
    }

#ifdef PRISMATIC_BUILDING_GUI
    pars.parent_thread->passPotentialToParent(pars.pot);
#endif

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
		pars.meta.aberrations = updateAberrations(pars.meta.aberrations, pars.meta.probeDefocus, pars.meta.C3, pars.meta.C5, pars.lambda);
		//need to use current shift-- so the matrix doesn't get refocused out to oblivion
		if(pars.meta.matrixRefocus)
		{
			refocus(pars);
		}
		PRISM03_calcOutput(pars);

		if(i == 0 and fpNum == 0) createScratchFile(pars);
		updateScratchData(pars);
	}
	pars.outputFile.close();

};

} // namespace Prismatic
