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


#include "meta.h"
#include "params.h"
#include "ArrayND.h"
#include "configure.h"
#include <algorithm>
#include "HRTEM_entry.h"
#include "PRISM01_calcPotential.h"
#include "PRISM02_calcSMatrix.h"
#include "PRISM03_calcOutput.h"
#include "fileIO.h"

namespace Prismatic
{

Parameters<PRISMATIC_FLOAT_PRECISION> HRTEM_entry(Metadata<PRISMATIC_FLOAT_PRECISION>& meta)
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
    prismatic_pars.outputFile.close();

    // compute projected potentials
    prismatic_pars.fpFlag = 0;

    if(prismatic_pars.meta.importPotential) configureImportFP(prismatic_pars);
    prismatic_pars.scale = 1.0;

	Array3D<PRISMATIC_FLOAT_PRECISION> net_output;
	if(not prismatic_pars.meta.saveComplexOutputWave)
	{
		//integrate output here
		net_output = zeros_ND<3, PRISMATIC_FLOAT_PRECISION>(prismatic_pars.Scompact.get_dimarr());
		for(auto i = 0; i < net_output.size(); i++)
		{
			//scale before integration
			prismatic_pars.Scompact *= prismatic_pars.Scompact.get_dimi()*prismatic_pars.Scompact.get_dimj();
			net_output[i] += pow(std::abs(prismatic_pars.Scompact[i]), 2.0) / prismatic_pars.meta.numFP;
		}
	}

	//run multiple frozen phonons
	for(auto i = 0; i < prismatic_pars.meta.numFP; i++)
	{
		HRTEM_runFP(prismatic_pars, i);

		//apply aberrations
		if(prismatic_pars.meta.aberrations.size() > 0)
		{
			// setup necessary coordinates
			setupCoordinates_2(prismatic_pars);
			setupDetector(prismatic_pars);
			setupFourierCoordinates(prismatic_pars);
			transformIndices(prismatic_pars);
		
			//now apply aberrations to each beam
			apply_aberrations(prismatic_pars);
		}

		if(prismatic_pars.meta.saveComplexOutputWave)
		{			
			//save FP individually
			prismatic_pars.outputFile = H5::H5File(prismatic_pars.meta.filenameOutput.c_str(), H5F_ACC_RDWR);
			std::cout << "Writing HRTEM data to output file." << std::endl;
			sortHRTEMbeams(prismatic_pars);
			setupHRTEMOutput(prismatic_pars);
			setupHRTEMOutput_virtual(prismatic_pars);
			saveHRTEM(prismatic_pars, net_output);
			prismatic_pars.outputFile.close();
		}
		else
		{
			//integrate output
			PRISMATIC_FLOAT_PRECISION scale = prismatic_pars.Scompact.get_dimj() * prismatic_pars.Scompact.get_dimi();
			for(auto i = 0; i < net_output.size(); i++)
			{
				net_output[i] += pow(std::abs(prismatic_pars.Scompact[i]*scale), 2.0) / prismatic_pars.meta.numFP;
			}
		}
	}


	//save data
	prismatic_pars.outputFile = H5::H5File(prismatic_pars.meta.filenameOutput.c_str(), H5F_ACC_RDWR);
	if(not prismatic_pars.meta.saveComplexOutputWave)
	{
		std::cout << "Writing HRTEM data to output file." << std::endl;
		sortHRTEMbeams(prismatic_pars);
		setupHRTEMOutput(prismatic_pars);
		setupHRTEMOutput_virtual(prismatic_pars);
		saveHRTEM(prismatic_pars, net_output);
	};
	
    std::cout << "Calculation complete.\n" << std::endl;
	writeMetadata(prismatic_pars);
	prismatic_pars.outputFile.close();
    return prismatic_pars;
};

void HRTEM_runFP(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, size_t fpNum)
{
	pars.meta.randomSeed = rand() % 100000;
	pars.meta.fpNum = fpNum;
	std::cout << "Frozen Phonon #" << fpNum << std::endl;
	pars.meta.toString();

	pars.outputFile = H5::H5File(pars.meta.filenameOutput.c_str(), H5F_ACC_RDWR);
	pars.fpFlag = fpNum;
	
	//update aberrations here to maintain consistency between runFP calls
	pars.meta.aberrations = updateAberrations(pars.meta.aberrations, pars.meta.probeDefocus, pars.meta.C3, pars.meta.C5);
	if(pars.meta.importPotential)
	{
		std::cout << "Using precalculated potential from " << pars.meta.importFile << std::endl;
		PRISM01_importPotential(pars);
	}
	else
	{
		PRISM01_calcPotential(pars);
	}

	PRISM02_calcSMatrix(pars);
};
    
}