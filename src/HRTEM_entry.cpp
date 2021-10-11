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

void HRTEM_entry(Metadata<PRISMATIC_FLOAT_PRECISION> &meta)
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

	HRTEM_entry_pars(pars);
}

void HRTEM_entry_pars(Parameters<PRISMATIC_FLOAT_PRECISION>& pars)
{

    pars.meta.toString();
    
    pars.outputFile = H5::H5File(pars.meta.filenameOutput.c_str(), H5F_ACC_TRUNC);
    setupOutputFile(pars);
    pars.outputFile.close();

    // compute projected potentials
    pars.fpFlag = 0;

    if(pars.meta.importPotential) configureImportFP(pars);
    pars.scale = 1.0;

	Array3D<PRISMATIC_FLOAT_PRECISION> net_output;

	//run multiple frozen phonons
	for(auto i = 0; i < pars.meta.numFP; i++)
	{
		HRTEM_runFP(pars, i);

		if(pars.meta.saveComplexOutputWave)
		{			
			//save FP individually
			pars.outputFile = H5::H5File(pars.meta.filenameOutput.c_str(), H5F_ACC_RDWR);
			std::cout << "Writing HRTEM data to output file." << std::endl;
			sortHRTEMbeams(pars);
			setupHRTEMOutput(pars);
			setupHRTEMOutput_virtual(pars);
			saveHRTEM(pars, net_output);
			pars.outputFile.close();
		}
		else
		{
			if(i == 0)
			{
				net_output = zeros_ND<3, PRISMATIC_FLOAT_PRECISION>({{pars.Scompact.get_dimi(), pars.Scompact.get_dimj(), pars.Scompact.get_dimk()}});
			}

			//integrate output
            sortHRTEMbeams(pars);
			PRISMATIC_FLOAT_PRECISION scale = pars.Scompact.get_dimj() * pars.Scompact.get_dimi();
			for(auto kk = 0; kk < pars.Scompact.get_dimk(); kk++)
			{
				for(auto jj = 0; jj < pars.Scompact.get_dimj(); jj++)
				{
					for(auto ii = 0; ii < pars.Scompact.get_dimi(); ii++)
					{
						net_output.at(ii,jj,kk) += pow(std::abs(pars.Scompact.at(pars.HRTEMbeamOrder[kk],jj,ii)*scale), 2.0) / pars.meta.numFP;
					}
				}
			}

		}
	}


	//save data
	pars.outputFile = H5::H5File(pars.meta.filenameOutput.c_str(), H5F_ACC_RDWR);
	if(not pars.meta.saveComplexOutputWave)
	{
		std::cout << "Writing HRTEM data to output file." << std::endl;
		setupHRTEMOutput(pars);
		setupHRTEMOutput_virtual(pars);
		saveHRTEM(pars, net_output);
	};
	
    std::cout << "Calculation complete.\n" << std::endl;
	save_qArr(pars);
	writeMetadata(pars);
	pars.outputFile.close();
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
	pars.meta.aberrations = updateAberrations(pars.meta.aberrations, pars.meta.probeDefocus, pars.meta.C3, pars.meta.C5, pars.lambda);
	if(!pars.potentialReady){
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

	PRISM02_calcSMatrix(pars);
};
    
}