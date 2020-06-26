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

    // compute projected potentials
    prismatic_pars.fpFlag = 0;

    if(prismatic_pars.meta.importPotential)
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
		//update original object as prismatic_pars is recreated later
		meta.numFP = prismatic_pars.meta.numFP;
		PRISM01_importPotential(prismatic_pars);
	}
	else
	{
		PRISM01_calcPotential(prismatic_pars);
	}

    prismatic_pars.scale = 1.0;

    PRISM02_calcSMatrix(prismatic_pars);
    prismatic_pars.outputFile.close();

	Array3D<PRISMATIC_FLOAT_PRECISION> net_output;
	Array3D<std::complex<PRISMATIC_FLOAT_PRECISION>> net_output_c;
	if(prismatic_pars.meta.saveComplexOutputWave)
	{
		//scale to mean intensity and average by FP
		net_output_c = prismatic_pars.Scompact*prismatic_pars.Scompact.get_dimi()*prismatic_pars.Scompact.get_dimj()/ prismatic_pars.meta.numFP;
	}
	else
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
	if(prismatic_pars.meta.numFP > 1)
	{
		for(auto fp_num = 1; fp_num < prismatic_pars.meta.numFP; ++fp_num)
		{
			meta.randomSeed = rand() % 100000;
			++meta.fpNum;
			Parameters<PRISMATIC_FLOAT_PRECISION> prismatic_pars(meta);
			std::cout << "Frozen Phonon #" << fp_num << std::endl;
			prismatic_pars.meta.toString();

			prismatic_pars.outputFile = H5::H5File(prismatic_pars.meta.filenameOutput.c_str(), H5F_ACC_RDWR);
			prismatic_pars.fpFlag = fp_num;

			if(prismatic_pars.meta.importPotential)
			{
				std::cout << "Using precalculated potential from " << prismatic_pars.meta.importFile << std::endl;
				PRISM01_importPotential(prismatic_pars);
			}
			else
			{
				PRISM01_calcPotential(prismatic_pars);
			}

			PRISM02_calcSMatrix(prismatic_pars);

			if(prismatic_pars.meta.saveComplexOutputWave)
			{
				net_output_c += prismatic_pars.Scompact / prismatic_pars.meta.numFP;
			}
			else
			{
				for(auto i = 0; i < net_output.size(); i++)
				{
					net_output[i] += pow(std::abs(prismatic_pars.Scompact[i]), 2.0) / prismatic_pars.meta.numFP;
				}
			}
			
		}

	}

	//save data
	prismatic_pars.outputFile = H5::H5File(prismatic_pars.meta.filenameOutput.c_str(), H5F_ACC_RDWR);
	std::cout << "Writing HRTEM data to output file." << std::endl;
	sortHRTEMbeams(prismatic_pars);
	setupHRTEMOutput(prismatic_pars);
	setupHRTEMOutput_virtual(prismatic_pars);
	saveHRTEM(prismatic_pars, net_output_c, net_output);
	
    std::cout << "Calculation complete.\n" << std::endl;
    return prismatic_pars;
};
    
}