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

	//run multiple frozen phonons

	//convert beam indices into qx, qy indices by sorted tilts
	size_t N = prismatic_pars.numberBeams;
	std::vector<std::pair<PRISMATIC_FLOAT_PRECISION, PRISMATIC_FLOAT_PRECISION>> tilts(N);
	for(auto i = 0; i < N; i++) tilts[i] = std::make_pair(prismatic_pars.xTilts_tem[i], prismatic_pars.yTilts_tem[i]);
	std::vector<size_t> indices(N);
	for(auto i = 0; i < N; i++) indices[i] = i;
    std::sort(indices.begin(), indices.end(), [&](int i, int j){return tilts[i]<tilts[j];} );
	
	//sort tilt arrays to keep the same order in dim writing
	std::vector<PRISMATIC_FLOAT_PRECISION> xTilts_tmp(prismatic_pars.xTilts_tem);
	std::vector<PRISMATIC_FLOAT_PRECISION> yTilts_tmp(prismatic_pars.yTilts_tem);
	for(auto i = 0; i < N; i++)
	{
		prismatic_pars.xTilts_tem[i] = xTilts_tmp[indices[i]];
		prismatic_pars.yTilts_tem[i] = yTilts_tmp[indices[i]];
	}

	//save data
	prismatic_pars.outputFile = H5::H5File(prismatic_pars.meta.filenameOutput.c_str(), H5F_ACC_RDWR);
	std::cout << "Writing HRTEM data to output file." << std::endl;
	setupHRTEMOutput(prismatic_pars);
	
	//first scale the Smatrix back to mean intensity
	prismatic_pars.Scompact *= prismatic_pars.Scompact.get_dimi()*prismatic_pars.Scompact.get_dimj();
	
	//open group and write tilt images incrementally in loop
	//this avoids restriding S-matrix array all at once, which is memory intensive
	H5::DataSet hrtem_ds = prismatic_pars.outputFile.openDataSet("4DSTEM_simulation/data/realslices/HRTEM/realslice");

	std::array<size_t, 2> dims_in = {prismatic_pars.Scompact.get_dimi(), prismatic_pars.Scompact.get_dimj()};
	std::array<size_t, 2> rorder = {1, 0};
	size_t strides = prismatic_pars.Scompact.get_dimi()*prismatic_pars.Scompact.get_dimj();
	hsize_t offset[3] = {0,0,0};
	hsize_t mdims[3] = {prismatic_pars.Scompact.get_dimi(), prismatic_pars.Scompact.get_dimj(), 1};

	for(auto i = 0; i < N; i++) std::cout << indices[i] << std::endl;
	H5::DataSpace mspace(3, mdims);
	for(auto i = 0; i < N; i++)
	{
		offset[2] = i;
		Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> tmp_output = zeros_ND<2, std::complex<PRISMATIC_FLOAT_PRECISION>>({{prismatic_pars.Scompact.get_dimi(), prismatic_pars.Scompact.get_dimj()}});
		std::copy(&prismatic_pars.Scompact[indices[i]*strides], &prismatic_pars.Scompact[(indices[i]+1)*strides], tmp_output.begin());
		// std::copy(&prismatic_pars.Scompact[i*strides], &prismatic_pars.Scompact[(i+1)*strides], tmp_output.begin());
		tmp_output = restride(tmp_output, dims_in, rorder);
		if(prismatic_pars.meta.saveComplexOutputWave)
		{
			H5::DataSpace fspace = hrtem_ds.getSpace();
			fspace.selectHyperslab(H5S_SELECT_SET, mdims, offset);
			hrtem_ds.write(&tmp_output[0], hrtem_ds.getDataType(), mspace, fspace);
		}
		else
		{
			Array2D<PRISMATIC_FLOAT_PRECISION> tmp_output_int = zeros_ND<2, PRISMATIC_FLOAT_PRECISION>({{prismatic_pars.Scompact.get_dimi(), prismatic_pars.Scompact.get_dimj()}});
			for(auto j = 0; j < tmp_output.size(); j++) tmp_output_int[j] = pow(std::abs(tmp_output[j]), 2.0);
			H5::DataSpace fspace = hrtem_ds.getSpace();
			fspace.selectHyperslab(H5S_SELECT_SET, mdims, offset);
			hrtem_ds.write(&tmp_output_int[0], hrtem_ds.getDataType(), mspace, fspace);
		}
	}
	
    std::cout << "Calculation complete.\n" << std::endl;
    return prismatic_pars;
};
    
}