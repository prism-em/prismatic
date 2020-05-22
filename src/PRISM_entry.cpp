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
		PRISM01_importPotential(prismatic_pars);
	}
	else
	{
		PRISM01_calcPotential(prismatic_pars);
	}
	

	// compute compact S-matrix
	if(prismatic_pars.meta.importSMatrix)
	{
		PRISM02_importSMatrix(prismatic_pars);
	}
	else
	{
		PRISM02_calcSMatrix(prismatic_pars);
	}
	

	//		Array3D<PRISMATIC_FLOAT_PRECISION> tmp = zeros_ND<3, PRISMATIC_FLOAT_PRECISION>({{prismatic_pars.Scompact.get_dimk(),prismatic_pars.Scompact.get_dimj(),prismatic_pars.Scompact.get_dimi()}});
	//		Array3D<PRISMATIC_FLOAT_PRECISION> tmp_r = zeros_ND<3, PRISMATIC_FLOAT_PRECISION>({{prismatic_pars.Scompact.get_dimk(),prismatic_pars.Scompact.get_dimj(),prismatic_pars.Scompact.get_dimi()}});
	//		Array3D<PRISMATIC_FLOAT_PRECISION> tmp_i = zeros_ND<3, PRISMATIC_FLOAT_PRECISION>({{prismatic_pars.Scompact.get_dimk(),prismatic_pars.Scompact.get_dimj(),prismatic_pars.Scompact.get_dimi()}});
	//		auto tmp_ptr = tmp.begin();
	//		auto tmp_r_ptr = tmp_r.begin();
	//		auto tmp_i_ptr = tmp_i.begin();
	//		for (auto&i : prismatic_pars.Scompact)*tmp_ptr++ = abs(i);
	//		for (auto&i : prismatic_pars.Scompact)*tmp_r_ptr++ = i.real();
	//		for (auto&i : prismatic_pars.Scompact)*tmp_i_ptr++ = i.imag();
	//		std::complex<PRISMATIC_FLOAT_PRECISION> ssum{0,0};
	//		for (auto& i:prismatic_pars.Scompact)ssum+=i;
	//		cout <<"S compact sum = " << ssum << endl;
	//		tmp.toMRC_f("debug_scompact.mrc");
	//		tmp_r.toMRC_f("debug_scompact_r.mrc");
	//		tmp_i.toMRC_f("debug_scompact_i.mrc");

	// compute final output
	PRISM03_calcOutput(prismatic_pars);
	prismatic_pars.outputFile.close();

	// calculate remaining frozen phonon configurations
	if (prismatic_pars.meta.numFP > 1)
	{
		// run the rest of the frozen phonons
		Array4D<PRISMATIC_FLOAT_PRECISION> net_output(prismatic_pars.output);
		Array4D<PRISMATIC_FLOAT_PRECISION> DPC_CoM_output;
		if (prismatic_pars.meta.saveDPC_CoM)
			DPC_CoM_output = prismatic_pars.DPC_CoM;
		for (auto fp_num = 1; fp_num < prismatic_pars.meta.numFP; ++fp_num)
		{
			meta.randomSeed = rand() % 100000;
			++meta.fpNum;
			Parameters<PRISMATIC_FLOAT_PRECISION> prismatic_pars(meta);
			cout << "Frozen Phonon #" << fp_num << endl;
			prismatic_pars.meta.toString();

			prismatic_pars.outputFile = H5::H5File(prismatic_pars.meta.filenameOutput.c_str(), H5F_ACC_RDWR);
			prismatic_pars.fpFlag = fp_num;

			PRISM01_calcPotential(prismatic_pars);
			PRISM02_calcSMatrix(prismatic_pars);
			PRISM03_calcOutput(prismatic_pars);
			net_output += prismatic_pars.output;
			if (meta.saveDPC_CoM)
				DPC_CoM_output += prismatic_pars.DPC_CoM;
			prismatic_pars.outputFile.close();
		}
		// divide to take average
		for (auto &i : net_output)
			i /= prismatic_pars.meta.numFP;
		prismatic_pars.output = net_output;

		if (prismatic_pars.meta.saveDPC_CoM)
		{
			for (auto &j : DPC_CoM_output)
				j /= prismatic_pars.meta.numFP; //since squared intensities are used to calculate DPC_CoM, this is incoherent averaging
			prismatic_pars.DPC_CoM = DPC_CoM_output;
		}
	}

	prismatic_pars.outputFile = H5::H5File(prismatic_pars.meta.filenameOutput.c_str(), H5F_ACC_RDWR);

	if (prismatic_pars.meta.save3DOutput)
	{
		PRISMATIC_FLOAT_PRECISION dummy = 1.0;
		setupVDOutput(prismatic_pars, prismatic_pars.output.get_diml(), dummy);
		Array3D<PRISMATIC_FLOAT_PRECISION> output_image = zeros_ND<3, PRISMATIC_FLOAT_PRECISION>({{prismatic_pars.output.get_dimj(), prismatic_pars.output.get_dimk(), prismatic_pars.output.get_dimi()}});

		std::stringstream nameString;
		nameString << "4DSTEM_simulation/data/realslices/virtual_detector_depth" << getDigitString(0);
		H5::Group dataGroup = prismatic_pars.outputFile.openGroup(nameString.str());

		std::string dataSetName = "realslice";
		H5::DataSet VD_data = dataGroup.openDataSet(dataSetName);
		hsize_t mdims[3] = {prismatic_pars.xp.size(), prismatic_pars.yp.size(), prismatic_pars.Ndet};

		for (auto b = 0; b < prismatic_pars.Ndet; b++)
		{
			for (auto y = 0; y < prismatic_pars.output.get_dimk(); ++y)
			{
				for (auto x = 0; x < prismatic_pars.output.get_dimj(); ++x)
				{
					output_image.at(x, y, b) = prismatic_pars.output.at(0, y, x, b);
				}
			}
		}
		writeDatacube3D(VD_data, &output_image[0], mdims);
		VD_data.close();
		dataGroup.close();
	}




	if (prismatic_pars.meta.save2DOutput)
	{
		size_t lower = std::max((size_t)0, (size_t)(prismatic_pars.meta.integrationAngleMin / prismatic_pars.meta.detectorAngleStep));
		size_t upper = std::min((size_t)prismatic_pars.detectorAngles.size(), (size_t)(prismatic_pars.meta.integrationAngleMax / prismatic_pars.meta.detectorAngleStep));
		Array2D<PRISMATIC_FLOAT_PRECISION> prism_image;
		prism_image = zeros_ND<2, PRISMATIC_FLOAT_PRECISION>(
			{{prismatic_pars.output.get_dimj(), prismatic_pars.output.get_dimk()}});
		PRISMATIC_FLOAT_PRECISION dummy = 1.0;
		setup2DOutput(prismatic_pars, prismatic_pars.output.get_diml(), dummy);

		std::stringstream nameString;
		nameString << "4DSTEM_simulation/data/realslices/annular_detector_depth" << getDigitString(0);
		H5::Group dataGroup = prismatic_pars.outputFile.openGroup(nameString.str());
		H5::DataSet AD_data = dataGroup.openDataSet("realslice");
		hsize_t mdims[2] = {prismatic_pars.xp.size(), prismatic_pars.yp.size()};

		for (auto y = 0; y < prismatic_pars.output.get_dimk(); ++y)
		{
			for (auto x = 0; x < prismatic_pars.output.get_dimj(); ++x)
			{
				for (auto b = lower; b < upper; ++b)
				{
					prism_image.at(x, y) += prismatic_pars.output.at(0, y, x, b);
                }
            }
        }

		writeRealSlice(AD_data, &prism_image[0], mdims);
		AD_data.close();
		dataGroup.close();
	}


	if (prismatic_pars.meta.saveDPC_CoM)
	{
		PRISMATIC_FLOAT_PRECISION dummy = 1.0;
		setupDPCOutput(prismatic_pars, prismatic_pars.output.get_diml(), dummy);

		//create dummy array to pass to
		Array3D<PRISMATIC_FLOAT_PRECISION> DPC_slice;

		std::stringstream nameString;
		nameString << "4DSTEM_simulation/data/realslices/DPC_CoM_depth" << getDigitString(0);
		H5::Group dataGroup = prismatic_pars.outputFile.openGroup(nameString.str());
		hsize_t mdims[3] = {prismatic_pars.xp.size(), prismatic_pars.yp.size(), 2};

		DPC_slice = zeros_ND<3, PRISMATIC_FLOAT_PRECISION>({{prismatic_pars.DPC_CoM.get_dimj(), prismatic_pars.DPC_CoM.get_dimk(), 2}});
		std::string dataSetName = "realslice";
		H5::DataSet DPC_data = dataGroup.openDataSet(dataSetName);

		for (auto b = 0; b < prismatic_pars.DPC_CoM.get_dimi(); ++b)
		{
			for (auto y = 0; y < prismatic_pars.DPC_CoM.get_dimk(); ++y)
			{
				for (auto x = 0; x < prismatic_pars.DPC_CoM.get_dimj(); ++x)
				{
					DPC_slice.at(x, y, b) = prismatic_pars.DPC_CoM.at(0, y, x, b);
				}
			}
		}

		writeDatacube3D(DPC_data, &DPC_slice[0], mdims);
		DPC_data.close();
		dataGroup.close();
	}

	PRISMATIC_FLOAT_PRECISION dummy = 1.0;
	writeMetadata(prismatic_pars, dummy);
	prismatic_pars.outputFile.close();

#ifdef PRISMATIC_ENABLE_GPU
	cout << "peak GPU memory usage = " << prismatic_pars.maxGPUMem << '\n';
#endif //PRISMATIC_ENABLE_GPU
	std::cout << "PRISM Calculation complete.\n"
			  << std::endl;
	return prismatic_pars;
}
} // namespace Prismatic
