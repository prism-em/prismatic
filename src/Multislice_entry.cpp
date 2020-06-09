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
#include "Multislice_calcOutput.h"
#include "PRISM01_calcPotential.h"
#include "PRISM02_calcSMatrix.h"
#include <algorithm>
#include "utility.h"
#include "fileIO.h"

namespace Prismatic
{
Parameters<PRISMATIC_FLOAT_PRECISION> Multislice_entry(Metadata<PRISMATIC_FLOAT_PRECISION> &meta)
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
	// compute final output
	Multislice_calcOutput(prismatic_pars);
	prismatic_pars.outputFile.close();

	// calculate remaining frozen phonon configurations
	//TODO: Clarify the scope issues occuring here. Extraneous copy of prismatic_pars structure?
	if (prismatic_pars.meta.numFP > 1)
	{
		// run the rest of the frozen phonons
		Array4D<PRISMATIC_FLOAT_PRECISION> net_output(prismatic_pars.output);
		Array4D<std::complex<PRISMATIC_FLOAT_PRECISION>> net_output_c(prismatic_pars.output_c);
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
			prismatic_pars.scale = 1.0;

			if(prismatic_pars.meta.importPotential)
			{
				std::cout << "Using precalculated potential from " << prismatic_pars.meta.importFile << std::endl;
				PRISM01_importPotential(prismatic_pars);
			}
			else
			{
				PRISM01_calcPotential(prismatic_pars);
			}

			Multislice_calcOutput(prismatic_pars);
			net_output += prismatic_pars.output;
			net_output_c += prismatic_pars.output_c;
			if (meta.saveDPC_CoM)
				DPC_CoM_output += prismatic_pars.DPC_CoM;
			prismatic_pars.outputFile.close();
		}
		// divide to take average
		for (auto &i : net_output)
			i /= prismatic_pars.meta.numFP;
		for (auto &i : net_output_c)
			i /= prismatic_pars.meta.numFP;
		prismatic_pars.output = net_output;
		prismatic_pars.output_c = net_output_c;

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
		if(prismatic_pars.meta.saveComplexOutputWave)
		{
			setupVDOutput(prismatic_pars, prismatic_pars.output_c.get_diml());
			//create dummy array to pass to
			Array3D<std::complex<PRISMATIC_FLOAT_PRECISION>> slice_image;
			slice_image = zeros_ND<3, std::complex<PRISMATIC_FLOAT_PRECISION>>({{prismatic_pars.output_c.get_dimj(), prismatic_pars.output_c.get_dimk(), prismatic_pars.output_c.get_dimi()}});

			for (auto j = 0; j < prismatic_pars.output_c.get_diml(); j++)
			{
				std::stringstream nameString;
				nameString << "4DSTEM_simulation/data/realslices/virtual_detector_depth" << getDigitString(j);
				H5::Group dataGroup = prismatic_pars.outputFile.openGroup(nameString.str());
				hsize_t mdims[3] = {prismatic_pars.xp.size(), prismatic_pars.yp.size(), prismatic_pars.Ndet};

				std::string dataSetName = "realslice";
				// H5::DataSet VD_data = dataGroup.openDataSet(dataSetName);
				for (auto b = 0; b < prismatic_pars.output_c.get_dimi(); ++b)
				{
					for (auto y = 0; y < prismatic_pars.output_c.get_dimk(); ++y)
					{
						for (auto x = 0; x < prismatic_pars.output_c.get_dimj(); ++x)
						{
							slice_image.at(x, y, b) = prismatic_pars.output_c.at(j, y, x, b);
						}
					}
				}

				writeComplexDataSet(dataGroup, "realslice", &slice_image[0], mdims, 3);
				dataGroup.close();
			}
		}
		else
		{
			setupVDOutput(prismatic_pars, prismatic_pars.output.get_diml());
			//create dummy array to pass to
			Array3D<PRISMATIC_FLOAT_PRECISION> slice_image;
			slice_image = zeros_ND<3, PRISMATIC_FLOAT_PRECISION>({{prismatic_pars.output.get_dimj(), prismatic_pars.output.get_dimk(), prismatic_pars.output.get_dimi()}});

			for (auto j = 0; j < prismatic_pars.output.get_diml(); j++)
			{
				std::stringstream nameString;
				nameString << "4DSTEM_simulation/data/realslices/virtual_detector_depth" << getDigitString(j);
				H5::Group dataGroup = prismatic_pars.outputFile.openGroup(nameString.str());
				hsize_t mdims[3] = {prismatic_pars.xp.size(), prismatic_pars.yp.size(), prismatic_pars.Ndet};

				std::string dataSetName = "realslice";
				H5::DataSet VD_data = dataGroup.openDataSet(dataSetName);
				for (auto b = 0; b < prismatic_pars.output.get_dimi(); ++b)
				{
					for (auto y = 0; y < prismatic_pars.output.get_dimk(); ++y)
					{
						for (auto x = 0; x < prismatic_pars.output.get_dimj(); ++x)
						{
							slice_image.at(x, y, b) = prismatic_pars.output.at(j, y, x, b);
						}
					}
				}

				writeDatacube3D(VD_data, &slice_image[0], mdims);
				VD_data.close();
				dataGroup.close();
			}
		}
	}

	if (prismatic_pars.meta.save2DOutput)
	{
		size_t lower = std::max((size_t)0, (size_t)(prismatic_pars.meta.integrationAngleMin / prismatic_pars.meta.detectorAngleStep));
		size_t upper = std::min(prismatic_pars.detectorAngles.size(), (size_t)(prismatic_pars.meta.integrationAngleMax / prismatic_pars.meta.detectorAngleStep));
		if(prismatic_pars.meta.saveComplexOutputWave)
		{
			setup2DOutput(prismatic_pars, prismatic_pars.output_c.get_diml());
			Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> prism_image;
			for (auto j = 0; j < prismatic_pars.output_c.get_diml(); j++)
			{
				//need to initiliaze output image at each slice to prevent overflow of value
				prism_image = zeros_ND<2, std::complex<PRISMATIC_FLOAT_PRECISION>>(
					{{prismatic_pars.output_c.get_dimj(), prismatic_pars.output_c.get_dimk()}});

				for (auto y = 0; y < prismatic_pars.output_c.get_dimk(); ++y)
				{
					for (auto x = 0; x < prismatic_pars.output_c.get_dimj(); ++x)
					{
						for (auto b = lower; b < upper; ++b)
						{
							prism_image.at(x, y) += prismatic_pars.output_c.at(j, y, x, b);
						}
					}
				}

				std::stringstream nameString;
				nameString << "4DSTEM_simulation/data/realslices/annular_detector_depth" << getDigitString(j);
				H5::Group dataGroup = prismatic_pars.outputFile.openGroup(nameString.str());
				hsize_t mdims[2] = {prismatic_pars.xp.size(), prismatic_pars.yp.size()};

				writeComplexDataSet(dataGroup, "realslice", &prism_image[0], mdims, 2);
				dataGroup.close();
			}
		}
		else
		{
			setup2DOutput(prismatic_pars, prismatic_pars.output.get_diml());
			Array2D<PRISMATIC_FLOAT_PRECISION> prism_image;
			for (auto j = 0; j < prismatic_pars.output.get_diml(); j++)
			{
				//need to initiliaze output image at each slice to prevent overflow of value
				prism_image = zeros_ND<2, PRISMATIC_FLOAT_PRECISION>(
					{{prismatic_pars.output.get_dimj(), prismatic_pars.output.get_dimk()}});

				for (auto y = 0; y < prismatic_pars.output.get_dimk(); ++y)
				{
					for (auto x = 0; x < prismatic_pars.output.get_dimj(); ++x)
					{
						for (auto b = lower; b < upper; ++b)
						{
							prism_image.at(x, y) += prismatic_pars.output.at(j, y, x, b);
						}
					}
				}
				//if (prismatic_pars.meta.numSlices != 0){
				//	image_filename = prismatic_pars.meta.outputFolder +  (std::string("multislice_2Doutput_slice") + std::to_string(j) + std::string("_")) + prismatic_pars.meta.filenameOutput;
				//}
				//prism_image.toMRC_f(image_filename.c_str());
				std::stringstream nameString;
				nameString << "4DSTEM_simulation/data/realslices/annular_detector_depth" << getDigitString(j);
				H5::Group dataGroup = prismatic_pars.outputFile.openGroup(nameString.str());
				H5::DataSet AD_data = dataGroup.openDataSet("realslice");
				hsize_t mdims[2] = {prismatic_pars.xp.size(), prismatic_pars.yp.size()};

				writeRealSlice(AD_data, &prism_image[0], mdims);
				AD_data.close();
				dataGroup.close();
			}
		}
	}

	if (prismatic_pars.meta.saveDPC_CoM and not prismatic_pars.meta.saveComplexOutputWave)
	{
		setupDPCOutput(prismatic_pars, prismatic_pars.output.get_diml());

		//create dummy array to pass to
		Array3D<PRISMATIC_FLOAT_PRECISION> DPC_slice;
		DPC_slice = zeros_ND<3, PRISMATIC_FLOAT_PRECISION>({{prismatic_pars.DPC_CoM.get_dimj(), prismatic_pars.DPC_CoM.get_dimk(), 2}});
		hsize_t mdims[3] = {prismatic_pars.xp.size(), prismatic_pars.yp.size(), 2};

		for (auto j = 0; j < prismatic_pars.output.get_diml(); j++)
		{
			std::stringstream nameString;
			nameString << "4DSTEM_simulation/data/realslices/DPC_CoM_depth" << getDigitString(j);
			H5::Group dataGroup = prismatic_pars.outputFile.openGroup(nameString.str());

			std::string dataSetName = "realslice";

			H5::DataSet DPC_data = dataGroup.openDataSet(dataSetName);

			for (auto b = 0; b < prismatic_pars.DPC_CoM.get_dimi(); ++b)
			{
				for (auto y = 0; y < prismatic_pars.DPC_CoM.get_dimk(); ++y)
				{
					for (auto x = 0; x < prismatic_pars.DPC_CoM.get_dimj(); ++x)
					{
						DPC_slice.at(x, y, b) = prismatic_pars.DPC_CoM.at(j, y, x, b);
					}
				}
			}

			writeDatacube3D(DPC_data, &DPC_slice[0], mdims);
			DPC_data.close();
			dataGroup.close();
		}
	}

	writeMetadata(prismatic_pars);
	prismatic_pars.outputFile.close();

#ifdef PRISMATIC_ENABLE_GPU
	cout << "peak GPU memory usage = " << prismatic_pars.maxGPUMem << '\n';
#endif //PRISMATIC_ENABLE_GPU
	std::cout << "Calculation complete.\n"
			  << std::endl;
	return prismatic_pars;
}
} // namespace Prismatic