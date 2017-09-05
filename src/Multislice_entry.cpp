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


namespace Prismatic{
	Parameters<PRISMATIC_FLOAT_PRECISION> Multislice_entry(Metadata<PRISMATIC_FLOAT_PRECISION>& meta){
		Parameters<PRISMATIC_FLOAT_PRECISION> prismatic_pars;
		try { // read atomic coordinates
			prismatic_pars = Parameters<PRISMATIC_FLOAT_PRECISION>(meta);
		} catch(...){
			std::cout << "Terminating" << std::endl;
			exit(1);
		}

		// compute projected potentials
		PRISM01_calcPotential(prismatic_pars);

		// compute final output
		Multislice_calcOutput(prismatic_pars);

		// calculate remaining frozen phonon configurations
		if (prismatic_pars.meta.numFP > 1) {
			// run the rest of the frozen phonons
			Array3D<PRISMATIC_FLOAT_PRECISION> net_output(prismatic_pars.output);
			for (auto fp_num = 1; fp_num < prismatic_pars.meta.numFP; ++fp_num){
				meta.randomSeed = rand() % 100000;
				++meta.fpNum;
				Parameters<PRISMATIC_FLOAT_PRECISION> prismatic_pars(meta);
				cout << "Frozen Phonon #" << fp_num << endl;
				prismatic_pars.meta.toString();
				PRISM01_calcPotential(prismatic_pars);
				Multislice_calcOutput(prismatic_pars);
				net_output += prismatic_pars.output;
			}
			// divide to take average
			for (auto&i:net_output) i/=prismatic_pars.meta.numFP;
			prismatic_pars.output = net_output;
		}

		if (prismatic_pars.meta.save3DOutput)prismatic_pars.output.toMRC_f(prismatic_pars.meta.filenameOutput.c_str());

		if (prismatic_pars.meta.save2DOutput) {
			size_t lower = std::max((size_t)0, (size_t)(prismatic_pars.meta.integrationAngleMin / prismatic_pars.meta.detectorAngleStep));
			size_t upper = std::min(prismatic_pars.detectorAngles.size(), (size_t) (prismatic_pars.meta.integrationAngleMax / prismatic_pars.meta.detectorAngleStep));
			Array2D<PRISMATIC_FLOAT_PRECISION> prism_image;
			prism_image = zeros_ND<2, PRISMATIC_FLOAT_PRECISION>(
					{{prismatic_pars.output.get_dimk(), prismatic_pars.output.get_dimj()}});
			for (auto y = 0; y < prismatic_pars.output.get_dimk(); ++y) {
				for (auto x = 0; x < prismatic_pars.output.get_dimj(); ++x) {
					for (auto b = lower; b < upper; ++b) {
						prism_image.at(y, x) += prismatic_pars.output.at(y, x, b);
					}
				}
			}
			std::string image_filename = std::string("multislice_2Doutput_") + prismatic_pars.meta.filenameOutput;
			prism_image.toMRC_f(image_filename.c_str());
		}
#ifdef PRISMATIC_ENABLE_GPU
		cout << "peak GPU memory usage = " << prismatic_pars.maxGPUMem << '\n';
#endif //PRISMATIC_ENABLE_GPU
		std::cout << "Calculation complete.\n" << std::endl;
		return prismatic_pars;
	}
}