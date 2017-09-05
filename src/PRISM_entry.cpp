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


namespace Prismatic{
	using namespace std;
	Parameters<PRISMATIC_FLOAT_PRECISION> PRISM_entry(Metadata<PRISMATIC_FLOAT_PRECISION>& meta){
		Parameters<PRISMATIC_FLOAT_PRECISION> prismatic_pars;
		try { // read atomic coordinates
			prismatic_pars = Parameters<PRISMATIC_FLOAT_PRECISION>(meta);
		} catch(...){
			std::cout << "Terminating" << std::endl;
			exit(1);
		}
		prismatic_pars.meta.toString();

//        to_xyz(prismatic_pars.atoms, "/Users/ajpryor/Documents/MATLAB/multislice/PRISM/build/test.XYZ", "comment", 5.43,5.43,5.43);

		// compute projected potentials
		PRISM01_calcPotential(prismatic_pars);

//		prismatic_pars.pot.toMRC_f("debug_potential.mrc");

		// compute compact S-matrix
		PRISM02_calcSMatrix(prismatic_pars);

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
	        	PRISM02_calcSMatrix(prismatic_pars);
                PRISM03_calcOutput(prismatic_pars);
                net_output += prismatic_pars.output;
            }
            // divide to take average
            for (auto&i:net_output) i/=prismatic_pars.meta.numFP;
	        prismatic_pars.output = net_output;
        }
        if (prismatic_pars.meta.save3DOutput)prismatic_pars.output.toMRC_f(prismatic_pars.meta.filenameOutput.c_str());

		if (prismatic_pars.meta.save2DOutput) {
			size_t lower = std::max((size_t)0, (size_t)(prismatic_pars.meta.integrationAngleMin / prismatic_pars.meta.detectorAngleStep));
			size_t upper = std::min((size_t)prismatic_pars.detectorAngles.size(), (size_t)(prismatic_pars.meta.integrationAngleMax / prismatic_pars.meta.detectorAngleStep));
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
			std::string image_filename = std::string("prism_2Doutput_") + prismatic_pars.meta.filenameOutput;
			prism_image.toMRC_f(image_filename.c_str());
		}

#ifdef PRISMATIC_ENABLE_GPU
		cout << "peak GPU memory usage = " << prismatic_pars.maxGPUMem << '\n';
#endif //PRISMATIC_ENABLE_GPU
        std::cout << "PRISM Calculation complete.\n" << std::endl;
		return prismatic_pars;
	}
}
