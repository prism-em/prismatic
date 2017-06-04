// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

#include "meta.h"
#include "params.h"
#include "ArrayND.h"
#include "configure.h"
#include "Multislice_calcOutput.h"
#include "PRISM01_calcPotential.h"
#include "PRISM02_calcSMatrix.h"
#include <algorithm>


namespace PRISM{
	Parameters<PRISM_FLOAT_PRECISION> Multislice_entry(Metadata<PRISM_FLOAT_PRECISION>& meta){
		Parameters<PRISM_FLOAT_PRECISION> prism_pars;
		try { // read atomic coordinates
			prism_pars = Parameters<PRISM_FLOAT_PRECISION>(meta);
		} catch(const std::runtime_error &e){
			std::cout << "Terminating" << std::endl;
			exit(1);
		}

		// compute projected potentials
		PRISM01_calcPotential(prism_pars);

		// compute final output
		Multislice_calcOutput(prism_pars);

		// calculate remaining frozen phonon configurations
		if (prism_pars.meta.numFP > 1) {
			// run the rest of the frozen phonons
			Array3D<PRISM_FLOAT_PRECISION> net_output(prism_pars.output);
			for (auto fp_num = 1; fp_num < prism_pars.meta.numFP; ++fp_num){
				meta.random_seed = rand() % 100000;
				Parameters<PRISM_FLOAT_PRECISION> prism_pars(meta);
				cout << "Frozen Phonon #" << fp_num << endl;
				prism_pars.meta.toString();
				PRISM01_calcPotential(prism_pars);
				Multislice_calcOutput(prism_pars);
				net_output += prism_pars.output;
			}
			// divide to take average
			for (auto&i:net_output) i/=prism_pars.meta.numFP;
			prism_pars.output = net_output;
		}

		if (prism_pars.meta.save3DOutput)prism_pars.output.toMRC_f(prism_pars.meta.filename_output.c_str());

		if (prism_pars.meta.save2DOutput) {
			size_t lower = std::max((size_t)0, (size_t)(prism_pars.meta.integration_angle_min / prism_pars.meta.detector_angle_step));
			size_t upper = std::min(prism_pars.detectorAngles.size(), (size_t) (prism_pars.meta.integration_angle_max / prism_pars.meta.detector_angle_step));
			Array2D<PRISM_FLOAT_PRECISION> prism_image;
			prism_image = zeros_ND<2, PRISM_FLOAT_PRECISION>(
					{{prism_pars.output.get_dimk(), prism_pars.output.get_dimj()}});
			for (auto y = 0; y < prism_pars.output.get_dimk(); ++y) {
				for (auto x = 0; x < prism_pars.output.get_dimj(); ++x) {
					for (auto b = lower; b < upper; ++b) {
						prism_image.at(y, x) += prism_pars.output.at(y, x, b);
					}
				}
			}
			std::string image_filename = std::string("multislice_2Doutput_") + prism_pars.meta.filename_output;
			prism_image.toMRC_f(image_filename.c_str());
		}
#ifdef PRISM_ENABLE_GPU
		cout << "peak GPU memory usage = " << prism_pars.max_mem << '\n';
#endif //PRISM_ENABLE_GPU
		std::cout << "Calculation complete.\n" << std::endl;
		return prism_pars;
	}
}