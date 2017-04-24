// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

#include "meta.h"
#include "params.h"
#include "ArrayND.h"
#include "configure.h"
#include "Multislice.h"
#include "PRISM01.h"
#include "PRISM02.h"
#include <algorithm>


namespace PRISM{
	Parameters<PRISM_FLOAT_PRECISION> Multislice_entry(Metadata<PRISM_FLOAT_PRECISION>& meta){
		Parameters<PRISM_FLOAT_PRECISION> prism_pars(meta);
		PRISM01(prism_pars);
		Multislice(prism_pars);

		if (prism_pars.meta.numFP == 1) {
			prism_pars.output.toMRC_f(prism_pars.meta.filename_output.c_str());
		} else {
			// run the rest of the frozen phonons
			++prism_pars.meta.fpNum;
			Array3D<PRISM_FLOAT_PRECISION> net_output(prism_pars.output);
			for (auto fp_num = 1; fp_num < prism_pars.meta.numFP; ++fp_num){
				cout << "Frozen Phonon #" << fp_num << endl;
				Multislice(prism_pars);
				net_output += prism_pars.output;
			}
			// divide to take average
			for (auto&i:net_output) i/=prism_pars.meta.numFP;
			net_output.toMRC_f(prism_pars.meta.filename_output.c_str());
		}
#ifdef PRISM_ENABLE_GPU
		cout << "peak GPU memory usage = " << prism_pars.max_mem << '\n';
#endif //PRISM_ENABLE_GPU
		std::cout << "Calculation complete.\n" << std::endl;
		return prism_pars;
	}
}
