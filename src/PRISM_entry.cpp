// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

#ifndef PRISM_PRISM_ENTRY_H
#define PRISM_PRISM_ENTRY_H
#include "PRISM_entry.h"
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include "configure.h"
#include "ArrayND.h"
#include "PRISM01.h"
#include "PRISM02.h"
#include "PRISM03.h"
#include "params.h"
#include <vector>

namespace PRISM{
	using namespace std;
	Parameters<PRISM_FLOAT_PRECISION> PRISM_entry(Metadata<PRISM_FLOAT_PRECISION>& meta){
		Parameters<PRISM_FLOAT_PRECISION> prism_pars(meta);
		PRISM01(prism_pars);
		prism_pars.pot.toMRC_f("debug_potential.mrc");
		PRISM02(prism_pars);

		Array3D<PRISM_FLOAT_PRECISION> tmp = zeros_ND<3, PRISM_FLOAT_PRECISION>({{prism_pars.Scompact.get_dimk(),prism_pars.Scompact.get_dimj(),prism_pars.Scompact.get_dimi()}});
		Array3D<PRISM_FLOAT_PRECISION> tmp_r = zeros_ND<3, PRISM_FLOAT_PRECISION>({{prism_pars.Scompact.get_dimk(),prism_pars.Scompact.get_dimj(),prism_pars.Scompact.get_dimi()}});
		Array3D<PRISM_FLOAT_PRECISION> tmp_i = zeros_ND<3, PRISM_FLOAT_PRECISION>({{prism_pars.Scompact.get_dimk(),prism_pars.Scompact.get_dimj(),prism_pars.Scompact.get_dimi()}});
		auto tmp_ptr = tmp.begin();
		auto tmp_r_ptr = tmp_r.begin();
		auto tmp_i_ptr = tmp_i.begin();
		for (auto&i : prism_pars.Scompact)*tmp_ptr++ = abs(i);
		for (auto&i : prism_pars.Scompact)*tmp_r_ptr++ = i.real();
		for (auto&i : prism_pars.Scompact)*tmp_i_ptr++ = i.imag();
		std::complex<PRISM_FLOAT_PRECISION> ssum{0,0};
		for (auto& i:prism_pars.Scompact)ssum+=i;
		cout <<"S compact sum = " << ssum << endl;
		tmp.toMRC_f("debug_scompact.mrc");
		tmp_r.toMRC_f("debug_scompact_r.mrc");
		tmp_i.toMRC_f("debug_scompact_i.mrc");

		PRISM03(prism_pars);

        if (prism_pars.meta.numFP > 1) {
            // run the rest of the frozen phonons
            Array3D<PRISM_FLOAT_PRECISION> net_output(prism_pars.output);
            for (auto fp_num = 1; fp_num < prism_pars.meta.numFP; ++fp_num){
				Parameters<PRISM_FLOAT_PRECISION> prism_pars(meta);
                cout << "Frozen Phonon #" << fp_num << endl;
//                ++prism_pars.meta.fpNum;
	        	meta.random_seed = rand() % 1000;
				prism_pars.meta = meta;
	        	PRISM01(prism_pars);
	        	PRISM02(prism_pars);
                PRISM03(prism_pars);
                net_output += prism_pars.output;
            }
            // divide to take average
            for (auto&i:net_output) i/=prism_pars.meta.numFP;
	    prism_pars.output = net_output;
        }
	        if (prism_pars.meta.save3DOutput)prism_pars.output.toMRC_f(prism_pars.meta.filename_output.c_str());

		if (prism_pars.meta.save2DOutput) {
//			size_t lower = 0;
//			size_t upper = 1;
			size_t lower = std::max((size_t)0, (size_t)(prism_pars.meta.integration_angle_min / prism_pars.meta.detector_angle_step));
			size_t upper = std::min((size_t)prism_pars.detectorAngles.size(), (size_t)(prism_pars.meta.integration_angle_max / prism_pars.meta.detector_angle_step));
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
			std::string image_filename = std::string("prism_image_") + prism_pars.meta.filename_output;
			prism_image.toMRC_f(image_filename.c_str());
		}

#ifdef PRISM_ENABLE_GPU
		cout << "peak GPU memory usage = " << prism_pars.max_mem << '\n';
#endif //PRISM_ENABLE_GPU
        std::cout << "PRISM Calculation complete.\n" << std::endl;
		return prism_pars;
	}
}
#endif //PRISM_PRISM_ENTRY_H
