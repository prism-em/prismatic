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
//		prism_pars.pot.toMRC_f("debug_potential.mrc");
		PRISM02(prism_pars);

//		Array3D<PRISM_FLOAT_PRECISION> tmp = zeros_ND<3, PRISM_FLOAT_PRECISION>({{prism_pars.Scompact.get_dimk(),prism_pars.Scompact.get_dimj(),prism_pars.Scompact.get_dimi()}});
//		auto tmp_ptr = tmp.begin();
//		for (auto&i : prism_pars.Scompact)*tmp_ptr++ = abs(i);
//		tmp.toMRC_f("debug_scompact.mrc");

		PRISM03(prism_pars);

//		size_t lower = 13;
//		size_t upper = 18;
//		Array2D<PRISM_FLOAT_PRECISION> prism_image;
//		prism_image = zeros_ND<2, PRISM_FLOAT_PRECISION>({{prism_pars.output.get_dimk(), prism_pars.output.get_dimj()}});
//		for (auto y = 0; y < prism_pars.output.get_dimk(); ++y){
//			for (auto x = 0; x < prism_pars.output.get_dimj(); ++x){
//				for (auto b = lower; b < upper; ++b){
//					prism_image.at(y,x) += prism_pars.output.at(y,x,b);
//				}
//			}
//		}
//		prism_image.toMRC_f("prism_image.mrc");
        if (prism_pars.meta.numFP == 1) {
            prism_pars.output.toMRC_f(prism_pars.meta.filename_output.c_str());
        } else {

            // run the rest of the frozen phonons
            Array3D<PRISM_FLOAT_PRECISION> net_output(prism_pars.output);
            for (auto fp_num = 1; fp_num < prism_pars.meta.numFP; ++fp_num){
                cout << "Frozen Phonon #" << fp_num << endl;
                PRISM03(prism_pars);
                net_output += prism_pars.output;
            }
            // divide to take average
            for (auto&i:net_output) i/=prism_pars.meta.numFP;
            net_output.toMRC_f(prism_pars.meta.filename_output.c_str());
        }
#ifdef PRISM_ENABLE_GPU
		cout << "peak GPU memory usage = " << prism_pars.max_mem << '\n';
#endif //PRISM_ENABLE_GPU
        std::cout << "PRISM Calculation complete.\n" << std::endl;
		return prism_pars;
	}
}
#endif //PRISM_PRISM_ENTRY_H
