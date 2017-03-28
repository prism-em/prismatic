//
// Created by AJ Pryor on 3/5/17.
//

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
//		prism_pars.pot.toMRC_f("DEBUG.mrc");
		Multislice(prism_pars);
		size_t lower = 13;
		size_t upper = 18;
		Array2D<PRISM_FLOAT_PRECISION> prism_image;
		prism_image = zeros_ND<2, PRISM_FLOAT_PRECISION>({{prism_pars.stack.get_diml(), prism_pars.stack.get_dimk()}});

		for (auto y = 0; y < prism_pars.stack.get_diml(); ++y){
			for (auto x = 0; x < prism_pars.stack.get_dimk(); ++x){
				for (auto b = lower; b < upper; ++b){
					prism_image.at(y,x) += prism_pars.stack.at(y,x,b,0);
				}
			}
		}

		prism_image.toMRC_f(prism_pars.meta.filename_output.c_str());
		std::cout << "Calculation complete.\n" << std::endl;
		return prism_pars;
	}
}
