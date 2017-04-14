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
		Multislice(prism_pars);

		// for now stack is a 4D array so here I convert to 3D to save it
		PRISM::Array3D<PRISM_FLOAT_PRECISION> reshaped_output = PRISM::zeros_ND<3, PRISM_FLOAT_PRECISION>(
				{{prism_pars.stack.get_diml(), prism_pars.stack.get_dimk(), prism_pars.stack.get_dimj()}});
		auto ptr = reshaped_output.begin();
		for (auto &i:prism_pars.stack)*ptr++=i;
		reshaped_output.toMRC_f(prism_pars.meta.filename_output.c_str());

		std::cout << "Calculation complete.\n" << std::endl;
		return prism_pars;
	}
}
