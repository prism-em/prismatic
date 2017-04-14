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
		prism_pars.stack.toMRC_f(prism_pars.meta.filename_output.c_str());
		std::cout << "Calculation complete.\n" << std::endl;
		return prism_pars;
	}
}
