//
// Created by AJ Pryor on 3/13/17.
//

//
// Created by AJ Pryor on 3/2/17.
//

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
		cout << prism_pars.pot.at(0,0,0);
//		prism_pars.pot.toMRC_f("test.mrc");
		PRISM02(prism_pars);

//		Array3D<PRISM_FLOAT_PRECISION> tmp = zeros_ND<3, PRISM_FLOAT_PRECISION>({{prism_pars.Scompact.get_dimk(),prism_pars.Scompact.get_dimj(),prism_pars.Scompact.get_dimi()}});
//		auto tmp_ptr = tmp.begin();
//		for (auto&i : prism_pars.Scompact)*tmp_ptr++ = abs(i);
//		tmp.toMRC_f("debug_scompact.mrc");

		PRISM03(prism_pars);

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
		std::cout << "PRISM Calculation complete.\n" << std::endl;
//		return 0;
		return prism_pars;
	}

}
#endif //PRISM_PRISM_ENTRY_H
