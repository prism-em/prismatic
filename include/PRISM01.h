//
// Created by AJ Pryor on 2/21/17.
//

#ifndef PRISM_PRISM01_H
#define PRISM_PRISM01_H
#include "params.h"
#include <iostream>
#include <algorithm>
#include <numeric>
#include <cstring>
#include <map>
#include <random>
#include <thread>
#include "params.h"
#include "ArrayND.h"
#include "projPot.h"
namespace PRISM {

	 void fetch_potentials(Array3D<PRISM_FLOAT_PRECISION>& potentials,
	                      const vector<size_t>& atomic_species,
	                      const Array1D<PRISM_FLOAT_PRECISION>& xr,
	                      const Array1D<PRISM_FLOAT_PRECISION>& yr);

	 vector<size_t> get_unique_atomic_species(Parameters<PRISM_FLOAT_PRECISION>& pars);

	 void generateProjectedPotentials(Parameters<PRISM_FLOAT_PRECISION>& pars,
	                                 const Array3D<PRISM_FLOAT_PRECISION>& potLookup,
	                                 const vector<size_t>& unique_species,
	                                 const Array1D<long>& xvec,
	                                 const Array1D<long>& yvec,
	                                 const Array1D<PRISM_FLOAT_PRECISION>& uLookup);

	 void PRISM01(Parameters<PRISM_FLOAT_PRECISION>& pars);
}
#endif //PRISM_PRISM01_H
