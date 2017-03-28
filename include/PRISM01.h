//
// Created by AJ Pryor on 2/21/17.
//

#ifndef PRISM_PRISM01_H
#define PRISM_PRISM01_H
#include "defines.h"
#include "params.h"
#include <iostream>
#include <algorithm>
#include <numeric>
#include <cstring>
#include <vector>
#include <map>
#include <random>
#include <thread>
#include "params.h"
#include "ArrayND.h"
#include "projectedPotential.h"

//#ifdef PRISM_BUILDING_GUI
//#include "prism_progressbar.h"
//#endif

 namespace PRISM {
	 void fetch_potentials(Array3D<PRISM_FLOAT_PRECISION>& potentials,
	                      const std::vector<size_t>& atomic_species,
	                      const Array1D<PRISM_FLOAT_PRECISION>& xr,
	                      const Array1D<PRISM_FLOAT_PRECISION>& yr);

	 std::vector<size_t> get_unique_atomic_species(Parameters<PRISM_FLOAT_PRECISION>& pars);

	 void generateProjectedPotentials(Parameters<PRISM_FLOAT_PRECISION>& pars,
	                                 const Array3D<PRISM_FLOAT_PRECISION>& potLookup,
	                                 const std::vector<size_t>& unique_species,
	                                 const Array1D<long>& xvec,
	                                 const Array1D<long>& yvec,
	                                 const Array1D<PRISM_FLOAT_PRECISION>& uLookup);


#ifdef PRISM_BUILDING_GUI
	void PRISM01(Parameters<PRISM_FLOAT_PRECISION>& pars, prism_progressbar *progressbar=NULL);
#else
	void PRISM01(Parameters<PRISM_FLOAT_PRECISION>& pars);
#endif //PRISM_ENABLE_GPU
}
#endif //PRISM_PRISM01_H
