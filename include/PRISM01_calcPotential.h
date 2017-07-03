// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// Prismatic is distributed under the GNU General Public License (GPL)
// If you use Prismatic, we kindly ask that you cite the following papers:

// 1. Ophus, C.: A fast image simulation algorithm for scanning
//    transmission electron microscopy. Advanced Structural and
//    Chemical Imaging 3(1), 13 (2017)

// 2. Pryor, Jr., A., Ophus, C., and Miao, J.: A Streaming Multi-GPU
//    Implementation of Image Simulation Algorithms for Scanning
//	  Transmission Electron Microscopy. arXiv:1706.08563 (2017)

#ifndef PRISMATIC_PRISM01_H
#define PRISMATIC_PRISM01_H
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


 namespace Prismatic {
	 void fetch_potentials(Array3D<PRISMATIC_FLOAT_PRECISION>& potentials,
	                      const std::vector<size_t>& atomic_species,
	                      const Array1D<PRISMATIC_FLOAT_PRECISION>& xr,
	                      const Array1D<PRISMATIC_FLOAT_PRECISION>& yr);

	 std::vector<size_t> get_unique_atomic_species(Parameters<PRISMATIC_FLOAT_PRECISION>& pars);

	 void generateProjectedPotentials(Parameters<PRISMATIC_FLOAT_PRECISION>& pars,
	                                 const Array3D<PRISMATIC_FLOAT_PRECISION>& potLookup,
	                                 const std::vector<size_t>& unique_species,
	                                 const Array1D<long>& xvec,
	                                 const Array1D<long>& yvec,
	                                 const Array1D<PRISMATIC_FLOAT_PRECISION>& uLookup);


//#ifdef PRISMATIC_BUILDING_GUI
//	void PRISM01_calcPotential(Parameters<PRISMATIC_FLOAT_PRECISION>& pars, prism_progressbar *progressbar=NULL);
//#else
	void PRISM01_calcPotential(Parameters<PRISMATIC_FLOAT_PRECISION>& pars);
//#endif //PRISMATIC_ENABLE_GPU
}
#endif //PRISMATIC_PRISM01_H
