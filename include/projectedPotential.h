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

#ifndef PRISM_PROJPOT_H
#define PRISM_PROJPOT_H
#include <vector>
#include "ArrayND.h"
#include <math.h>
#include <iostream>
#include <algorithm>
#include <numeric>
#include "kirkland_params.h"
#include "configure.h"
namespace Prismatic {

	PRISMATIC_FLOAT_PRECISION get_potMin(const Array2D<PRISMATIC_FLOAT_PRECISION>& pot,
                                     const Array1D<PRISMATIC_FLOAT_PRECISION>& xr,
                                     const Array1D<PRISMATIC_FLOAT_PRECISION>& yr);

	Array2D<PRISMATIC_FLOAT_PRECISION> projPot(const size_t &Z,
	                                       const Array1D<PRISMATIC_FLOAT_PRECISION> &xr,
	                                       const Array1D<PRISMATIC_FLOAT_PRECISION> &yr);

}
#endif //PRISM_PROJPOT_H
