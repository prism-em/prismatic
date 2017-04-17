// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

#ifndef PRISM_PROJPOT_H
#define PRISM_PROJPOT_H
#include <vector>
#include "ArrayND.h"
#include "boost/math/special_functions/bessel.hpp"
#include <math.h>
#include <iostream>
#include <algorithm>
#include <numeric>
#include "fparams.h"
#include "configure.h"
namespace PRISM {

	PRISM_FLOAT_PRECISION get_potMin(const Array2D<PRISM_FLOAT_PRECISION>& pot,
                                     const Array1D<PRISM_FLOAT_PRECISION>& xr,
                                     const Array1D<PRISM_FLOAT_PRECISION>& yr);

	Array2D<PRISM_FLOAT_PRECISION> projPot(const size_t &Z,
	                                       const Array1D<PRISM_FLOAT_PRECISION> &xr,
	                                       const Array1D<PRISM_FLOAT_PRECISION> &yr);

}
#endif //PRISM_PROJPOPRISM_FLOAT_PRECISION_H
