// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

#ifndef PRISM_PRISM03_H
#define PRISM_PRISM03_H
#include "params.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <thread>
#include <mutex>
#include <numeric>
#include "fftw3.h"
#include "utility.h"

namespace PRISM {
	Array2D<PRISM_FLOAT_PRECISION> array2D_subset(const Array2D<PRISM_FLOAT_PRECISION> &arr,
	                          const size_t &starty, const size_t &stepy, const size_t &stopy,
	                          const size_t &startx, const size_t &stepx, const size_t &stopx);
	inline void setupCoordinates(Parameters<PRISM_FLOAT_PRECISION> &pars);
	inline void setupDetector(Parameters<PRISM_FLOAT_PRECISION> &pars);
	inline void setupBeams(Parameters<PRISM_FLOAT_PRECISION> &pars);
	inline void createStack_integrate(Parameters<PRISM_FLOAT_PRECISION> &pars);
	inline void setupFourierCoordinates(Parameters<PRISM_FLOAT_PRECISION> &pars);

	void buildSignal_CPU(Parameters<PRISM_FLOAT_PRECISION> &pars,
						 const size_t &ay,
						 const size_t &ax,
						 PRISM_FFTW_PLAN& plan,
						 Array2D<std::complex<PRISM_FLOAT_PRECISION> >& psi);

	void buildPRISMOutput_CPUOnly(Parameters<PRISM_FLOAT_PRECISION>& pars);
	void PRISM03(Parameters<PRISM_FLOAT_PRECISION> &pars);
}
#endif //PRISM_PRISM03_H
