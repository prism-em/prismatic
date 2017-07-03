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

#ifndef PRISMATIC_PRISM03_H
#define PRISMATIC_PRISM03_H
#include "params.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <thread>
#include <mutex>
#include <numeric>
#include "fftw3.h"
#include "utility.h"

namespace Prismatic {
	Array2D<PRISMATIC_FLOAT_PRECISION> array2D_subset(const Array2D<PRISMATIC_FLOAT_PRECISION> &arr,
	                          const size_t &starty, const size_t &stepy, const size_t &stopy,
	                          const size_t &startx, const size_t &stepx, const size_t &stopx);

	void setupCoordinates_2(Parameters<PRISMATIC_FLOAT_PRECISION> &pars);
	void setupDetector(Parameters<PRISMATIC_FLOAT_PRECISION> &pars);
	void setupBeams_2(Parameters<PRISMATIC_FLOAT_PRECISION> &pars);
	void createStack_integrate(Parameters<PRISMATIC_FLOAT_PRECISION> &pars);
	void setupFourierCoordinates(Parameters<PRISMATIC_FLOAT_PRECISION> &pars);
	void transformIndices(Parameters<PRISMATIC_FLOAT_PRECISION> &pars);
	void initializeProbes(Parameters<PRISMATIC_FLOAT_PRECISION> &pars);
	std::pair<Prismatic::Array2D< std::complex<PRISMATIC_FLOAT_PRECISION> >, Prismatic::Array2D< std::complex<PRISMATIC_FLOAT_PRECISION> > >
	getSinglePRISMProbe_CPU(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, const PRISMATIC_FLOAT_PRECISION xp, const PRISMATIC_FLOAT_PRECISION yp);
	void buildSignal_CPU(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
						 const size_t &ay,
						 const size_t &ax,
						 PRISMATIC_FFTW_PLAN& plan,
						 Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> >& psi);

	void buildPRISMOutput_CPUOnly(Parameters<PRISMATIC_FLOAT_PRECISION>& pars);

	void PRISM03_calcOutput(Parameters<PRISMATIC_FLOAT_PRECISION> &pars);
}
#endif //PRISMATIC_PRISM03_H
