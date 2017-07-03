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

#ifndef PRISMATIC_MULTISLICE_H
#define PRISMATIC_MULTISLICE_H
#include <iostream>
#include <thread>
#include <vector>
#include <mutex>
#include <numeric>
#include "configure.h"
#include <numeric>
#include "meta.h"
#include "ArrayND.h"
#include "params.h"
#include "utility.h"
#include "fftw3.h"
#include "WorkDispatcher.h"

namespace Prismatic{
	using namespace std;
	void setupCoordinates_multislice(Parameters<PRISMATIC_FLOAT_PRECISION>& pars);

	void setupDetector_multislice(Parameters<PRISMATIC_FLOAT_PRECISION>& pars);

	void setupProbes_multislice(Parameters<PRISMATIC_FLOAT_PRECISION>& pars);

	void createTransmission(Parameters<PRISMATIC_FLOAT_PRECISION>& pars);

	void createStack(Parameters<PRISMATIC_FLOAT_PRECISION>& pars);

	void formatOutput_CPU_integrate(Parameters<PRISMATIC_FLOAT_PRECISION>& pars,
	                                Array2D< complex<PRISMATIC_FLOAT_PRECISION> >& psi,
	                                const Array2D<PRISMATIC_FLOAT_PRECISION> &alphaInd,
	                                const size_t ay,
	                                const size_t ax);

	void formatOutput_CPU_integrate_batch(Parameters<PRISMATIC_FLOAT_PRECISION>& pars,
	                                      Array1D< complex<PRISMATIC_FLOAT_PRECISION> >& psi_stack,
	                                      const Array2D<PRISMATIC_FLOAT_PRECISION> &alphaInd,
	                                      const size_t ay,
	                                      const size_t ax);

	std::pair<Prismatic::Array2D< std::complex<PRISMATIC_FLOAT_PRECISION> >, Prismatic::Array2D< std::complex<PRISMATIC_FLOAT_PRECISION> > >
	getSingleMultisliceProbe_CPU(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, const PRISMATIC_FLOAT_PRECISION xp, const PRISMATIC_FLOAT_PRECISION yp);
	void getMultisliceProbe_CPU_batch(Parameters<PRISMATIC_FLOAT_PRECISION>& pars,
	                                  const size_t Nstart,
	                                  const size_t Nstop,
	                                  PRISMATIC_FFTW_PLAN& plan_forward,
	                                  PRISMATIC_FFTW_PLAN& plan_inverse,
	                                  Array1D<complex<PRISMATIC_FLOAT_PRECISION> >& psi_stack);
	void getMultisliceProbe_CPU(Parameters<PRISMATIC_FLOAT_PRECISION>& pars,
	                            const size_t ay,
	                            const size_t ax,
	                            PRISMATIC_FFTW_PLAN& plan_forward,
	                            PRISMATIC_FFTW_PLAN& plan_inverse,
	                            Array2D<complex<PRISMATIC_FLOAT_PRECISION> >& psi);
	void buildMultisliceOutput_CPUOnly(Parameters<PRISMATIC_FLOAT_PRECISION>& pars);


	void Multislice_calcOutput(Parameters<PRISMATIC_FLOAT_PRECISION>& pars);
}
#endif //PRISMATIC_MULTISLICE_H
