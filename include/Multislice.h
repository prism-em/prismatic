// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

#ifndef PRISM_MULTISLICE_H
#define PRISM_MULTISLICE_H
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

namespace PRISM{
	using namespace std;
	void setupCoordinates_multislice(Parameters<PRISM_FLOAT_PRECISION>& pars);

	void setupDetector_multislice(Parameters<PRISM_FLOAT_PRECISION>& pars);

	void setupProbes_multislice(Parameters<PRISM_FLOAT_PRECISION>& pars);

	void createTransmission(Parameters<PRISM_FLOAT_PRECISION>& pars);

	void createStack(Parameters<PRISM_FLOAT_PRECISION>& pars);

	void formatOutput_CPU_integrate(Parameters<PRISM_FLOAT_PRECISION>& pars,
	                             Array2D< complex<PRISM_FLOAT_PRECISION> >& psi,
	                             const Array2D<PRISM_FLOAT_PRECISION> &alphaInd,
	                             const size_t ay,
	                             const size_t ax);
	std::pair<PRISM::Array2D< std::complex<PRISM_FLOAT_PRECISION> >, PRISM::Array2D< std::complex<PRISM_FLOAT_PRECISION> > >
	getSingleMultisliceProbe_CPU(Parameters<PRISM_FLOAT_PRECISION> &pars, const PRISM_FLOAT_PRECISION xp, const PRISM_FLOAT_PRECISION yp);
	void getMultisliceProbe_CPU(Parameters<PRISM_FLOAT_PRECISION>& pars,
                                const size_t ay,
                                const size_t ax,
								PRISM_FFTW_PLAN& plan_forward,
								PRISM_FFTW_PLAN& plan_inverse,
								Array2D<complex<PRISM_FLOAT_PRECISION> >& psi);

	void buildMultisliceOutput_CPUOnly(Parameters<PRISM_FLOAT_PRECISION>& pars);


	void Multislice(Parameters<PRISM_FLOAT_PRECISION>& pars);
}
#endif //PRISM_MULTISLICE_H
