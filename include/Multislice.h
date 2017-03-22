//
// Created by AJ Pryor on 3/6/17.
//

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
#include "getWorkID.h"
namespace PRISM{
	using namespace std;
	void setupCoordinates_multislice(Parameters<PRISM_FLOAT_PRECISION>& pars);

	void setupDetector_multislice(Parameters<PRISM_FLOAT_PRECISION>& pars);

	void setupProbes_multislice(Parameters<PRISM_FLOAT_PRECISION>& pars);

	void formatOutput_CPU_integrate(Parameters<PRISM_FLOAT_PRECISION>& pars,
	                             Array2D< complex<PRISM_FLOAT_PRECISION> >& psi,
	                             const Array2D<PRISM_FLOAT_PRECISION> &alphaInd,
	                             const size_t& ay,
	                             const size_t& ax);

	void getMultisliceProbe_CPU(Parameters<PRISM_FLOAT_PRECISION>& pars,
                                Array3D<complex<PRISM_FLOAT_PRECISION> >& trans,
                                const Array2D<complex<PRISM_FLOAT_PRECISION> >& PsiProbeInit,
                                const size_t& ay,
                                const size_t& ax,
                                const Array2D<PRISM_FLOAT_PRECISION> &alphaInd);

	void buildMultisliceOutput_CPUOnly(Parameters<PRISM_FLOAT_PRECISION>& pars,
                                       Array3D<complex<PRISM_FLOAT_PRECISION> >& trans,
                                       Array2D<complex<PRISM_FLOAT_PRECISION> >& PsiProbeInit,
                                       Array2D<PRISM_FLOAT_PRECISION> &alphaInd);


	void Multislice(Parameters<PRISM_FLOAT_PRECISION>& pars);
}
#endif //PRISM_MULTISLICE_H
