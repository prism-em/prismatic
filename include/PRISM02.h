//
// Created by AJ Pryor on 2/24/17.
//

#ifndef PRISM_PRISM02_H
#define PRISM_PRISM02_H
#include <mutex>
#include <complex>
#include "params.h"
#include "fftw3.h"
#include "configure.h"
#include "defines.h"

namespace PRISM {
	inline void setupCoordinates(Parameters<PRISM_FLOAT_PRECISION>& pars);
	inline void setupBeams(Parameters<PRISM_FLOAT_PRECISION>& pars);
	inline void setupSMatrixCoordinates(Parameters<PRISM_FLOAT_PRECISION>& pars);
	inline void downsampleFourierComponents(Parameters<PRISM_FLOAT_PRECISION> &pars);
	void propagatePlaneWave_CPU(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                            size_t a0,
	                            Array2D<std::complex<PRISM_FLOAT_PRECISION> > &psi,
	                            const PRISM_FFTW_PLAN &plan_forward,
	                            const PRISM_FFTW_PLAN &plan_inverse,
	                            std::mutex& fftw_plan_lock);

	void fill_Scompact_CPUOnly(Parameters<PRISM_FLOAT_PRECISION> &pars);

	void PRISM02(Parameters<PRISM_FLOAT_PRECISION>& pars);

}
#endif //PRISM_PRISM02_H
