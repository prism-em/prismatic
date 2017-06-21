// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

#ifndef PRISMATIC_PRISM02_H
#define PRISMATIC_PRISM02_H
#include <mutex>
#include <complex>
#include "params.h"
#include "fftw3.h"
#include "configure.h"
#include "defines.h"

namespace Prismatic {
	inline void setupCoordinates(Parameters<PRISMATIC_FLOAT_PRECISION>& pars);
	inline void setupBeams(Parameters<PRISMATIC_FLOAT_PRECISION>& pars);
	inline void setupSMatrixCoordinates(Parameters<PRISMATIC_FLOAT_PRECISION>& pars);
	inline void downsampleFourierComponents(Parameters<PRISMATIC_FLOAT_PRECISION> &pars);

	void propagatePlaneWave_CPU(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
	                            size_t currentBeam,
	                            Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> > &psi,
	                            const PRISMATIC_FFTW_PLAN &plan_forward,
	                            const PRISMATIC_FFTW_PLAN &plan_inverse,
	                            std::mutex& fftw_plan_lock);

	void propagatePlaneWave_CPU_batch(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
	                                  size_t currentBeam,
	                                  size_t stopBeam,
	                                  Array1D<std::complex<PRISMATIC_FLOAT_PRECISION> > &psi_stack,
	                                  const PRISMATIC_FFTW_PLAN &plan_forward,
	                                  const PRISMATIC_FFTW_PLAN &plan_inverse,
	                                  std::mutex &fftw_plan_lock);

	void fill_Scompact_CPUOnly(Parameters<PRISMATIC_FLOAT_PRECISION> &pars);

	void PRISM02_calcSMatrix(Parameters<PRISMATIC_FLOAT_PRECISION>& pars);

}
#endif //PRISMATIC_PRISM02_H
