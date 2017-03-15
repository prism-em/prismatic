#ifndef PRISM_PRISM03_CUH
#define PRISM_PRISM03_CUH
#include "configure.h"
namespace PRISM {
	void buildPRISMOutput_GPU(Parameters<PRISM_FLOAT_PRECISION>& pars,
                              const PRISM_FLOAT_PRECISION xTiltShift,
                              const PRISM_FLOAT_PRECISION yTiltShift,
                              const Array2D<PRISM_FLOAT_PRECISION>& alphaInd,
                              const Array2D<std::complex<PRISM_FLOAT_PRECISION> >& PsiProbeInit);
	void buildSignal_GPU(Parameters<PRISM_FLOAT_PRECISION>&  pars,const size_t& ay,
	                     const size_t& ax,
	                     const PRISM_FLOAT_PRECISION& yTiltShift,
	                     const PRISM_FLOAT_PRECISION& xTiltShift,
	                     const Array2D<PRISM_FLOAT_PRECISION>& alphaInd,
	                     const Array2D<std::complex<PRISM_FLOAT_PRECISION> >& PsiProbeInit);
}
#endif //PRISM_PRISM03_CUH