#ifndef PRISM_PRISM03_CUH
#define PRISM_PRISM03_CUH
#include "configure.h"
#include "cufft.h"
namespace PRISM {
	void buildPRISMOutput_GPU(Parameters<PRISM_FLOAT_PRECISION>& pars,
                              const PRISM_FLOAT_PRECISION xTiltShift,
                              const PRISM_FLOAT_PRECISION yTiltShift,
                              const Array2D<PRISM_FLOAT_PRECISION>& alphaInd,
                              const Array2D<std::complex<PRISM_FLOAT_PRECISION> >& PsiProbeInit);

	void buildSignal_GPU(Parameters<PRISM_FLOAT_PRECISION>&  pars,
	                     const size_t& ay,
	                     const size_t& ax,
	                     const PRISM_FLOAT_PRECISION& yTiltShift,
	                     const PRISM_FLOAT_PRECISION& xTiltShift,
	                     const PRISM_CUDA_COMPLEX_FLOAT *permuted_Scompact_d,
	                     const PRISM_CUDA_COMPLEX_FLOAT *PsiProbeInit_d,
	                     const PRISM_FLOAT_PRECISION *qxaReduce_d,
	                     const PRISM_FLOAT_PRECISION *qyaReduce_d,
	                     const size_t *yBeams_d,
	                     const size_t *xBeams_d,
	                     const PRISM_FLOAT_PRECISION *alphaInd_d,
	                     PRISM_CUDA_COMPLEX_FLOAT *psi_ds,
	                     PRISM_CUDA_COMPLEX_FLOAT *phaseCoeffs_ds,
	                     PRISM_FLOAT_PRECISION *psi_intensity_ds,
	                     size_t *y_ds,
	                     size_t *x_ds,
	                     PRISM_FLOAT_PRECISION *output_ph,
	                     const cufftHandle &cufft_plan,
	                     const cudaStream_t& stream);
}
#endif //PRISM_PRISM03_CUH