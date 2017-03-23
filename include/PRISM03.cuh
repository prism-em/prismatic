#ifndef PRISM_PRISM03_CUH
#define PRISM_PRISM03_CUH
#include "configure.h"
#include "cufft.h"
namespace PRISM {
	void buildPRISMOutput_GPU_singlexfer(Parameters<PRISM_FLOAT_PRECISION>& pars);

	void buildPRISMOutput_GPU_streaming(Parameters<PRISM_FLOAT_PRECISION>& pars);

	void buildSignal_GPU_singlexfer(Parameters<PRISM_FLOAT_PRECISION>&  pars,
	                                const size_t& ay,
	                                const size_t& ax,
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
	                                long *y_ds,
	                                long *x_ds,
	                                PRISM_FLOAT_PRECISION *output_ph,
	                                PRISM_FLOAT_PRECISION *integratedOutput_ds,
	                                const cufftHandle &cufft_plan,
	                                const cudaStream_t& stream);

	void buildSignal_GPU_streaming(Parameters<PRISM_FLOAT_PRECISION>&  pars,
	                               const size_t& ay,
	                               const size_t& ax,
	                               PRISM_CUDA_COMPLEX_FLOAT *permuted_Scompact_ds,
	                               const std::complex<PRISM_FLOAT_PRECISION> *permuted_Scompact_ph,
	                               const PRISM_CUDA_COMPLEX_FLOAT *PsiProbeInit_d,
	                               const PRISM_FLOAT_PRECISION *qxaReduce_d,
	                               const PRISM_FLOAT_PRECISION *qyaReduce_d,
	                               const size_t *yBeams_d,
	                               const size_t *xBeams_d,
	                               const PRISM_FLOAT_PRECISION *alphaInd_d,
	                               PRISM_CUDA_COMPLEX_FLOAT *psi_ds,
	                               PRISM_CUDA_COMPLEX_FLOAT *phaseCoeffs_ds,
	                               PRISM_FLOAT_PRECISION *psi_intensity_ds,
	                               long *y_ds,
	                               long *x_ds,
	                               PRISM_FLOAT_PRECISION *output_ph,
	                               PRISM_FLOAT_PRECISION *integratedOutput_ds,
	                               const cufftHandle &cufft_plan,
	                               const cudaStream_t& stream);
}
#endif //PRISM_PRISM03_CUH