// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

#ifndef PRISMATIC_PRISM03_CUH
#define PRISMATIC_PRISM03_CUH
#include "configure.h"
#include "cufft.h"
#include "params.cuh"
namespace Prismatic {
	void buildPRISMOutput_GPU_singlexfer(Parameters<PRISMATIC_FLOAT_PRECISION>& pars);

	void buildPRISMOutput_GPU_streaming(Parameters<PRISMATIC_FLOAT_PRECISION>& pars);

	void buildSignal_GPU_singlexfer(Parameters<PRISMATIC_FLOAT_PRECISION>&  pars,
	                                const size_t& ay,
	                                const size_t& ax,
	                                const PRISMATIC_CUDA_COMPLEX_FLOAT *permuted_Scompact_d,
	                                const PRISMATIC_CUDA_COMPLEX_FLOAT *PsiProbeInit_d,
	                                const PRISMATIC_FLOAT_PRECISION *qxaReduce_d,
	                                const PRISMATIC_FLOAT_PRECISION *qyaReduce_d,
	                                const size_t *yBeams_d,
	                                const size_t *xBeams_d,
	                                const PRISMATIC_FLOAT_PRECISION *alphaInd_d,
	                                PRISMATIC_CUDA_COMPLEX_FLOAT *psi_ds,
	                                PRISMATIC_CUDA_COMPLEX_FLOAT *phaseCoeffs_ds,
	                                PRISMATIC_FLOAT_PRECISION *psi_intensity_ds,
	                                long  *y_ds,
	                                long  *x_ds,
	                                PRISMATIC_FLOAT_PRECISION *output_ph,
	                                PRISMATIC_FLOAT_PRECISION *integratedOutput_ds,
	                                const cufftHandle &cufft_plan,
	                                const cudaStream_t& stream,
	                                CudaParameters<PRISMATIC_FLOAT_PRECISION>& cuda_pars);

	void buildSignal_GPU_streaming(Parameters<PRISMATIC_FLOAT_PRECISION>&  pars,
	                               const size_t& ay,
	                               const size_t& ax,
	                               PRISMATIC_CUDA_COMPLEX_FLOAT *permuted_Scompact_ds,
	                               const std::complex<PRISMATIC_FLOAT_PRECISION> *permuted_Scompact_ph,
	                               const PRISMATIC_CUDA_COMPLEX_FLOAT *PsiProbeInit_d,
	                               const PRISMATIC_FLOAT_PRECISION *qxaReduce_d,
	                               const PRISMATIC_FLOAT_PRECISION *qyaReduce_d,
	                               const size_t *yBeams_d,
	                               const size_t *xBeams_d,
	                               const PRISMATIC_FLOAT_PRECISION *alphaInd_d,
	                               PRISMATIC_CUDA_COMPLEX_FLOAT *psi_ds,
	                               PRISMATIC_CUDA_COMPLEX_FLOAT *phaseCoeffs_ds,
	                               PRISMATIC_FLOAT_PRECISION *psi_intensity_ds,
	                               long *y_ds,
	                               long *x_ds,
	                               PRISMATIC_FLOAT_PRECISION *output_ph,
	                               PRISMATIC_FLOAT_PRECISION *integratedOutput_ds,
	                               const cufftHandle &cufft_plan,
	                               const cudaStream_t& stream,
	                               CudaParameters<PRISMATIC_FLOAT_PRECISION>& cuda_pars);
}
#endif //PRISMATIC_PRISM03_CUH