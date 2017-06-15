// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

#ifndef PRISM_PRISM03_CUH
#define PRISM_PRISM03_CUH
#include "configure.h"
#include "cufft.h"
#include "params.cuh"
namespace PRISM {
	void buildPRISMOutput_GPU_singlexfer(Parameters<PRISMATIC_FLOAT_PRECISION>& pars);

	void buildPRISMOutput_GPU_streaming(Parameters<PRISMATIC_FLOAT_PRECISION>& pars);

	void buildSignal_GPU_singlexfer(Parameters<PRISMATIC_FLOAT_PRECISION>&  pars,
	                                const size_t& ay,
	                                const size_t& ax,
	                                const PRISM_CUDA_COMPLEX_FLOAT *permuted_Scompact_d,
	                                const PRISM_CUDA_COMPLEX_FLOAT *PsiProbeInit_d,
	                                const PRISMATIC_FLOAT_PRECISION *qxaReduce_d,
	                                const PRISMATIC_FLOAT_PRECISION *qyaReduce_d,
	                                const size_t *yBeams_d,
	                                const size_t *xBeams_d,
	                                const PRISMATIC_FLOAT_PRECISION *alphaInd_d,
	                                PRISM_CUDA_COMPLEX_FLOAT *psi_ds,
	                                PRISM_CUDA_COMPLEX_FLOAT *phaseCoeffs_ds,
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
	                               PRISM_CUDA_COMPLEX_FLOAT *permuted_Scompact_ds,
	                               const std::complex<PRISMATIC_FLOAT_PRECISION> *permuted_Scompact_ph,
	                               const PRISM_CUDA_COMPLEX_FLOAT *PsiProbeInit_d,
	                               const PRISMATIC_FLOAT_PRECISION *qxaReduce_d,
	                               const PRISMATIC_FLOAT_PRECISION *qyaReduce_d,
	                               const size_t *yBeams_d,
	                               const size_t *xBeams_d,
	                               const PRISMATIC_FLOAT_PRECISION *alphaInd_d,
	                               PRISM_CUDA_COMPLEX_FLOAT *psi_ds,
	                               PRISM_CUDA_COMPLEX_FLOAT *phaseCoeffs_ds,
	                               PRISMATIC_FLOAT_PRECISION *psi_intensity_ds,
	                               long *y_ds,
	                               long *x_ds,
	                               PRISMATIC_FLOAT_PRECISION *output_ph,
	                               PRISMATIC_FLOAT_PRECISION *integratedOutput_ds,
	                               const cufftHandle &cufft_plan,
	                               const cudaStream_t& stream,
	                               CudaParameters<PRISMATIC_FLOAT_PRECISION>& cuda_pars);
}
#endif //PRISM_PRISM03_CUH