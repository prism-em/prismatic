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
	                                const PRISMATIC_CUDA_COMPLEX_FLOAT *permutedScompact_d,
	                                const PRISMATIC_CUDA_COMPLEX_FLOAT *PsiProbeInit_d,
	                                const PRISMATIC_FLOAT_PRECISION *qxaReduce_d,
	                                const PRISMATIC_FLOAT_PRECISION *qyaReduce_d,
	                                const size_t *yBeams_d,
	                                const size_t *xBeams_d,
	                                const PRISMATIC_FLOAT_PRECISION *alphaInd_d,
	                                PRISMATIC_CUDA_COMPLEX_FLOAT *psi_ds,
	                                PRISMATIC_CUDA_COMPLEX_FLOAT *phaseCoeffs_ds,
	                                PRISMATIC_FLOAT_PRECISION *psiIntensity_ds,
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
	                               PRISMATIC_CUDA_COMPLEX_FLOAT *permutedScompact_ds,
	                               const std::complex<PRISMATIC_FLOAT_PRECISION> *permutedScompact_ph,
	                               const PRISMATIC_CUDA_COMPLEX_FLOAT *PsiProbeInit_d,
	                               const PRISMATIC_FLOAT_PRECISION *qxaReduce_d,
	                               const PRISMATIC_FLOAT_PRECISION *qyaReduce_d,
	                               const size_t *yBeams_d,
	                               const size_t *xBeams_d,
	                               const PRISMATIC_FLOAT_PRECISION *alphaInd_d,
	                               PRISMATIC_CUDA_COMPLEX_FLOAT *psi_ds,
	                               PRISMATIC_CUDA_COMPLEX_FLOAT *phaseCoeffs_ds,
	                               PRISMATIC_FLOAT_PRECISION *psiIntensity_ds,
	                               long *y_ds,
	                               long *x_ds,
	                               PRISMATIC_FLOAT_PRECISION *output_ph,
	                               PRISMATIC_FLOAT_PRECISION *integratedOutput_ds,
	                               const cufftHandle &cufft_plan,
	                               const cudaStream_t& stream,
	                               CudaParameters<PRISMATIC_FLOAT_PRECISION>& cuda_pars);
}
#endif //PRISMATIC_PRISM03_CUH