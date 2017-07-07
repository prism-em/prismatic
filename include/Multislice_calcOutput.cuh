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

#ifndef PRISMATIC_MULTISLICE_CUH
#define PRISMATIC_MULTISLICE_CUH
#include <cuda_runtime.h>
#include "params.h"
#include "configure.h"
#include <complex>
#ifdef PRISMATIC_BUILDING_GUI
#include "prism_progressbar.h"
#endif
namespace Prismatic {

	 void buildMultisliceOutput_GPU_singlexfer(Parameters <PRISMATIC_FLOAT_PRECISION> &pars);

	 void buildMultisliceOutput_GPU_streaming(Parameters <PRISMATIC_FLOAT_PRECISION> &pars);

	 void getMultisliceProbe_GPU_singlexfer(Parameters<PRISMATIC_FLOAT_PRECISION>& pars,
	                                       PRISMATIC_CUDA_COMPLEX_FLOAT* trans_d,
	                                       PRISMATIC_CUDA_COMPLEX_FLOAT* PsiProbeInit_d,
	                                       PRISMATIC_CUDA_COMPLEX_FLOAT* psi_ds,
	                                       PRISMATIC_FLOAT_PRECISION* output_ph,
	                                       PRISMATIC_FLOAT_PRECISION* psiIntensity_ds,
	                                       PRISMATIC_FLOAT_PRECISION* integratedOutput_ds,
	                                       const PRISMATIC_FLOAT_PRECISION* qya_d,
	                                       const PRISMATIC_FLOAT_PRECISION* qxa_d,
	                                       const PRISMATIC_CUDA_COMPLEX_FLOAT* prop_d,
	                                       const size_t ay,
	                                       const size_t ax,
	                                       const size_t dimj,
	                                       const size_t dimi,
	                                       const PRISMATIC_FLOAT_PRECISION* alphaInd_d,
	                                       const cufftHandle& plan,
	                                       cudaStream_t& stream);

	 void getMultisliceProbe_GPU_singlexfer_batch(Parameters<PRISMATIC_FLOAT_PRECISION>& pars,
	                                             PRISMATIC_CUDA_COMPLEX_FLOAT* trans_d,
	                                             PRISMATIC_CUDA_COMPLEX_FLOAT* PsiProbeInit_d,
	                                             PRISMATIC_CUDA_COMPLEX_FLOAT* psi_ds,
	                                             PRISMATIC_FLOAT_PRECISION* output_ph,
	                                             PRISMATIC_FLOAT_PRECISION* psiIntensity_ds,
	                                             PRISMATIC_FLOAT_PRECISION* integratedOutput_ds,
	                                             const PRISMATIC_FLOAT_PRECISION* qya_d,
	                                             const PRISMATIC_FLOAT_PRECISION* qxa_d,
	                                             const PRISMATIC_CUDA_COMPLEX_FLOAT* prop_d,
	                                             const size_t Nstart,
	                                             const size_t Nstop,
	                                             const size_t dimj,
	                                             const size_t dimi,
	                                             const PRISMATIC_FLOAT_PRECISION* alphaInd_d,
	                                             const cufftHandle& plan,
	                                             cudaStream_t& stream);

	 void getMultisliceProbe_GPU_streaming(Parameters<PRISMATIC_FLOAT_PRECISION>& pars,
	                                      PRISMATIC_CUDA_COMPLEX_FLOAT* trans_d,
	                                      const std::complex<PRISMATIC_FLOAT_PRECISION>* trans_ph,
	                                      PRISMATIC_CUDA_COMPLEX_FLOAT* PsiProbeInit_d,
	                                      PRISMATIC_CUDA_COMPLEX_FLOAT* psi_ds,
	                                      PRISMATIC_FLOAT_PRECISION* output_ph,
	                                      PRISMATIC_FLOAT_PRECISION* psiIntensity_ds,
	                                      PRISMATIC_FLOAT_PRECISION* integratedOutput_ds,
	                                      const PRISMATIC_FLOAT_PRECISION* qya_d,
	                                      const PRISMATIC_FLOAT_PRECISION* qxa_d,
	                                      const PRISMATIC_CUDA_COMPLEX_FLOAT* prop_d,
	                                      const size_t& ay,
	                                      const size_t& ax,
	                                      const size_t dimj,
	                                      const size_t dimi,
	                                      const PRISMATIC_FLOAT_PRECISION* alphaInd_d,
	                                      const cufftHandle& plan,
	                                      cudaStream_t& stream);

	 void getMultisliceProbe_GPU_streaming_batch(Parameters<PRISMATIC_FLOAT_PRECISION>& pars,
	                                            PRISMATIC_CUDA_COMPLEX_FLOAT* trans_d,
	                                            const std::complex<PRISMATIC_FLOAT_PRECISION>* trans_ph,
	                                            PRISMATIC_CUDA_COMPLEX_FLOAT* PsiProbeInit_d,
	                                            PRISMATIC_CUDA_COMPLEX_FLOAT* psi_ds,
	                                            PRISMATIC_FLOAT_PRECISION* output_ph,
	                                            PRISMATIC_FLOAT_PRECISION* psiIntensity_ds,
	                                            PRISMATIC_FLOAT_PRECISION* integratedOutput_ds,
	                                            const PRISMATIC_FLOAT_PRECISION* qya_d,
	                                            const PRISMATIC_FLOAT_PRECISION* qxa_d,
	                                            const PRISMATIC_CUDA_COMPLEX_FLOAT* prop_d,
	                                            const size_t Nstart,
	                                            const size_t Nstop,
	                                            const size_t dimj,
	                                            const size_t dimi,
	                                            const PRISMATIC_FLOAT_PRECISION* alphaInd_d,
	                                            const cufftHandle& plan,
	                                            cudaStream_t& stream);
}
#endif //PRISMATIC_MULTISLICE_CUH
