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

#ifndef PRISMATIC_PRISM02_CUH
#define PRISMATIC_PRISM02_CUH
#include "params.h"
#include "configure.h"
#include <complex>
namespace Prismatic{
	 void fill_Scompact_GPU_singlexfer(Parameters<PRISMATIC_FLOAT_PRECISION> &pars);
	 void fill_Scompact_GPU_streaming(Parameters <PRISMATIC_FLOAT_PRECISION> &pars);
	 void propagatePlaneWave_GPU_singlexfer(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
	                                       PRISMATIC_CUDA_COMPLEX_FLOAT* trans_d,
	                                       PRISMATIC_CUDA_COMPLEX_FLOAT* psi_d,
	                                       PRISMATIC_CUDA_COMPLEX_FLOAT* psi_small_d,
	                                       std::complex<PRISMATIC_FLOAT_PRECISION>* Scompact_slice_ph,
	                                       const size_t* qyInd_d,
	                                       const size_t* qxInd_d,
	                                       const PRISMATIC_CUDA_COMPLEX_FLOAT* prop_d,
	                                       const size_t* beamsIndex,
	                                       const size_t beamNumber,
	                                       const cufftHandle& plan,
	                                       const cufftHandle& plan_small,
	                                       cudaStream_t& stream);

	 void propagatePlaneWave_GPU_singlexfer_batch(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
	                                             PRISMATIC_CUDA_COMPLEX_FLOAT* trans_d,
	                                             PRISMATIC_CUDA_COMPLEX_FLOAT* psi_d,
	                                             PRISMATIC_CUDA_COMPLEX_FLOAT* psi_small_d,
	                                             std::complex<PRISMATIC_FLOAT_PRECISION>* Scompact_slice_ph,
	                                             const size_t* qyInd_d,
	                                             const size_t* qxInd_d,
	                                             const PRISMATIC_CUDA_COMPLEX_FLOAT* prop_d,
	                                             const size_t* beamsIndex,
	                                             const size_t beamNumber,
	                                             const size_t stopBeam,
	                                             const cufftHandle& plan,
	                                             const cufftHandle& plan_small,
	                                             cudaStream_t& stream);

	 void propagatePlaneWave_GPU_streaming(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
	                                      PRISMATIC_CUDA_COMPLEX_FLOAT* trans_d,
	                                      const std::complex<PRISMATIC_FLOAT_PRECISION> *trans_ph,
	                                      PRISMATIC_CUDA_COMPLEX_FLOAT* psi_d,
	                                      PRISMATIC_CUDA_COMPLEX_FLOAT* psi_small_d,
	                                      std::complex<PRISMATIC_FLOAT_PRECISION>* Scompact_slice_ph,
	                                      const size_t* qyInd_d,
	                                      const size_t* qxInd_d,
	                                      const PRISMATIC_CUDA_COMPLEX_FLOAT* prop_d,
	                                      const size_t* beamsIndex,
	                                      const size_t beamNumber,
	                                      const cufftHandle& plan,
	                                      const cufftHandle& plan_small,
	                                      cudaStream_t& stream);

	 void propagatePlaneWave_GPU_streaming_batch(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
	                                            PRISMATIC_CUDA_COMPLEX_FLOAT* trans_d,
	                                            const std::complex<PRISMATIC_FLOAT_PRECISION> *trans_ph,
	                                            PRISMATIC_CUDA_COMPLEX_FLOAT* psi_d,
	                                            PRISMATIC_CUDA_COMPLEX_FLOAT* psi_small_d,
	                                            std::complex<PRISMATIC_FLOAT_PRECISION>* Scompact_slice_ph,
	                                            const size_t* qyInd_d,
	                                            const size_t* qxInd_d,
	                                            const PRISMATIC_CUDA_COMPLEX_FLOAT* prop_d,
	                                            const size_t* beamsIndex,
	                                            const size_t beamNumber,
	                                            const size_t stopBeam,
	                                            const cufftHandle& plan,
	                                            const cufftHandle& plan_small,
	                                            cudaStream_t& stream);

}
#endif //PRISMATIC_PRISM02_CUH