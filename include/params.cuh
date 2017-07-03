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

#ifndef PRISMATIC_PARAMS_CUH
#define PRISMATIC_PARAMS_CUH
#include "defines.h"
#include "configure.h"
namespace Prismatic {
	template<class T>
	class CudaParameters {
	public:
		PRISMATIC_CUDA_COMPLEX_FLOAT **permutedScompact_d;
		PRISMATIC_CUDA_COMPLEX_FLOAT **PsiProbeInit_d;
		PRISMATIC_FLOAT_PRECISION    **qxaReduce_d;
		PRISMATIC_FLOAT_PRECISION    **qyaReduce_d;
		size_t 						 **yBeams_d;
		size_t 						 **xBeams_d;
		PRISMATIC_FLOAT_PRECISION    **alphaInd_d;
		PRISMATIC_CUDA_COMPLEX_FLOAT **psi_ds;
		PRISMATIC_CUDA_COMPLEX_FLOAT **psiSmall_ds;
		PRISMATIC_CUDA_COMPLEX_FLOAT **phaseCoeffs_ds;
		PRISMATIC_FLOAT_PRECISION    **psiIntensity_ds;
		long						 **y_ds;
		long 						 **x_ds;
		PRISMATIC_FLOAT_PRECISION    **integratedOutput_ds;
		cufftHandle 			     *cufftPlans;
		cufftHandle					 *cufftPlansSmall;
		cudaStream_t 				 *streams;
		size_t 						 **qxInd_ph;
		size_t 					     **qyInd_ph;
		size_t 						 **beamsIndex_ph;
        PRISMATIC_CUDA_COMPLEX_FLOAT **trans_d;
        PRISMATIC_CUDA_COMPLEX_FLOAT **prop_d;
        size_t 					     **qxInd_d;
        size_t 						 **qyInd_d;
        size_t 					     **beamsIndex_d;


		// pinned memory buffers on host
		PRISMATIC_FLOAT_PRECISION               **output_ph;
		std::complex<PRISMATIC_FLOAT_PRECISION> *permutedScompact_ph;
		std::complex<PRISMATIC_FLOAT_PRECISION> *PsiProbeInit_ph;
		PRISMATIC_FLOAT_PRECISION               *qxaReduce_ph;
		PRISMATIC_FLOAT_PRECISION               **qxa_ph;
		PRISMATIC_FLOAT_PRECISION               **qxa_d;
		PRISMATIC_FLOAT_PRECISION               *qyaReduce_ph;
		PRISMATIC_FLOAT_PRECISION               **qya_ph;
		PRISMATIC_FLOAT_PRECISION               **qya_d;
		PRISMATIC_FLOAT_PRECISION               *alphaInd_ph;
		size_t                                  *xBeams_ph;
		size_t                                  *yBeams_ph;
		std::complex<PRISMATIC_FLOAT_PRECISION> **ScompactSlice_ph;
		std::complex<PRISMATIC_FLOAT_PRECISION> *trans_ph;
		std::complex<PRISMATIC_FLOAT_PRECISION> *prop_ph;

	};
}
#endif //PRISMATIC_PARAMS_CUH