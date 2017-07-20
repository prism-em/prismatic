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

#ifndef PRISMATIC_DEFINES_H
#define PRISMATIC_DEFINES_H


namespace Prismatic {

	enum class Algorithm {
		PRISM, Multislice
	};
}
#ifdef PRISMATIC_ENABLE_GPU
#include "cuComplex.h"
#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include "cufft.h"

//#define BLOCK_SIZE1D 1024
#define BLOCK_SIZE1D 512
#define PI 3.14159265359
// helpful function for checking CUDA errors.
// Source: http://stackoverflow.com/questions/14038589/what-is-the-canonical-way-to-check-for-errors-using-the-cuda-runtime-api
#define cudaErrchk(ans) { GPUAssert((ans), __FILE__, __LINE__); }
inline void GPUAssert(cudaError_t code, const char *file, int line, bool abort=true){
	if (code != cudaSuccess)
	{
		fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}

// helpful function for checking cuFFT errors
#define cufftErrchk(ans) { GPUAssert_cufft((ans), __FILE__, __LINE__); }
inline void GPUAssert_cufft(int code, const char *file, int line, bool abort=true){
	if (code != CUFFT_SUCCESS)
	{
		fprintf(stderr,"GPUassert: %i %s %d\n", code, file, line);
		if (abort) exit(code);
	}
}

#ifdef PRISMATIC_ENABLE_DOUBLE_PRECISION
typedef cuDoubleComplex PRISMATIC_CUDA_COMPLEX_FLOAT;
#define PRISMATIC_CUFFT_EXECUTE cufftExecZ2Z
#define PRISMATIC_CUFFT_PLAN_TYPE CUFFT_Z2Z
#define DEBUGMESSAGE "using double for cuda"
#define PRISMATIC_MAKE_CU_COMPLEX make_cuDoubleComplex
#else
#define DEBUGMESSAGE "using float for cuda"
typedef cuFloatComplex PRISMATIC_CUDA_COMPLEX_FLOAT;
#define PRISMATIC_CUFFT_EXECUTE cufftExecC2C
#define PRISMATIC_CUFFT_PLAN_TYPE CUFFT_C2C
#define PRISMATIC_MAKE_CU_COMPLEX make_cuFloatComplex
#endif //PRISMATIC_ENABLE_DOUBLE_PRECISION

#endif //PRISMATIC_ENABLE_GPU



//#define PRISMATIC_ENABLE_DOUBLE_PRECISION
#ifdef PRISMATIC_ENABLE_DOUBLE_PRECISION
#define MESSAGE "DOUBLE PRECISION"
	typedef double PRISMATIC_FLOAT_PRECISION;
#define PRISMATIC_FFTW_PLAN fftw_plan
#define PRISMATIC_FFTW_PLAN_DFT_2D fftw_plan_dft_2d
#define PRISMATIC_FFTW_PLAN_DFT_BATCH fftw_plan_many_dft
#define PRISMATIC_FFTW_EXECUTE fftw_execute
#define PRISMATIC_FFTW_DESTROY_PLAN fftw_destroy_plan
#define PRISMATIC_FFTW_COMPLEX fftw_complex
#define PRISMATIC_FFTW_INIT_THREADS fftw_init_threads
#define	PRISMATIC_FFTW_PLAN_WITH_NTHREADS fftw_plan_with_nthreads
#define PRISMATIC_FFTW_CLEANUP_THREADS fftw_cleanup_threads

#else
	typedef float PRISMATIC_FLOAT_PRECISION;
#define MESSAGE "FLOAT PRECISION"
#define PRISMATIC_FFTW_PLAN fftwf_plan
#define PRISMATIC_FFTW_PLAN_DFT_2D fftwf_plan_dft_2d
#define PRISMATIC_FFTW_PLAN_DFT_BATCH fftwf_plan_many_dft
#define PRISMATIC_FFTW_EXECUTE fftwf_execute
#define PRISMATIC_FFTW_DESTROY_PLAN fftwf_destroy_plan
#define PRISMATIC_FFTW_COMPLEX fftwf_complex
#define PRISMATIC_FFTW_INIT_THREADS fftwf_init_threads
#define	PRISMATIC_FFTW_PLAN_WITH_NTHREADS fftwf_plan_with_nthreads
#define PRISMATIC_FFTW_CLEANUP_THREADS fftwf_cleanup_threads
#endif //PRISMATIC_ENABLE_DOUBLE_PRECISION

//#ifdef PRISMATIC_BUILDING_GUI
//class prism_progressbar;
//#endif
#endif //PRISMATIC_DEFINES_H
