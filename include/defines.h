//
// Created by AJ Pryor on 3/11/17.
//

#ifndef PRISM_DEFINES_H
#define PRISM_DEFINES_H
namespace PRISM {

	enum class Algorithm {
		PRISM, Multislice
	};
}
#ifdef PRISM_ENABLE_GPU
#include "cuComplex.h"
#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include "cufft.h"

#define BLOCK_SIZE1D 1024
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

#ifdef PRISM_ENABLE_DOUBLE_PRECISION
typedef cuDoubleComplex PRISM_CUDA_COMPLEX_FLOAT;
#define PRISM_CUFFT_EXECUTE cufftExecZ2Z
#define PRISM_CUFFT_PLAN_TYPE CUFFT_Z2Z
#define DEBUGMESSAGE "using double for cuda"
#define PRISM_MAKE_CU_COMPLEX make_cuDoubleComplex
#else
#define DEBUGMESSAGE "using float for cuda"
typedef cuFloatComplex PRISM_CUDA_COMPLEX_FLOAT;
#define PRISM_CUFFT_EXECUTE cufftExecC2C
#define PRISM_CUFFT_PLAN_TYPE CUFFT_C2C
#define PRISM_MAKE_CU_COMPLEX make_cuFloatComplex
#endif //PRISM_ENABLE_DOUBLE_PRECISION


#endif //PRISM_ENABLE_GPU
//#define PRISM_ENABLE_DOUBLE_PRECISION
#ifdef PRISM_ENABLE_DOUBLE_PRECISION
#define MESSAGE "DOUBLE PRECISION"
	typedef double PRISM_FLOAT_PRECISION;
#define PRISM_FFTW_PLAN fftw_plan
#define PRISM_FFTW_PLAN_DFT_2D fftw_plan_dft_2d
#define PRISM_FFTW_EXECUTE fftw_execute
#define PRISM_FFTW_DESTROY_PLAN fftw_destroy_plan
#define PRISM_FFTW_COMPLEX fftw_complex
#define PRISM_FFTW_INIT_THREADS fftw_init_threads
#define	PRISM_FFTW_PLAN_WITH_NTHREADS fftw_plan_with_nthreads
#define PRISM_FFTW_CLEANUP_THREADS fftw_cleanup_threads

#else
	typedef float PRISM_FLOAT_PRECISION;
#define MESSAGE "FLOAT PRECISION"
#define PRISM_FFTW_PLAN fftwf_plan
#define PRISM_FFTW_PLAN_DFT_2D fftwf_plan_dft_2d
#define PRISM_FFTW_EXECUTE fftwf_execute
#define PRISM_FFTW_DESTROY_PLAN fftwf_destroy_plan
#define PRISM_FFTW_COMPLEX fftwf_complex
#define PRISM_FFTW_INIT_THREADS fftwf_init_threads
#define	PRISM_FFTW_PLAN_WITH_NTHREADS fftwf_plan_with_nthreads
#define PRISM_FFTW_CLEANUP_THREADS fftwf_cleanup_threads
//#define PRISM_FFTW_INIT_THREADS
//#define	PRISM_FFTW_PLAN_WITH_NTHREADS
//#define PRISM_FFTW_CLEANUP_THREADS
#endif //PRISM_ENABLE_DOUBLE_PRECISION
#endif //PRISM_DEFINES_H
