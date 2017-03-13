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

#include <cuda_runtime.h>

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
#endif //PRISM_ENABLE_DOUBLE_PRECISION
#endif //PRISM_DEFINES_H
