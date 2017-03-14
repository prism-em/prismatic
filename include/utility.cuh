//
// Created by AJ Pryor on 3/6/17.
//

#ifndef PRISM_UTILITY_CUH
#define PRISM_UTILITY_CUH
#include <cuda_runtime.h>
#include "cuComplex.h"
#include "defines.h"
// define some constants
//__device__ __constant__ PRISM_FLOAT_PRECISION pi       = PI;
//__device__ __constant__ PRISM_CUDA_COMPLEX_FLOAT i     = {0, 1};
//__device__ __constant__ PRISM_CUDA_COMPLEX_FLOAT pi_cx = {PI, 0};
//__device__ __constant__ PRISM_CUDA_COMPLEX_FLOAT minus_2pii = {0, -2*PI};

// computes exp(real(a) + i * imag(a))
__device__ __forceinline__ cuDoubleComplex exp_cx(const cuDoubleComplex a);
__device__ __forceinline__ cuFloatComplex exp_cx(const cuFloatComplex a);
// creates initial probe using existing GPU memory rather than streaming each probe
__global__ void initializePsi_oneNonzero(cuFloatComplex *psi_d, const size_t N, const size_t beamLoc);

__global__ void initializePsi_oneNonzero(cuDoubleComplex *psi_d, const size_t N, const size_t beamLoc);

// multiply two complex arrays
__global__ void multiply_inplace(PRISM_CUDA_COMPLEX_FLOAT* arr,
                                 const PRISM_CUDA_COMPLEX_FLOAT* other,
                                 const size_t N);

// divide two complex arrays
__global__ void divide_inplace(PRISM_CUDA_COMPLEX_FLOAT* arr,
                               const PRISM_FLOAT_PRECISION val,
                               const size_t N);

__global__ void array_subset(const PRISM_CUDA_COMPLEX_FLOAT* psi_d,
                             PRISM_CUDA_COMPLEX_FLOAT* psi_small_d,
                             const size_t* qyInd_d,
                             const size_t* qxInd_d,
                             const size_t dimi,
                             const size_t dimj_small,
                             const size_t dimi_small,
                             const size_t N);

__global__ void initializePsi(PRISM_CUDA_COMPLEX_FLOAT *psi_d,
                              const PRISM_CUDA_COMPLEX_FLOAT* PsiProbeInit_d,
                              const PRISM_FLOAT_PRECISION* qya_d,
                              const PRISM_FLOAT_PRECISION* qxa_d,
                              const size_t N,
                              const PRISM_FLOAT_PRECISION yp,
                              const PRISM_FLOAT_PRECISION xp);

__global__ void abs_squared(PRISM_FLOAT_PRECISION* arr,
                            const PRISM_CUDA_COMPLEX_FLOAT* other,
                            const size_t N);

__global__ void setAll(PRISM_FLOAT_PRECISION *data, PRISM_FLOAT_PRECISION val, size_t N);
#endif // PRISM_UTILITY_CUH