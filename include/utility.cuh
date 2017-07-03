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

#ifndef PRISMATIC_UTILITY_CUH
#define PRISMATIC_UTILITY_CUH
#define PI 3.14159265359
#include "cuComplex.h"
#include "configure.h"
#include "utility.h"
__device__ __forceinline__ cuDoubleComplex exp_cx(const cuDoubleComplex a);
__device__ __forceinline__ cuFloatComplex exp_cx(const cuFloatComplex a);

// creates initial probe using existing GPU memory rather than streaming each probe
__global__ void initializePsi_oneNonzero(cuFloatComplex *psi_d, const size_t N, const size_t beamLoc);


__global__ void initializePsi_oneNonzero(cuDoubleComplex *psi_d, const size_t N, const size_t beamLoc);


// multiply two complex arrays
__global__ void multiply_inplace(cuDoubleComplex* arr,
                                 const cuDoubleComplex* other,
                                 const size_t N);


// multiply two complex arrays
__global__ void multiply_inplace(cuFloatComplex* arr,
                                 const cuFloatComplex* other,
                                 const size_t N);


// multiply two complex arrays
__global__ void multiply_cx(cuDoubleComplex* arr,
                            const cuDoubleComplex* other,
                            const size_t N);


// multiply two complex arrays
__global__ void multiply_cx(cuFloatComplex* arr,
                            const cuFloatComplex* other,
                            const size_t N);


//// divide two complex arrays
//__global__ void divide_inplace(PRISMATIC_CUDA_COMPLEX_FLOAT* arr,
//                               const PRISMATIC_FLOAT_PRECISION val,
//                               const size_t N){
//	int idx = threadIdx.x + blockDim.x*blockIdx.x;
//	if (idx < N) {
//		arr[idx].x /= val;
//		arr[idx].y /= val;
//	}
//}

__global__ void divide_inplace(cuDoubleComplex* arr,
                               const cuDoubleComplex val,
                               const size_t N);


__global__ void divide_inplace(cuFloatComplex* arr,
                               const cuFloatComplex val,
                               const size_t N);


// set all array values to val
__global__ void setAll(double *data, double val, size_t N);


// set all array values to val
__global__ void setAll(float *data, float val, size_t N);


// creates initial probe using existing GPU memory rather than streaming each probe
__global__ void initializePsi(cuDoubleComplex *psi_d,
                              const cuDoubleComplex* PsiProbeInit_d,
                              const double* qya_d,
                              const double* qxa_d,
                              const size_t N,
                              const double yp,
                              const double xp);


// creates initial probe using existing GPU memory rather than streaming each probe
__global__ void initializePsi(cuFloatComplex *psi_d,
                              const cuFloatComplex* PsiProbeInit_d,
                              const float* qya_d,
                              const float* qxa_d,
                              const size_t N,
                              const float yp,
                              const float xp);



// compute modulus squared of other and store in arr
__global__ void abs_squared(double* arr,
                            const cuDoubleComplex* other,
                            const size_t N);


// compute modulus squared of other and store in arr
__global__ void abs_squared(float* arr,
                            const cuFloatComplex* other,
                            const size_t N);


__global__ void array_subset(const cuDoubleComplex* psi_d,
                             cuDoubleComplex* psi_small_d,
                             const size_t* qyInd_d,
                             const size_t* qxInd_d,
                             const size_t dimi,
                             const size_t dimj_small,
                             const size_t dimi_small);

__global__ void array_subset(const cuFloatComplex* psi_d,
                             cuFloatComplex* psi_small_d,
                             const size_t* qyInd_d,
                             const size_t* qxInd_d,
                             const size_t dimi,
                             const size_t dimj_small,
                             const size_t dimi_small);

__global__ void shiftIndices(long* vec_out, const long by, const long imageSize, const long N);
__global__ void zeroIndices(long* vec_out, const long N);
__global__ void resetIndices(long* vec_out, const long N);
__global__ void computePhaseCoeffs(PRISMATIC_CUDA_COMPLEX_FLOAT* phaseCoeffs,
                                   const PRISMATIC_CUDA_COMPLEX_FLOAT *PsiProbeInit_d,
                                   const PRISMATIC_FLOAT_PRECISION * qyaReduce_d,
                                   const PRISMATIC_FLOAT_PRECISION * qxaReduce_d,
                                   const size_t *yBeams_d,
                                   const size_t *xBeams_d,
                                   const PRISMATIC_FLOAT_PRECISION yp,
                                   const PRISMATIC_FLOAT_PRECISION xp,
                                   const PRISMATIC_FLOAT_PRECISION yTiltShift,
                                   const PRISMATIC_FLOAT_PRECISION xTiltShift,
                                   const size_t dimi,
                                   const size_t numBeams);

//template <size_t BlockSizeX>
//__global__ void scaleReduceS(const PRISMATIC_CUDA_COMPLEX_FLOAT *permutedScompact_d,
//                             const PRISMATIC_CUDA_COMPLEX_FLOAT *phaseCoeffs_ds,
//                             PRISMATIC_CUDA_COMPLEX_FLOAT *psi_ds,
//                             const long *y_ds,
//                             const long* x_ds,
//                             const size_t numberBeams,
//                             const size_t dimj_S,
//                             const size_t dimk_S,
//                             const size_t dimj_psi,
//                             const size_t dimi_psi);

__global__ void integrateDetector(const float* psiIntensity_ds,
                       const float* alphaInd_d,
                       float* integratedOutput,
                       const size_t N,
                       const size_t num_integration_bins);

void formatOutput_GPU_integrate(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
                                PRISMATIC_FLOAT_PRECISION *psiIntensity_ds,
                                const PRISMATIC_FLOAT_PRECISION *alphaInd_d,
                                PRISMATIC_FLOAT_PRECISION *stack_ph,
                                PRISMATIC_FLOAT_PRECISION *integratedOutput_ds,
                                const size_t ay,
                                const size_t ax,
                                const size_t& dimj,
                                const size_t& dimi,
                                const cudaStream_t& stream,
                                const long& scale = 1);

__global__ void multiply_cxarr_scalar(cuDoubleComplex* arr,
                                      const cuDoubleComplex val,
                                      const size_t N);

__global__ void multiply_cxarr_scalar(cuFloatComplex* arr,
                                      const cuFloatComplex val,
                                      const size_t N);

__global__ void multiply_arr_scalar(double* arr,
                                    const double val,
                                    const size_t N);

__global__ void multiply_arr_scalar(float* arr,
                                    const float val,
                                    const size_t N);

//size_t getNextPower2(const double& val);
//size_t getNextPower2(const float& val);
size_t getNextPower2(const size_t& val);

#if __CUDA_ARCH__ < 600
__device__  double atomicAdd_double(double* address, const double val);
#endif
//__global__ void shiftIndices(size_t* vec, const size_t* vec_in, double by, const size_t N);
#endif // PRISMATIC_UTILITY_CUH