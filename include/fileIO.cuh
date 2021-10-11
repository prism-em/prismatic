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

#ifndef PRISMATIC_FILEIO_CUH
#define PRISMATIC_FILEIO_CUH
#include "cuComplex.h"
#include "configure.h"
#include "utility.h"
__device__ __forceinline__ cuDoubleComplex exp_cx(const cuDoubleComplex a);
__device__ __forceinline__ cuFloatComplex exp_cx(const cuFloatComplex a);

void formatOutput_GPU_integrate(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
    PRISMATIC_FLOAT_PRECISION *psiIntensity_ds,
    const PRISMATIC_FLOAT_PRECISION *alphaInd_d,
    PRISMATIC_FLOAT_PRECISION *output_ph,
    PRISMATIC_FLOAT_PRECISION *integratedOutput_ds,
    const PRISMATIC_FLOAT_PRECISION* qya_d,
    const PRISMATIC_FLOAT_PRECISION* qxa_d,
    const size_t currentSlice,
    const size_t ay,
    const size_t ax,
    const size_t& dimj,
    const size_t& dimi,
    const cudaStream_t& stream,
    const long& scale = 1);

void formatOutput_GPU_c_integrate(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
    PRISMATIC_CUDA_COMPLEX_FLOAT *psi,
    PRISMATIC_FLOAT_PRECISION *psiIntensity_ds,
    const PRISMATIC_FLOAT_PRECISION *alphaInd_d,
    PRISMATIC_FLOAT_PRECISION *output_ph,
    PRISMATIC_FLOAT_PRECISION *integratedOutput_ds,
    const PRISMATIC_FLOAT_PRECISION* qya_d,
    const PRISMATIC_FLOAT_PRECISION* qxa_d,
    const size_t currentSlice,
    const size_t ay,
    const size_t ax,
    const size_t& dimj,
    const size_t& dimi,
    const cudaStream_t& stream,
    const long& scale = 1);

#endif