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

// This header file is for the shared library cuPrismatic, which is used for exposing the CUDA code in Prismatic to its Python package, PyPrismatic.
// The cuPrismatic library is necessary for Python because distutils does not work with nvcc, so CMake is used to create an additional shared library that can
// then be linked against in the Python package.
#ifndef CUPRISMATIC_H
#define CUPRISMATIC_H

//#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && defined(BUILDING_CUPRISMATIC)
//    #ifdef CUPRISMATIC_EXPORT
//        #define PRISMATIC_API __declspec(dllexport)
//    #else
//        #define PRISMATIC_API __declspec(dllimport)
//    #endif //CUPRISMATIC_BUILDING
//#else
//    #define PRISMATIC_API
//#endif

//#include "PRISM02_calcSMatrix.cuh"
//#include "PRISM03_calcOutput.cuh"
//#include "Multislice_calcOutput.cuh"
//#include "params.cuh"
//#include "utility.cuh"
#include "go.h"

#endif // CUPRISMATIC_H