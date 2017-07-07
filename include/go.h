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

#ifndef PRISM_GO_H
#define PRISM_GO_H
#include "configure.h"
#include "params.h"

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && defined(BUILDING_CUPRISMATIC)
#ifdef CUPRISMATIC_EXPORT
        #define PRISMATIC_API __declspec(dllexport)
    #else
        #define PRISMATIC_API __declspec(dllimport)
    #endif //CUPRISMATIC_BUILDING
#else
#define PRISMATIC_API
#endif

namespace Prismatic{
	PRISMATIC_API void go(Metadata<PRISMATIC_FLOAT_PRECISION> meta);
}
#endif //PRISM_GO_H
