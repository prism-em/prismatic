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

#ifndef PRISM_PROBE_H
#define PRISM_PROBE_H
#include <string>
#include <vector>
#include "defines.h"

namespace Prismatic
{

void readProbes(const std::string &filename, 
                std::vector<PRISMATIC_FLOAT_PRECISION> &xprobes,
                std::vector<PRISMATIC_FLOAT_PRECISION> &yprobes);

} // namespace Prismatic
#endif //PRISM_PROBE_H