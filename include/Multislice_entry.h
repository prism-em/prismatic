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

#ifndef PRISM_MULTISLICE_ENTRY_H
#define PRISM_MULTISLICE_ENTRY_H
#include "meta.h"
#include "params.h"
#include "ArrayND.h"
#include "configure.h"
#include "Multislice_calcOutput.h"
#include "PRISM01_calcPotential.h"
#include "PRISM02_calcSMatrix.h"
#include <algorithm>


namespace Prismatic{
	void Multislice_entry(Metadata<PRISMATIC_FLOAT_PRECISION>& meta);
    
    void Multislice_entry_pars(Parameters<PRISMATIC_FLOAT_PRECISION> &pars);

	void Multislice_runFP(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, size_t fpNum);

	void Multislice_series_runFP(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, size_t fpNum);
}
#endif //PRISM_MULTISLICE_ENTRY_H
