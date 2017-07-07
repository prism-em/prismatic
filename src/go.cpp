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

//#include <iostream>
//#include <stdlib.h>
//#include <algorithm>
//#include <vector>
#include "PRISM_entry.h"
#include "Multislice_entry.h"
//#include "configure.h"
//#include "ArrayND.h"
#include "params.h"
#include "go.h"


namespace Prismatic{
	void go(Metadata<PRISMATIC_FLOAT_PRECISION> meta){
		// configure simulation behavior
		Prismatic::configure(meta);

		// execute simulation
		Prismatic::execute_plan(meta);
	}
}
