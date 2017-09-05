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

#include "params.h"
#include "go.h"
#include "configure.h"
#include "parseInput.h"

using namespace std;
int main(int argc, const char** argv) {
	Prismatic::Metadata<PRISMATIC_FLOAT_PRECISION> meta;

	// parse command line options
	if (!Prismatic::parseInputs(meta, argc, &argv))return 1;

	Prismatic::printHeader();

	// print metadata
        meta.toString();

	// execute simulation
	Prismatic::go(meta);

	return 0;
}
