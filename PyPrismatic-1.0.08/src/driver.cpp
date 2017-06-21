// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

#include "params.h"
#include "configure.h"
#include "parseInput.h"

using namespace std;
int main(int argc, const char** argv) {
	Prismatic::Metadata<PRISMATIC_FLOAT_PRECISION> meta;

	// parse command line options
	if (!Prismatic::parseInputs(meta, argc, &argv))return 1;

	// print metadata
    meta.toString();

	// configure simulation behavior
	Prismatic::configure(meta);

	// execute simulation
	Prismatic::execute_plan(meta);
	return 0;
}
