// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

//#include "ArrayND.h"
#include "params.h"
//#include "meta.h"
#include "configure.h"
//#include "PRISM_entry.h"

using namespace std;
int main(int argc, const char** argv) {
	PRISM::Metadata<PRISM_FLOAT_PRECISION> prism_meta;
	prism_meta.interpolationFactor = (argc>2) ? (size_t)atoi(argv[2]) : 50;
	prism_meta.algorithm = PRISM::Algorithm::PRISM;
	if (argc > 3){
		prism_meta.algorithm = string(argv[3]) == "m"?PRISM::Algorithm::Multislice:PRISM::Algorithm::PRISM;
	}
	if (argc > 4){
		prism_meta.dxy = (PRISM_FLOAT_PRECISION)atof(argv[4]);
	}
	if (argc > 5){
		prism_meta.realspace_pixelSize = (PRISM_FLOAT_PRECISION)atof(argv[5]);
	}
	if (argc > 6){
		prism_meta.sliceThickness = (PRISM_FLOAT_PRECISION)atof(argv[6]);
	}
	if (argc > 7){
		prism_meta.NUM_THREADS = atoi(argv[7]);
	}
	if (argc > 8){
		prism_meta.NUM_GPUS = atoi(argv[8]);
	}
	std::string filename;
	if (argc>1) {
		prism_meta.filename_atoms = std::string(argv[1]);
	} else{
		cout << "PRISM: Correct syntax is prism filename [interpolation factor] [algorithm] [probe spacing]  [pixel size] [slice thickness] [num threads] [num gpus]" << endl;
		return 0;
	}
	prism_meta.filename_output = "stack.mrc";
	prism_meta.toString();
	PRISM::configure(prism_meta);
	PRISM::execute_plan(prism_meta);
	return 0;
}
