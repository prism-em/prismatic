//#include "ArrayND.h"
#include "params.h"
//#include "meta.h"
#include "configure.h"
//#include "PRISM_entry.h"

using namespace std;
int main(int argc, const char** argv) {
	PRISM::Metadata<PRISM_FLOAT_PRECISION> prism_meta;
	prism_meta.interpolationFactor = (argc>2) ? (size_t)atoi(argv[2]) : 50;


	prism_meta.potBound = 1.0;
	prism_meta.numFP = 8.0 / 8.0;
	prism_meta.sliceThickness = 2;
        prism_meta.cellDim = vector<size_t>{80,100,100}; // this is z,y,x format
	prism_meta.realspace_pixelSize = 100.0/1000.0/1;
//prism_meta.cellDim = vector<size_t>{40,16,16}; // this is z,y,x format
//prism_meta.cellDim = vector<size_t>{37,12,14}; // this is z,y,x format
//	prism_meta.cellDim = vector<size_t>{37,14,14}; // this is z,y,x format
//	prism_meta.cellDim = vector<size_t>{10,32,32}; // this is z,y,x format
//	prism_meta.cellDim = vector<size_t>{2*40,2*16,2*16}; // this is z,y,x format
//	prism_meta.cellDim = vector<size_t>{80,4*16,4*16}; // this is z,y,x format
	//prism_meta.cellDim = vector<size_t>{1000,8,8}; // this is z,y,x format
//	prism_meta.cellDim = vector<size_t>{20,16*16,16*16}; // this is z,y,x format
//	prism_meta.cellDim = vector<size_t>{40,64,64}; // this is z,y,x format
	prism_meta.E0 = 80e3;
	prism_meta.alphaBeamMax = 24 / 1000.0;
  int nDevices;
#ifdef PRISM_ENABLE_GPU
  cudaGetDeviceCount(&nDevices);
cout << "ndevices = " << nDevices << endl;
#endif	
	prism_meta.algorithm = PRISM::Algorithm::PRISM;
//	prism_meta.algorithm = PRISM::Algorithm::Multislice;

	std::string filename;
	if (argc>1) {
		prism_meta.filename_atoms = std::string(argv[1]);
	} else{
		cout << "PRISM: Correct syntax is prism filename [interpolation factor]" << endl;
		return 0;
	}
	prism_meta.filename_output = "prism_image.mrc";

	PRISM::configure(prism_meta);
	PRISM::execute_plan(prism_meta);
	return 0;
}
