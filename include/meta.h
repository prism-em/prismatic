//
// Created by AJ Pryor on 3/2/17.
//

#ifndef PRISM_META_H
#define PRISM_META_H
#include <vector>
#include <string>
#include <cstddef>
namespace PRISM{
	template <class T>
	class Metadata{
	public:
		Metadata(){
			interpolationFactor = 5;
			filename_atoms      = "";
			filename_output     = "";
			realspace_pixelSize = 100.0 / 1000.0;
			potBound = 1.0;
			numFP = 8.0 / 8.0;
			sliceThickness = 2;
			cellDim = std::vector<size_t>{1,1,1}; // this is z,y,x format
			E0 = 80e3;
			alphaBeamMax = 24 / 1000.0;
			NUM_GPUS = 1;
			NUM_STREAMS_PER_GPU = 2;
			NUM_THREADS = 12;
			algorithm = 0; // 0 PRISM; 1 Multislice
			cpu_gpu_ratio = 0.1;
		}
		size_t interpolationFactor;
		std::string filename_atoms;
		std::string filename_output;
		T realspace_pixelSize;

		T potBound;
		size_t numFP;
		T sliceThickness;
		T cpu_gpu_ratio; // what fraction of computation to do on the cpu vs gpu
		std::vector<size_t> cellDim; // this is z,y,x format
		T E0;
		T alphaBeamMax;
		size_t NUM_THREADS;
		size_t NUM_GPUS;
		size_t NUM_STREAMS_PER_GPU;
		int algorithm;
	};
}
#endif //PRISM_META_H
