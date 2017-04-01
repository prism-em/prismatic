//
// Created by AJ Pryor on 3/2/17.
//

#ifndef PRISM_META_H
#define PRISM_META_H
#include <vector>
#include <string>
#include <cstddef>
#include <iostream>
#include "defines.h"
namespace PRISM{

	template <class T>
	class Metadata{
	public:
		void toString();
		Metadata(){
			interpolationFactor = 5;
			filename_atoms      = "/path/to/atoms.txt";
			filename_output     = "/path/to/output.mrc";
			realspace_pixelSize = 100.0 / 1000.0;
			potBound = 1.0;
			numFP = 8.0 / 8.0;
			sliceThickness = 2;
			cellDim = std::vector<size_t>{20, 20, 20}; // this is z,y,x format
			E0 = 80e3;
			alphaBeamMax = 24 / 1000.0;
			NUM_GPUS = 4;
			NUM_STREAMS_PER_GPU = 3;
			NUM_THREADS = 48;
			algorithm = Algorithm::PRISM; // 0 PRISM; 1 Multislice
//			also_do_CPU_work = true;
			also_do_CPU_work = false;
			//cpu_gpu_ratio = 1 / 15; // relative speed of job completion between cpu and gpu, used to determine early stopping point for cpu work
			cpu_gpu_ratio = 0.05; // relative speed of job completion between cpu and gpu, used to determine early stopping point for cpu work
//			stream_data = true;
			stream_data = false;
			dxy = 0.25;
		}
		size_t interpolationFactor; // PRISM f parameter
		std::string filename_atoms; // filename of txt file containing atoms (x,y,z,Z CSV format -- one atom per line)
		std::string filename_output;// filename of output image
		T realspace_pixelSize; // pixel size
		T potBound; // bounding integration radius for potential calculation
		size_t numFP; // number of frozen phonon configurations to compute
		T sliceThickness; // thickness of slice in Z
		T dxy;
		bool also_do_CPU_work; // what fraction of computation to do on the cpu vs gpu
		bool stream_data;
		std::vector<size_t> cellDim; // this is z,y,x format
		T cpu_gpu_ratio;
		T E0; // electron energy
		T alphaBeamMax; // max semi angle for probe
		size_t NUM_THREADS; // number of CPU threads to use
		size_t NUM_GPUS; // number of GPUs to use
		size_t NUM_STREAMS_PER_GPU; // number of CUDA streams to use per GPU
		Algorithm algorithm;
	};

	template <class T>
	void Metadata<T>::toString(){
		std::cout << "interpolationFactor = 5" << interpolationFactor << std::endl;
		std::cout << "filename_atoms = " <<  filename_atoms     << std::endl;
		std::cout << "filename_output = " << filename_output  << std::endl;
		std::cout << "realspace_pixelSize = " << realspace_pixelSize<< std::endl;
		std::cout << "potBound = " << potBound << std::endl;
		std::cout << "numFP = " << numFP << std::endl;
		std::cout << "sliceThickness = " << sliceThickness<< std::endl;
		std::cout << "E0 = " << E0 << std::endl;
		std::cout << "alphaBeamMax = " << alphaBeamMax << std::endl;
		std::cout << "NUM_GPUS = " << NUM_GPUS<< std::endl;
		std::cout << "NUM_STREAMS_PER_GPU = " << NUM_STREAMS_PER_GPU<< std::endl;
		std::cout << "NUM_THREADS = " << NUM_THREADS<< std::endl;
		std::cout << "also_do_CPU_work = " << also_do_CPU_work << std::endl;
		std::cout << "cpu_gpu_ratio = " << cpu_gpu_ratio  << std::endl;
		std::cout << "stream_data = " << stream_data  << std::endl;
		std::cout << "dxy = " << dxy << std::endl;
		std::cout << "cellDim[0] = " << cellDim[0] << std::endl;
		std::cout << "cellDim[1] = " << cellDim[1] << std::endl;
		std::cout << "cellDim[2] = " << cellDim[2] << std::endl;
	}
}
#endif //PRISM_META_H
