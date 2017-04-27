// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

#ifndef PRISM_META_H
#define PRISM_META_H
#include <vector>
#include <string>
#include <cstddef>
#include <iostream>
#include "defines.h"
namespace PRISM{

    enum class StreamingMode{Stream, SingleXfer, Auto};
	template <class T>
	class Metadata{
	public:
		void toString();
		Metadata(){
			interpolationFactor   = 5;
			filename_atoms        = "/path/to/atoms.txt";
			filename_output       = "output.mrc";
			realspace_pixelSize   = 0.1;
			potBound              = 1.0;
			numFP                 = 1;
            fpNum                 = 1;
			sliceThickness        = 2.0;
			cellDim               = std::vector<size_t>{20, 20, 20}; // this is z,y,x format
			E0                    = 80e3;
			alphaBeamMax          = 24 / 1000.0;
			NUM_GPUS              = 4;
			NUM_STREAMS_PER_GPU   = 3;
			NUM_THREADS           = 4;
			gpu_cpu_ratio         = 20; // relative speed of job completion between gpu and cpu, used to determine early stopping point for cpu work
            probe_step            = 0.25;

			//add to cli
			detector_angle_step   = 2.5 / 1000;
			probeXtilt            = 0;
			probeYtilt            = 0;
			probeDefocus          = 0.0;
			probeSemiangle        = 20.0 / 1000;
			algorithm             = Algorithm::PRISM;
			also_do_CPU_work      = true;
			save2DOutput          = false;
			save3DOutput          = true;
			save4DOutput          = false;
			integration_angle_min = 0;
			integration_angle_max = detector_angle_step;
			transfer_mode         = StreamingMode::Auto;
		}
		size_t interpolationFactor; // PRISM f parameter
		std::string filename_atoms; // filename of txt file containing atoms (x,y,z,Z CSV format -- one atom per line)
		std::string filename_output;// filename of output image
		T realspace_pixelSize; // pixel size
		T potBound; // bounding integration radius for potential calculation
		size_t numFP; // number of frozen phonon configurations to compute
		size_t fpNum; // current frozen phonon number
		T sliceThickness; // thickness of slice in Z
		T probe_step;
		std::vector<size_t> cellDim; // this is z,y,x format
		T gpu_cpu_ratio;
		T E0; // electron energy
		T alphaBeamMax; // max semi angle for probe
		T detector_angle_step;
		T probeDefocus;
		T probeSemiangle;
		T probeXtilt;
		T probeYtilt;
		size_t NUM_THREADS; // number of CPU threads to use
		size_t NUM_GPUS; // number of GPUs to use
		size_t NUM_STREAMS_PER_GPU; // number of CUDA streams to use per GPU
		Algorithm algorithm;
		bool also_do_CPU_work; // what fraction of computation to do on the cpu vs gpu
		bool save2DOutput;
		T integration_angle_min;
		T integration_angle_max;
		bool save3DOutput;
		bool save4DOutput;
		StreamingMode transfer_mode;
	};

	template <class T>
	void Metadata<T>::toString(){
		std::cout << "interpolationFactor = " << interpolationFactor << std::endl;
		std::cout << "filename_atoms = " <<  filename_atoms     << std::endl;
		std::cout << "filename_output = " << filename_output  << std::endl;
		std::cout << "realspace_pixelSize = " << realspace_pixelSize<< std::endl;
		std::cout << "potBound = " << potBound << std::endl;
		std::cout << "numFP = " << numFP << std::endl;
		std::cout << "sliceThickness = " << sliceThickness<< std::endl;
		std::cout << "E0 = " << E0 << std::endl;
		std::cout << "alphaBeamMax = " << alphaBeamMax << std::endl;
		std::cout << "NUM_THREADS = " << NUM_THREADS<< std::endl;
		std::cout << "probe_step = " << probe_step << std::endl;
		std::cout << "cellDim[0] = " << cellDim[0] << std::endl;
		std::cout << "cellDim[1] = " << cellDim[1] << std::endl;
		std::cout << "cellDim[2] = " << cellDim[2] << std::endl;
        if (algorithm == PRISM::Algorithm::PRISM){
            std::cout << "Algorithm: PRISM" << std::endl;
        } else {
			std::cout << "Algorithm: Multislice" << std::endl;
		}

#ifdef PRISM_ENABLE_GPU
        std::cout << "NUM_GPUS = " << NUM_GPUS<< std::endl;
        std::cout << "NUM_STREAMS_PER_GPU = " << NUM_STREAMS_PER_GPU<< std::endl;
        std::cout << "also_do_CPU_work = " << also_do_CPU_work << std::endl;
        std::cout << "gpu_cpu_ratio = " << gpu_cpu_ratio  << std::endl;
		if (transfer_mode == PRISM::StreamingMode::Auto){
			std::cout << "Data Transfer Mode : Auto" << std::endl;
		} else if (transfer_mode == PRISM::StreamingMode::SingleXfer){
			std::cout << "Data Transfer : Single Transfer" << std::endl;
		} else {
			std::cout << "Data Transfer : Streaming" << std::endl;
		}
#endif // PRISM_ENABLE_GPU
	}

}
#endif //PRISM_META_H
