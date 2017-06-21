// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

#ifndef PRISMATIC_META_H
#define PRISMATIC_META_H
#include <vector>
#include <string>
#include <cstddef>
#include <iostream>
#include "defines.h"
namespace Prismatic{

    enum class StreamingMode{Stream, SingleXfer, Auto};
	template <class T>
	class Metadata{
	public:
		void toString();
		Metadata(){
			interpolationFactorY  = 4;
			interpolationFactorX  = 4;
			filename_atoms        = "/path/to/atoms.txt";
			filename_output       = "output.mrc";
			realspace_pixelSize[0]= 0.1;
			realspace_pixelSize[1]= 0.1;
			potBound              = 1.0;
			numFP                 = 1;
            fpNum                 = 1;
			sliceThickness        = 2.0;
			cellDim               = std::vector<T>{20.0, 20.0, 20.0}; // this is z,y,x format
			tileX                 = 1;
			tileY                 = 1;
			tileZ                 = 1;
			E0                    = 80e3;
			alphaBeamMax          = 24 / 1000.0;
			NUM_GPUS              = 4;
			NUM_STREAMS_PER_GPU   = 3;
			NUM_THREADS           = 12;
			batch_size_target_CPU = 1;
			batch_size_target_GPU = 1;
			batch_size_CPU 		  = batch_size_target_CPU;
			batch_size_GPU        = batch_size_target_GPU;
			gpu_cpu_ratio         = 100; // relative speed of job completion between gpu and cpu, used to determine early stopping point for cpu work
            probe_stepX           = 0.25;
			probe_stepY           = 0.25;
			probeDefocus          = 0.0;
			C3                    = 0.0;
			C5                    = 0.0;
			probeSemiangle        = 20.0 / 1000;
			detector_angle_step   = 1.0 / 1000;
			probeXtilt            = 0;
			probeYtilt            = 0;
			scanWindowXMin        = 0.0;
			scanWindowXMax        = 1.0;
			scanWindowYMin        = 0.0;
			scanWindowYMax        = 1.0;
			srand(time(0));
			random_seed           = rand() % 100000;
			algorithm             = Algorithm::PRISM;
			include_thermal_effects = true;
			also_do_CPU_work      = true;
			save2DOutput          = false;
			save3DOutput          = true;
			save4DOutput          = false;
			user_specified_celldims = false;
			integration_angle_min = 0;
			integration_angle_max = detector_angle_step;
			transfer_mode         = StreamingMode::Auto;
		}
		size_t interpolationFactorY; // PRISM f_y parameter
		size_t interpolationFactorX; // PRISM f_x parameter
		std::string filename_atoms; // filename of txt file containing atoms (x,y,z,Z CSV format -- one atom per line)
		std::string filename_output;// filename of output image
		T realspace_pixelSize[2]; // pixel size
		T potBound; // bounding integration radius for potential calculation
		size_t numFP; // number of frozen phonon configurations to compute
		size_t fpNum; // current frozen phonon number
		T sliceThickness; // thickness of slice in Z
		T probe_stepX;
		T probe_stepY;
		std::vector<T> cellDim; // this is z,y,x format
		size_t tileX, tileY, tileZ; // how many unit cells to repeat in x,y,z
		size_t batch_size_target_CPU; // desired number of probes/beams to propagate simultaneously for CPU
		size_t batch_size_target_GPU; // desired number of probes/beams to propagate simultaneously for GPU
		size_t batch_size_CPU; // actual number of probes/beams to propagate simultaneously for CPU
		size_t batch_size_GPU; // actual number of probes/beams to propagate simultaneously for GPU
		T gpu_cpu_ratio;
		T E0; // electron energy
		T alphaBeamMax; // max semi angle for probe
		T detector_angle_step;
		T probeDefocus;
		T C3;
		T C5;
		T probeSemiangle;
		T probeXtilt;
		T probeYtilt;
		T scanWindowXMin;
		T scanWindowXMax;
		T scanWindowYMin;
		T scanWindowYMax;
		T random_seed;
		size_t NUM_THREADS; // number of CPU threads to use
		size_t NUM_GPUS; // number of GPUs to use
		size_t NUM_STREAMS_PER_GPU; // number of CUDA streams to use per GPU
		Algorithm algorithm;
		bool include_thermal_effects;
		bool also_do_CPU_work; // what fraction of computation to do on the cpu vs gpu
		bool save2DOutput;
		T integration_angle_min;
		T integration_angle_max;
		bool save3DOutput;
		bool save4DOutput;
		bool user_specified_celldims;
		StreamingMode transfer_mode;
	};

	template <class T>
	void Metadata<T>::toString(){
		std::cout << "interpolationFactorX = " << interpolationFactorX << std::endl;
		std::cout << "interpolationFactorY = " << interpolationFactorY << std::endl;
		std::cout << "filename_atoms = " <<  filename_atoms     << std::endl;
		std::cout << "filename_output = " << filename_output  << std::endl;
		std::cout << "realspace_pixelSize[0] = " << realspace_pixelSize[0]<< std::endl;
		std::cout << "realspace_pixelSize[1] = " << realspace_pixelSize[1]<< std::endl;
		std::cout << "potBound = " << potBound << std::endl;
		std::cout << "numFP = " << numFP << std::endl;
		std::cout << "sliceThickness = " << sliceThickness<< std::endl;
		std::cout << "E0 = " << E0 << std::endl;
		std::cout << "alphaBeamMax = " << alphaBeamMax << std::endl;
		std::cout << "NUM_THREADS = " << NUM_THREADS<< std::endl;
		std::cout << "batch_size_target_CPU = " << batch_size_target_CPU<< std::endl;
		std::cout << "batch_size_target_GPU = " << batch_size_target_GPU<< std::endl;
		std::cout << "probe_stepX = " << probe_stepX << std::endl;
		std::cout << "probe_stepY = " << probe_stepY << std::endl;
		std::cout << "cellDim[0] = " << cellDim[0] << std::endl;
		std::cout << "cellDim[1] = " << cellDim[1] << std::endl;
		std::cout << "cellDim[2] = " << cellDim[2] << std::endl;
		std::cout << "tileX = " << tileX << std::endl;
		std::cout << "tileY = " << tileY << std::endl;
		std::cout << "tileZ = " << tileZ << std::endl;
		std::cout << "probeDefocus = " << probeDefocus<< std::endl;
		std::cout << "C3 = " << C3<< std::endl;
		std::cout << "C5 = " << C5<< std::endl;
		std::cout << "probeSemiangle = " << probeSemiangle<< std::endl;
		std::cout << "detector_angle_step = " << detector_angle_step<< std::endl;
		std::cout << "probeXtilt = " << probeXtilt<< std::endl;
		std::cout << "probeYtilt = " << probeYtilt<< std::endl;
		std::cout << "scanWindowXMin = " << scanWindowXMin<< std::endl;
		std::cout << "scanWindowXMax = " << scanWindowXMax<< std::endl;
		std::cout << "scanWindowYMin = " << scanWindowYMin<< std::endl;
		std::cout << "scanWindowYMax = " << scanWindowYMax<< std::endl;
		std::cout << "integration_angle_min = " << integration_angle_min<< std::endl;
		std::cout << "integration_angle_max = " << integration_angle_max<< std::endl;
		std::cout << "random_seed = " << random_seed << std::endl;

		if (include_thermal_effects) {
			std::cout << "include_thermal_effects = true" << std::endl;
		} else {
			std::cout << "include_thermal_effects = false" << std::endl;
		}
		if (also_do_CPU_work) {
			std::cout << "also_do_CPU_work = true" << std::endl;
		} else {
			std::cout << "also_do_CPU_work = false" << std::endl;
		}
		if (save2DOutput) {
			std::cout << "save2DOutput = true" << std::endl;
		} else {
			std::cout << "save2DOutput = false" << std::endl;
		}
		if (save3DOutput) {
			std::cout << "save3DOutput = true" << std::endl;
		} else {
			std::cout << "save3DOutput = false" << std::endl;
		}
		if (save4DOutput) {
			std::cout << "save4DOutput = true" << std::endl;
		} else {
			std::cout << "save4DOutput = false" << std::endl;
		}
        if (algorithm == Prismatic::Algorithm::PRISM){
            std::cout << "Algorithm: PRISM" << std::endl;
        } else {
			std::cout << "Algorithm: Multislice" << std::endl;
		}

#ifdef PRISMATIC_ENABLE_GPU
        std::cout << "NUM_GPUS = " << NUM_GPUS<< std::endl;
        std::cout << "NUM_STREAMS_PER_GPU = " << NUM_STREAMS_PER_GPU<< std::endl;
        std::cout << "also_do_CPU_work = " << also_do_CPU_work << std::endl;
        std::cout << "gpu_cpu_ratio = " << gpu_cpu_ratio  << std::endl;
		if (transfer_mode == Prismatic::StreamingMode::Auto){
			std::cout << "Data Transfer Mode : Auto" << std::endl;
		} else if (transfer_mode == Prismatic::StreamingMode::SingleXfer){
			std::cout << "Data Transfer : Single Transfer" << std::endl;
		} else {
			std::cout << "Data Transfer : Streaming" << std::endl;
		}
#endif // PRISMATIC_ENABLE_GPU
	}

}
#endif //PRISMATIC_META_H
