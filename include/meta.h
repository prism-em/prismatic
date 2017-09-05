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
		bool operator==(const Metadata<T> other);
	
		Metadata(){
			interpolationFactorY  = 4;
			interpolationFactorX  = 4;
			filenameAtoms        = "/path/to/atoms.txt";
			filenameOutput       = "output.mrc";
			realspacePixelSize[0] = 0.1;
			realspacePixelSize[1] = 0.1;
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
			numGPUs               = 4;
			numStreamsPerGPU      = 3;
			numThreads            = 12;
			batchSizeTargetCPU    = 1;
			batchSizeTargetGPU    = 1;
			batchSizeCPU 		  = batchSizeTargetCPU;
			batchSizeGPU          = batchSizeTargetGPU;
			earlyCPUStopCount     = 100; // relative speed of job completion between gpu and cpu, used to determine early stopping point for cpu work
            probeStepX            = 0.25;
			probeStepY            = 0.25;
			probeDefocus          = 0.0;
			C3                    = 0.0;
			C5                    = 0.0;
			probeSemiangle        = 20.0 / 1000;
			detectorAngleStep     = 1.0 / 1000;
			probeXtilt            = 0;
			probeYtilt            = 0;
			scanWindowXMin        = 0.0;
			scanWindowXMax        = 0.99999;
			scanWindowYMin        = 0.0;
			scanWindowYMax        = 0.99999;
			srand(time(0));
			randomSeed            = rand() % 100000;
			algorithm             = Algorithm::PRISM;
			includeThermalEffects = true;
			includeOccupancy      = true;
			alsoDoCPUWork         = true;
			save2DOutput          = false;
			save3DOutput          = true;
			save4DOutput          = false;
			userSpecifiedCelldims = false;
			integrationAngleMin   = 0;
			integrationAngleMax   = detectorAngleStep;
			transferMode          = StreamingMode::Auto;
		}
		size_t interpolationFactorY; // PRISM f_y parameter
		size_t interpolationFactorX; // PRISM f_x parameter
		std::string filenameAtoms; // filename of txt file containing atoms (x,y,z,Z CSV format -- one atom per line)
		std::string filenameOutput;// filename of output image
		T realspacePixelSize[2]; // pixel size
		T potBound; // bounding integration radius for potential calculation
		size_t numFP; // number of frozen phonon configurations to compute
		size_t fpNum; // current frozen phonon number
		T sliceThickness; // thickness of slice in Z
		T probeStepX;
		T probeStepY;
		std::vector<T> cellDim; // this is z,y,x format
		size_t tileX, tileY, tileZ; // how many unit cells to repeat in x,y,z
		size_t batchSizeTargetCPU; // desired number of probes/beams to propagate simultaneously for CPU
		size_t batchSizeTargetGPU; // desired number of probes/beams to propagate simultaneously for GPU
		size_t batchSizeCPU; // actual number of probes/beams to propagate simultaneously for CPU
		size_t batchSizeGPU; // actual number of probes/beams to propagate simultaneously for GPU
		T earlyCPUStopCount;
		T E0; // electron energy
		T alphaBeamMax; // max semi angle for probe
		T detectorAngleStep;
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
		T randomSeed;
		size_t numThreads; // number of CPU threads to use
		size_t numGPUs; // number of GPUs to use
		size_t numStreamsPerGPU; // number of CUDA streams to use per GPU
		Algorithm algorithm;
		bool includeThermalEffects;
		bool includeOccupancy;
		bool alsoDoCPUWork; // what fraction of computation to do on the cpu vs gpu
		bool save2DOutput;
		T integrationAngleMin;
		T integrationAngleMax;
		bool save3DOutput;
		bool save4DOutput;
		bool userSpecifiedCelldims;
		StreamingMode transferMode;

	};

	template <class T>
	void Metadata<T>::toString(){
		std::cout << "\nSimulation parameters:" << std::endl;
		std::cout << "=====================\n" << std::endl;
		if (algorithm == Prismatic::Algorithm::PRISM){
			std::cout << "Algorithm: PRISM" << std::endl;
			std::cout << "interpolationFactorX = " << interpolationFactorX << std::endl;
			std::cout << "interpolationFactorY = " << interpolationFactorY << std::endl;
		} else {
			std::cout << "Algorithm: Multislice" << std::endl;
		}

		std::cout << "filenameAtoms = " <<  filenameAtoms     << std::endl;
		std::cout << "filenameOutput = " << filenameOutput  << std::endl;
		std::cout << "realspacePixelSize[0] = " << realspacePixelSize[0]<< std::endl;
		std::cout << "realspacePixelSize[1] = " << realspacePixelSize[1]<< std::endl;
		std::cout << "potBound = " << potBound << std::endl;
		std::cout << "numFP = " << numFP << std::endl;
		std::cout << "sliceThickness = " << sliceThickness<< std::endl;
		std::cout << "E0 = " << E0 << std::endl;
		std::cout << "alphaBeamMax = " << alphaBeamMax << std::endl;
		std::cout << "numThreads = " << numThreads<< std::endl;
		std::cout << "batchSizeTargetCPU = " << batchSizeTargetCPU<< std::endl;
		std::cout << "batchSizeTargetGPU = " << batchSizeTargetGPU<< std::endl;
		std::cout << "probeStepX = " << probeStepX << std::endl;
		std::cout << "probeStepY = " << probeStepY << std::endl;
		std::cout << "cellDim[0] = " << cellDim[0] << std::endl;
		std::cout << "cellDim[1] = " << cellDim[1] << std::endl;
		std::cout << "cellDim[2] = " << cellDim[2] << std::endl;
		std::cout << "tileX = " << tileX << std::endl;
		std::cout << "tileY = " << tileY << std::endl;
		std::cout << "tileZ = " << tileZ << std::endl;
		std::cout << "probeDefocus = " << probeDefocus<< std::endl;
		std::cout << "C3 = " << C3 << std::endl;
		std::cout << "C5 = " << C5 << std::endl;
		std::cout << "probeSemiangle = " << probeSemiangle<< std::endl;
		std::cout << "detectorAngleStep = " << detectorAngleStep<< std::endl;
		std::cout << "probeXtilt = " << probeXtilt<< std::endl;
		std::cout << "probeYtilt = " << probeYtilt<< std::endl;
		std::cout << "scanWindowXMin = " << scanWindowXMin<< std::endl;
		std::cout << "scanWindowXMax = " << scanWindowXMax<< std::endl;
		std::cout << "scanWindowYMin = " << scanWindowYMin<< std::endl;
		std::cout << "scanWindowYMax = " << scanWindowYMax<< std::endl;
		std::cout << "integrationAngleMin = " << integrationAngleMin<< std::endl;
		std::cout << "integrationAngleMax = " << integrationAngleMax<< std::endl;
		std::cout << "randomSeed = " << randomSeed << std::endl;

		if (includeOccupancy) {
			std::cout << "includeOccupancy = true" << std::endl;
		} else {
			std::cout << "includeOccupancy = false" << std::endl;
		}
		if (includeThermalEffects) {
			std::cout << "includeThermalEffects = true" << std::endl;
		} else {
			std::cout << "includeThermalEffects = false" << std::endl;
		}
		if (alsoDoCPUWork) {
			std::cout << "alsoDoCPUWork = true" << std::endl;
		} else {
			std::cout << "alsoDoCPUWork = false" << std::endl;
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


#ifdef PRISMATIC_ENABLE_GPU
        std::cout << "numGPUs = " << numGPUs<< std::endl;
        std::cout << "numStreamsPerGPU = " << numStreamsPerGPU<< std::endl;
        std::cout << "alsoDoCPUWork = " << alsoDoCPUWork << std::endl;
        std::cout << "earlyCPUStopCount = " << earlyCPUStopCount  << std::endl;
		if (transferMode == Prismatic::StreamingMode::Auto){
			std::cout << "Data Transfer Mode : Auto" << std::endl;
		} else if (transferMode == Prismatic::StreamingMode::SingleXfer){
			std::cout << "Data Transfer : Single Transfer" << std::endl;
		} else {
			std::cout << "Data Transfer : Streaming" << std::endl;
		}
#endif // PRISMATIC_ENABLE_GPU
	}

	template <class T>
	bool Metadata<T>::operator==(const Metadata<T> other){
		if(interpolationFactorY != other.interpolationFactorY)return false;
		if(interpolationFactorX != other.interpolationFactorX)return false;
		if(filenameAtoms != other.filenameAtoms)return false;
		if(filenameOutput != other.filenameOutput)return false;
		if(realspacePixelSize[0] != realspacePixelSize[0])return false;
		if(realspacePixelSize[1] != other.realspacePixelSize[1])return false;
		if(potBound != other.potBound)return false;
		if(numFP != other.numFP)return false;
		if(sliceThickness != other.sliceThickness)return false;
		if(cellDim[0] != other.cellDim[0])return false;
		if(cellDim[1] != other.cellDim[1])return false;
		if(cellDim[2] != other.cellDim[2])return false;
		if(tileX != other.tileX)return false;
		if(tileY != other.tileY)return false;
		if(tileZ != other.tileZ)return false;
		if(E0 != other.E0)return false;
		if(alphaBeamMax != other.alphaBeamMax)return false;
		if(probeStepX != other.probeStepX)return false;
		if(probeStepY != other.probeStepY)return false;
		if(probeSemiangle != other.probeSemiangle)return false;
		if(C3 != other.C3)return false;
		if(C5 != other.C5)return false;
		if(probeDefocus != other.probeDefocus)return false;
		if(detectorAngleStep != other.detectorAngleStep)return false;
		if(probeXtilt != other.probeXtilt)return false;
		if(probeYtilt != other.probeYtilt)return false;
		if(scanWindowXMin != other.scanWindowXMin)return false;
		if(scanWindowXMax != other.scanWindowXMax)return false;
		if(scanWindowYMin != other.scanWindowYMin)return false;
		if(scanWindowYMax != other.scanWindowYMax)return false;
		if(randomSeed != other.randomSeed)return false;
		if(includeThermalEffects != other.includeThermalEffects)return false;
		if(includeOccupancy != other.includeOccupancy)return false;
		if(alsoDoCPUWork != other.alsoDoCPUWork)return false;
		if(save2DOutput != other.save2DOutput)return false;
		if(save3DOutput != other.save3DOutput)return false;
		if(save4DOutput != other.save4DOutput)return false;
		if(userSpecifiedCelldims != other.userSpecifiedCelldims)return false;
		return true;
	}

}
#endif //PRISMATIC_META_H
