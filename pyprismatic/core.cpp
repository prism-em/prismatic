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

#include <Python.h>
#include "params.h"
#include "configure.h"
#include "parseInput.h"
#include "go.h"
#include <iostream>
#include <stdio.h>
#include "aberration.h"
#include "probe.h"
#ifdef PRISMATIC_ENABLE_GPU
#include "cuprismatic.h"
#endif //PRISMATIC_ENABLE_GPU

static PyObject *pyprismatic_core_go(PyObject *self, PyObject *args)
{
	Prismatic::Metadata<PRISMATIC_FLOAT_PRECISION> meta;
	int interpolationFactorY = 1;
	int interpolationFactorX = 1;
	int randomSeed;
	int numFP, batchSizeTargetCPU, batchSizeTargetGPU,
		tileX, tileY, tileZ,
		numGPUs, numStreamsPerGPU, numThreads, includeThermalEffects, includeOccupancy, alsoDoCPUWork,
		save2DOutput, save3DOutput, save4DOutput, saveDPC_CoM, savePotentialSlices, nyquistSampling, crop4DOutput,
		zSampling, numSlices, potential3D, saveSMatrix, importPotential, importSMatrix, saveComplexOutputWave, saveProbe, matrixRefocus,
        earlyCPUStopCount;
    unsigned long long int maxFileSize;
	char *filenameAtoms, *filenameOutput, *algorithm, *transferMode,
		 *aberrations_file, *probes_file, *importFile, *importPath;
	double realspacePixelSizeX, realspacePixelSizeY, potBound,
		sliceThickness, probeStepX, probeStepY,
		cellDimX, cellDimY, cellDimZ, E0, alphaBeamMax,
		detectorAngleStep, probeDefocus, C3,
		C5, probeSemiangle, probeXtilt,
		probeYtilt, scanWindowXMin, scanWindowXMax,
		scanWindowYMin, scanWindowYMax,
		integrationAngleMin, integrationAngleMax, crop4Damax, zStart,
		scanWindowXMin_r, scanWindowXMax_r, scanWindowYMin_r, scanWindowYMax_r,
		probeDefocus_min, probeDefocus_max, probeDefocus_step, probeDefocus_sigma,
		minXtilt, maxXtilt, minYtilt, maxYtilt, minRtilt, maxRtilt, xTiltOffset, yTiltOffset, xTiltStep, yTiltStep;
#ifdef PRISMATIC_ENABLE_GPU
	std::cout << "COMPILED FOR GPU" << std::endl;
#endif //PRISMATIC_ENABLE_GPU
	//i - integer
	//s - string
	//d - double
	//p - bool
	if (!PyArg_ParseTuple( 
			args, "iissdddidiiddddiiiddiiiiiidddddddddsddddddddddddddddddddddsispppppppddsppppdpppippssK",
			&interpolationFactorX,
			&interpolationFactorY,
			&filenameAtoms,
			&filenameOutput,
			&realspacePixelSizeX,
			&realspacePixelSizeY,
			&potBound,
			&numFP,
			&sliceThickness,
			&zSampling,
			&numSlices,
			&zStart,
			&cellDimX,
			&cellDimY,
			&cellDimZ,
			&tileX,
			&tileY,
			&tileZ,
			&E0,
			&alphaBeamMax,
			&numGPUs,
			&numStreamsPerGPU,
			&numThreads,
			&batchSizeTargetCPU,
			&batchSizeTargetGPU,
			&earlyCPUStopCount,
			&probeStepX,
			&probeStepY,
			&probeDefocus,
			&probeDefocus_min,
			&probeDefocus_max,
			&probeDefocus_step,
			&probeDefocus_sigma,
			&C3,
			&C5,
			&aberrations_file, //pass a string to read a file for now; file is created on python side as scratch
			&probeSemiangle,
			&detectorAngleStep,
			&probeXtilt,
			&probeYtilt,
			&minXtilt,
			&maxXtilt,
			&minYtilt,
			&maxYtilt,
			&minRtilt,
			&maxRtilt,
			&xTiltOffset,
			&yTiltOffset,
			&xTiltStep,
			&yTiltStep,
			&scanWindowXMin,
			&scanWindowXMax,
			&scanWindowYMin,
			&scanWindowYMax,
			&scanWindowXMin_r,
			&scanWindowXMax_r,
			&scanWindowYMin_r,
			&scanWindowYMax_r,
			&probes_file,
			&randomSeed,
			&algorithm,
			&potential3D,
			&includeThermalEffects,
			&includeOccupancy,
			&alsoDoCPUWork,
			&save2DOutput,
			&save3DOutput,
			&save4DOutput,
			&integrationAngleMin,
			&integrationAngleMax,
			&transferMode,
			&saveDPC_CoM,
			&savePotentialSlices,
			&saveSMatrix,
			&crop4DOutput,
			&crop4Damax,
			&nyquistSampling,
			&importPotential,
			&importSMatrix,
			&saveProbe,
			&saveComplexOutputWave,
			&matrixRefocus,
			&importFile,
			&importPath,
			&maxFileSize))
	{
		return NULL;
	}
	meta.interpolationFactorX = interpolationFactorX;
	meta.interpolationFactorY = interpolationFactorY;
	meta.filenameAtoms = filenameAtoms;
	meta.filenameOutput = filenameOutput;
	meta.realspacePixelSize[0] = realspacePixelSizeY;
	meta.realspacePixelSize[1] = realspacePixelSizeX;
	meta.potBound = potBound;
	meta.numFP = numFP;
	meta.sliceThickness = sliceThickness;
	meta.zSampling = zSampling;
	meta.numSlices = numSlices;
	meta.zStart = zStart;
	meta.cellDim[2] = cellDimX;
	meta.cellDim[1] = cellDimY;
	meta.cellDim[0] = cellDimZ;
	meta.tileX = tileX;
	meta.tileY = tileY;
	meta.tileZ = tileZ;
	meta.E0 = E0 * 1000;
	meta.alphaBeamMax = alphaBeamMax / 1000;
	meta.numGPUs = numGPUs;
	meta.numStreamsPerGPU = numStreamsPerGPU;
	meta.numThreads = numThreads;
	meta.batchSizeTargetCPU = batchSizeTargetCPU;
	meta.batchSizeTargetGPU = batchSizeTargetGPU;
	meta.earlyCPUStopCount = earlyCPUStopCount;
	meta.probeStepX = probeStepX;
	meta.probeStepY = probeStepY;
	meta.probeDefocus = probeDefocus;
	
	//set up sim series if values require it
	if(probeDefocus_min != 0.0 && probeDefocus_max > probeDefocus_min && probeDefocus_step != 0)
	{
		meta.probeDefocus_min = probeDefocus_min;
		meta.probeDefocus_max = probeDefocus_max;
		meta.probeDefocus_step = probeDefocus_step;
		meta.simSeries = true;
	}

	if(probeDefocus_sigma != 0.0)
	{
		meta.probeDefocus_sigma = probeDefocus_sigma;
		meta.simSeries = true;
	}
	
	meta.C3 = C3;
	meta.C5 = C5;
	if(std::string(aberrations_file).size() > 0)
	{
		meta.arbitraryAberrations = true;
		meta.aberrations = Prismatic::readAberrations(std::string(aberrations_file));
	}
	meta.probeSemiangle = probeSemiangle / 1000;
	meta.detectorAngleStep = detectorAngleStep / 1000;
	meta.probeXtilt = probeXtilt / 1000;
	meta.probeYtilt = probeYtilt / 1000;
	meta.minXtilt = minXtilt / 1000;
	meta.maxXtilt = maxXtilt / 1000;
	meta.xTiltStep = xTiltStep /1000;
	meta.xTiltOffset = xTiltOffset / 1000;
	meta.minYtilt = minYtilt / 1000;
	meta.maxYtilt = maxYtilt / 1000;
	meta.yTiltStep = yTiltStep /1000;
	meta.yTiltOffset = yTiltOffset / 1000;
	meta.minRtilt = minRtilt / 1000;
	meta.maxRtilt = maxRtilt / 1000;
	meta.scanWindowXMin = scanWindowXMin;
	meta.scanWindowXMax = scanWindowXMax;
	meta.scanWindowYMin = scanWindowYMin;
	meta.scanWindowYMax = scanWindowYMax;
	meta.scanWindowXMin_r = scanWindowXMin_r;
	meta.scanWindowXMax_r = scanWindowXMax_r;
	meta.scanWindowYMin_r = scanWindowYMin_r;
	meta.scanWindowYMax_r = scanWindowYMax_r;
	if(scanWindowXMax_r) meta.realSpaceWindow_x = true;
	if(scanWindowYMax_r) meta.realSpaceWindow_y = true;

	if(std::string(probes_file).size() > 0)
	{
		meta.arbitraryProbes = true;
		std::cout << "Reading probe positions from" << " " << std::string(probes_file) << std::endl;
		Prismatic::readProbes(std::string(probes_file), meta.probes_x, meta.probes_y);
	}
	meta.randomSeed = randomSeed;
	if (std::string(algorithm) == "multislice" || std::string(algorithm) == "m")
	{
		meta.algorithm = Prismatic::Algorithm::Multislice;
	}
	else if(std::string(algorithm) == "prism" || std::string(algorithm) == "p")
	{
		meta.algorithm = Prismatic::Algorithm::PRISM;
	}
	else if(std::string(algorithm) == "hrtem" || std::string(algorithm) == "t")
	{
		meta.algorithm = Prismatic::Algorithm::HRTEM;
	}
	meta.potential3D = potential3D;
	meta.includeThermalEffects = includeThermalEffects;
	meta.includeOccupancy = includeOccupancy;
	meta.alsoDoCPUWork = alsoDoCPUWork;
	meta.save2DOutput = save2DOutput;
	meta.save3DOutput = save3DOutput;
	meta.save4DOutput = save4DOutput;
	meta.crop4DOutput = crop4DOutput;
	meta.crop4Damax = crop4Damax / 1000;
	meta.saveDPC_CoM = saveDPC_CoM;
	meta.savePotentialSlices = savePotentialSlices;
	meta.saveSMatrix = saveSMatrix;
	meta.integrationAngleMin = integrationAngleMin / 1000;
	meta.integrationAngleMax = integrationAngleMax / 1000;
	meta.nyquistSampling = nyquistSampling;
	meta.importPotential = importPotential;
	meta.importSMatrix = importSMatrix;
	meta.saveComplexOutputWave = saveComplexOutputWave;
	meta.saveProbe = saveProbe;
    meta.saveProbeComplex = (saveProbe == 1) ? false : true;
	meta.maxFileSize = maxFileSize;
	meta.matrixRefocus = matrixRefocus;
	meta.importFile = std::string(importFile);
	meta.importPath = std::string(importPath);

	if (std::string(transferMode) == "singlexfer")
	{
		meta.transferMode = Prismatic::StreamingMode::SingleXfer;
	}
	else if (std::string(transferMode) == "streaming")
	{
		meta.transferMode = Prismatic::StreamingMode::Stream;
	}
	else
	{
		meta.transferMode = Prismatic::StreamingMode::Auto;
	}

	// print metadata
	//meta.toString();
	int scratch = Prismatic::writeParamFile(meta,"scratch_param.txt");

	Prismatic::Metadata<PRISMATIC_FLOAT_PRECISION> tmp_meta;
	if(Prismatic::parseParamFile(tmp_meta,"scratch_param.txt"))
	{
		Prismatic::go(meta);
	}else{
		std::cout << "Invalid parameters detected. Cancelling calculation, please check inputs." << std::endl;
	}

	Py_RETURN_NONE;
}

static PyMethodDef pyprismatic_core_methods[] = {
	{"go", (PyCFunction)pyprismatic_core_go, METH_VARARGS, "Execute Prismatic calculation"},
	{NULL, NULL, 0, NULL}};

static struct PyModuleDef module_def = {
	PyModuleDef_HEAD_INIT, "pypristmatic.core", "Python wrapper for Prismatic\
	 package for fast image simulation using the PRISM and multislice\
	 algorithms in Scanning Transmission Electron Microscopy (STEM)",
	-1, pyprismatic_core_methods};

PyMODINIT_FUNC PyInit_core()
{
	return PyModule_Create(&module_def);
}
