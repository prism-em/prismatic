#include <Python.h>
#include "params.h"
#include "configure.h"
#include "parseInput.h"
#ifdef PRISMATIC_ENABLE_GPU
#include "cuprismatic.h"
#endif //PRISMATIC_ENABLE_GPU


static PyObject* pyprismatic_core_go(PyObject *self, PyObject *args){
	Prismatic::Metadata<PRISMATIC_FLOAT_PRECISION> meta;
	int  interpolationFactorY = 1;
	int  interpolationFactorX = 1;
	int numFP, batchSizeTargetCPU, batchSizeTargetGPU, 
	tileX, tileY, tileZ, 
	numGPUs, numStreamsPerGPU, numThreads,includeThermalEffects, alsoDoCPUWork, save2DOutput,
	save3DOutput, save4DOutput;
	char *filenameAtoms, *filenameOutput, *algorithm, *transferMode;
	double realspacePixelSizeX, realspacePixelSizeY, potBound,
	sliceThickness, probeStepX, probeStepY,
	cellDimX, cellDimY, cellDimZ, earlyCPUStopCount, E0, alphaBeamMax,
	detectorAngleStep, probeDefocus, C3,
	C5, probeSemiangle, probeXtilt,
	probeYtilt, scanWindowXMin, scanWindowXMax,
	scanWindowYMin, scanWindowYMax, randomSeed, 
	integrationAngleMin, integrationAngleMax;
	#ifdef PRISMATIC_ENABLE_GPU
	std::cout <<"COMPILED FOR GPU" << std::endl;
	#endif //PRISMATIC_ENABLE_GPU

	if (!PyArg_ParseTuple(
		args, "iissdddiddddiiiddiiiiiddddddddddddddispppppdds",
	    &interpolationFactorX,
	    &interpolationFactorY,
	    &filenameAtoms,
	    &filenameOutput,
	    &realspacePixelSizeX,
	    &realspacePixelSizeY, 
	    &potBound,
	    &numFP,
	    &sliceThickness,
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
		&C3,
		&C5,
		&probeSemiangle,
		&detectorAngleStep,
		&probeXtilt,
		&probeYtilt,
		&scanWindowXMin,
		&scanWindowXMax,
		&scanWindowYMin,
		&scanWindowYMax,
		&randomSeed,	
		&algorithm,
		&includeThermalEffects,
		&alsoDoCPUWork,
		&save2DOutput,
		&save3DOutput,
		&save4DOutput,
		&integrationAngleMin,
		&integrationAngleMax,
		&transferMode)){
		return NULL;
	} 
	meta.interpolationFactorX 	 = interpolationFactorX;
	meta.interpolationFactorY 	 = interpolationFactorY;
	meta.filename_atoms       	 = filenameAtoms;
	meta.filename_output      	 = filenameOutput;
	meta.realspace_pixelSize[0]  = realspacePixelSizeY;
	meta.realspace_pixelSize[1]  = realspacePixelSizeX;
	meta.potBound       	 	 = potBound;
	meta.numFP      			 = numFP;
	meta.sliceThickness      	 = sliceThickness;
	meta.cellDim[2] 		  	 = cellDimX;
	meta.cellDim[1]  			 = cellDimY;
	meta.cellDim[0]  			 = cellDimZ;
	meta.tileX  				 = tileX;
	meta.tileY  				 = tileY;
	meta.tileZ  			     = tileZ;
	meta.E0  	     			 = E0; 
	meta.alphaBeamMax 			 = alphaBeamMax;
	meta.NUM_GPUS                = numGPUs;
    meta.NUM_STREAMS_PER_GPU     = numStreamsPerGPU;
    meta.NUM_THREADS             = numThreads;
	meta.batch_size_target_CPU   = batchSizeTargetCPU;
	meta.batch_size_target_GPU   = batchSizeTargetGPU;
	meta.gpu_cpu_ratio 			 = earlyCPUStopCount;
	meta.probe_stepX      	 	 = probeStepX;
	meta.probe_stepY      	 	 = probeStepY;
	meta.probeDefocus 			 = probeDefocus;
	meta.C3 				     = C3;
	meta.C5 				     = C5;
	meta.probeSemiangle 		 = probeSemiangle;
	meta.detector_angle_step 	 = detectorAngleStep;
	meta.probeXtilt 			 = probeXtilt;
	meta.probeYtilt 			 = probeYtilt;
	meta.scanWindowXMin 		 = scanWindowXMin;
	meta.scanWindowXMax 		 = scanWindowXMax;
	meta.scanWindowYMin 	     = scanWindowYMin;
	meta.scanWindowYMax 		 = scanWindowYMax;
	meta.random_seed 		     = randomSeed;
	if (std::string(algorithm) == "multislice"){
	} else {
		meta.algorithm 			 = Prismatic::Algorithm::PRISM;
	}

	meta.include_thermal_effects = includeThermalEffects;
	meta.also_do_CPU_work        = alsoDoCPUWork;
	meta.save2DOutput			 = save2DOutput;
	meta.save3DOutput			 = save3DOutput;
	meta.save4DOutput			 = save4DOutput;
	meta.integration_angle_min   = integrationAngleMin;
	meta.integration_angle_max   = integrationAngleMax;
	
	if (std::string(transferMode) == "singlexfer"){
		meta.transfer_mode 		 = Prismatic::StreamingMode::SingleXfer;
	} else if (std::string(transferMode) == "streaming"){
		meta.transfer_mode 		 = Prismatic::StreamingMode::Stream;
	} else {
		meta.transfer_mode 		 = Prismatic::StreamingMode::Auto;
	}

	meta.algorithm = Prismatic::Algorithm::PRISM;

	// print metadata
    meta.toString();

	// configure simulation behavior
	Prismatic::configure(meta);

	// execute simulation
	Prismatic::execute_plan(meta);

	Py_RETURN_NONE;

}

static PyMethodDef pyprismatic_core_methods[] = {
	{"go",(PyCFunction)pyprismatic_core_go, METH_VARARGS, "Execute Prismatic calculation"},
   	{NULL, NULL, 0, NULL}
};


static struct PyModuleDef module_def = {
	PyModuleDef_HEAD_INIT,"pypristmatic.core","Python wrapper for Prismatic\
	 package for fast image simulation using the PRISM and multislice\
	 algorithms in Scanning Transmission Electron Microscopy (STEM)",-1,pyprismatic_core_methods
};

PyMODINIT_FUNC PyInit_core(){
	return PyModule_Create(&module_def);
}
