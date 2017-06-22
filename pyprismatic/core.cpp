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
	int numFP,
	tileX, tileY, tileZ, 
	batch_size_target_CPU, batch_size_target_GPU, 
	NUM_GPUS, NUM_STREAMS_PER_GPU, NUM_THREADS;
	char *filename_atoms, *filename_output, *algorithm, *transfer_mode;
	double realspace_pixelSizeX, realspace_pixelSizeY, potBound,
	 sliceThickness, probe_stepX, probe_stepY,
	 cellDimX, cellDimY, cellDimZ, gpu_cpu_ratio, E0, alphaBeamMax,
	 detector_angle_step, probeDefocus, C3,
	 C5, probeSemiangle, probeXtilt,
	 probeYtilt, scanWindowXMin, scanWindowXMax,
	 scanWindowYMin, scanWindowYMax, random_seed, 
	 integration_angle_min, integration_angle_max;
	 int include_thermal_effects, also_do_CPU_work, save2DOutput,
	 save3DOutput, save4DOutput;
	std::cout << "before meta.interpolationFactorX= " << meta.interpolationFactorX<< std::endl;
	std::cout << "before meta.interpolationFactorY= " << meta.interpolationFactorY<< std::endl;
	#ifdef PRISMATIC_ENABLE_GPU
	std::cout <<"COMPILED FOR GPU" << std::endl;
	#endif //PRISMATIC_ENABLE_GPU

	// if (!PyArg_ParseTuple(
	// 	args, "ii",
	//     &interpolationFactorX,
	//     &interpolationFactorY)){
	// 	return NULL;
	// }

	if (!PyArg_ParseTuple(
		// args, "iissdddiddddiiiddiiiiiddddddddds",
		args, "iissdddiddddiiiddiiiiiddddddddddddddispppppdds",
	    &interpolationFactorX,
	    &interpolationFactorY,
	    &filename_atoms,
	    &filename_output,
	    &realspace_pixelSizeX,
	    &realspace_pixelSizeY, 
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
	    &NUM_GPUS,
	    &NUM_STREAMS_PER_GPU,
	    &NUM_THREADS,
	    &batch_size_target_CPU,
	    &batch_size_target_GPU,
	    &gpu_cpu_ratio,
	    &probe_stepX,
	    &probe_stepY, 
		&probeDefocus,
		&C3,
		&C5,
		&probeSemiangle,
		&detector_angle_step,
		&probeXtilt,
		&probeYtilt,
		&scanWindowXMin,
		&scanWindowXMax,
		&scanWindowYMin,
		&scanWindowYMax,
		&random_seed,	
		&algorithm,
		&include_thermal_effects,
		&also_do_CPU_work,
		&save2DOutput,
		&save3DOutput,
		&save4DOutput,
		&integration_angle_min,
		&integration_angle_max,
		&transfer_mode)){
		return NULL;
	} 
	meta.interpolationFactorX 	 = interpolationFactorX;
	meta.interpolationFactorY 	 = interpolationFactorY;
	std::cout << "interpolationFactorX= " << interpolationFactorX<< std::endl;
	std::cout << "interpolationFactorY= " << interpolationFactorY<< std::endl;
	std::cout << "meta.interpolationFactorX= " << meta.interpolationFactorX<< std::endl;
	meta.filename_atoms       	 = filename_atoms;
	meta.filename_output      	 = filename_output;
	meta.realspace_pixelSize[0]  = realspace_pixelSizeY;
	meta.realspace_pixelSize[1]  = realspace_pixelSizeX;
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
	meta.NUM_GPUS                = NUM_GPUS;
    meta.NUM_STREAMS_PER_GPU     = NUM_STREAMS_PER_GPU;
    meta.NUM_THREADS             = NUM_THREADS;
	meta.batch_size_target_CPU   = batch_size_target_CPU;
	meta.batch_size_target_GPU   = batch_size_target_GPU;
	meta.gpu_cpu_ratio 			 = gpu_cpu_ratio;
	meta.probe_stepX      	 	 = probe_stepX;
	meta.probe_stepY      	 	 = probe_stepY;
	meta.probeDefocus 			 = probeDefocus;
	meta.C3 				     = C3;
	meta.C5 				     = C5;
	meta.probeSemiangle 		 = probeSemiangle;
	meta.detector_angle_step 	 = detector_angle_step;
	meta.probeXtilt 			 = probeXtilt;
	meta.probeYtilt 			 = probeYtilt;
	meta.scanWindowXMin 		 = scanWindowXMin;
	meta.scanWindowXMax 		 = scanWindowXMax;
	meta.scanWindowYMin 	     = scanWindowYMin;
	meta.scanWindowYMax 		 = scanWindowYMax;
	meta.random_seed 		     = random_seed;
	if (std::string(algorithm) == "multislice"){
	} else {
		meta.algorithm 			 = Prismatic::Algorithm::PRISM;
	}

	meta.include_thermal_effects = include_thermal_effects;
	meta.also_do_CPU_work        = also_do_CPU_work;
	meta.save2DOutput			 = save2DOutput;
	meta.save3DOutput			 = save3DOutput;
	meta.save4DOutput			 = save4DOutput;
	meta.integration_angle_min   = integration_angle_min;
	meta.integration_angle_max   = integration_angle_max;
	
	if (std::string(transfer_mode) == "singlexfer"){
		meta.transfer_mode 		 = Prismatic::StreamingMode::SingleXfer;
	} else if (std::string(transfer_mode) == "streaming"){
		meta.transfer_mode 		 = Prismatic::StreamingMode::Stream;
	} else {
		meta.transfer_mode 		 = Prismatic::StreamingMode::Auto;
	}

	// // meta.filename_atoms = "/home/aj/hdd1/clion/PRISM/SI100.XYZ";
	// // meta.filename_output = "/home/aj/hdd1/clion/PRISM/output_python.mrc";
	meta.algorithm = Prismatic::Algorithm::PRISM;

	// print metadata
    meta.toString();

	// configure simulation behavior
	Prismatic::configure(meta);

	// execute simulation
	Prismatic::execute_plan(meta);

	// "interpolationFactorX",
	// "interpolationFactorY",
	// "filename_atoms",
	// "filename_output",
	// "realspace_pixelSize",
	// "potBound",
	// "numFP",
	// "sliceThickness",
	// "cellDim",
	// "tileX",
	// "tileY",
	// "tileZ",
	// "E0",
	// "alphaBeamMax",
	// "NUM_GPUS",
	// "NUM_STREAMS_PER_GPU",
	// "NUM_THREADS",
	// "batch_size_target_CPU",
	// "batch_size_target_GPU",
	// "batch_size_CPU",
	// "batch_size_GPU",
	// "gpu_cpu_ratio",
 //    "probe_stepX",
	// "probe_stepY",
	// "probeDefocus",
	// "C3",
	// "C5",
	// "probeSemiangle",
	// "detector_angle_step",
	// "probeXtilt",
	// "probeYtilt",
	// "scanWindowXMin",
	// "scanWindowXMax",
	// "scanWindowYMin",
	// "scanWindowYMax",
	// "random_seed",
	// "algorithm",
	// "include_thermal_effects",
	// "also_do_CPU_work",
	// "save2DOutput",
	// "save3DOutput",
	// "save4DOutput",
	// "integration_angle_min",
	// "integration_angle_max",
	// "transfer_mode"


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
