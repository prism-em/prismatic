//
// Created by AJ Pryor on 3/5/17.
//
#include "configure.h"
#include "PRISM_entry.h"
#include "Multislice_entry.h"
#include "Multislice.h"
#include <iostream>

//#define PRISM_ENABLE_GPU

#ifdef PRISM_ENABLE_GPU
#include "Multislice.cuh"
#endif //PRISM_ENABLE_GPU
namespace PRISM {
	entry_func execute_plan;
	ms_output_func buildMultisliceOutput;
	void configure(Metadata<PRISM_FLOAT_PRECISION>& meta) {
		if (meta.algorithm == 0) {
			std::cout << "Execution plan: PRISM w/ single FP configuration" << std::endl;
			execute_plan = PRISM_entry;
		} else {
			std::cout << "Execution plan: Multislice w/ single FP configuration" << std::endl;
			execute_plan = Multislice_entry;
#ifdef PRISM_ENABLE_GPU
			buildMultisliceOutput = buildMultisliceOutput_gpu;
#else
			buildMultisliceOutput = buildMultisliceOutput_cpuOnly;
#endif //PRISM_ENABLE_GPU
		}
	}

//	template void configure(Metadata<double>&);
//	template void configure(Metadata<float>&);
}