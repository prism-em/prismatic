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
	format_output_func formatOutput_cpu;
#ifdef PRISM_ENABLE_GPU
	format_output_func formatOutput_gpu;
#endif
	void configure(Metadata<PRISM_FLOAT_PRECISION>& meta) {
		formatOutput_cpu = formatOutput_cpu_integrate;
#ifdef PRISM_ENABLE_GPU
		format_output_gpu =formatOutput_gpu_integrate;
#endif
		if (meta.algorithm == 0) {
			std::cout << "Execution plan: PRISM w/ single FP configuration" << std::endl;
			execute_plan = PRISM_entry;
		} else {
			std::cout << "Execution plan: Multislice w/ single FP configuration" << std::endl;
			execute_plan = Multislice_entry;
#ifdef PRISM_ENABLE_GPU
			std::cout << "Using GPU codes" << std::endl;
			buildMultisliceOutput = buildMultisliceOutput_gpu;
#else
			buildMultisliceOutput = buildMultisliceOutput_cpuOnly;
#endif //PRISM_ENABLE_GPU
		}
	}

//	template void configure(Metadata<double>&);
//	template void configure(Metadata<float>&);
}