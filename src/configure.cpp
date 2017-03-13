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
#include "Multislice_entry.h"
#endif //PRISM_ENABLE_GPU
namespace PRISM {
	entry_func execute_plan;
	ms_output_func buildMultisliceOutput;
	format_output_func formatOutput_CPU;
	fill_Scompact_func fill_Scompact;
#ifdef PRISM_ENABLE_GPU
	format_output_func_GPU formatOutput_GPU;
#endif
	void configure(Metadata<PRISM_FLOAT_PRECISION>& meta) {
		formatOutput_CPU = formatOutput_CPU_integrate;
#ifdef PRISM_ENABLE_GPU
		formatOutput_GPU = formatOutput_GPU_integrate;
#endif
		if (meta.algorithm == Algorithm::PRISM) {
			std::cout << "Execution plan: PRISM w/ single FP configuration" << std::endl;
			execute_plan = PRISM_entry;
#ifdef PRISM_ENABLE_GPU
#else
			fill_Scompact = fill_Scompact_CPUOnly;
#endif //PRISM_ENABLE_GPU
		} else if (meta.algorithm == Algorithm::Multislice) {
			std::cout << "Execution plan: Multislice w/ single FP configuration" << std::endl;
			execute_plan = Multislice_entry;
#ifdef PRISM_ENABLE_GPU
			std::cout << "Using GPU codes" << std::endl;
			buildMultisliceOutput = buildMultisliceOutput_GPU;
#else
			buildMultisliceOutput = buildMultisliceOutput_CPUOnly;
#endif //PRISM_ENABLE_GPU
		}
	}
}