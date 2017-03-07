//
// Created by AJ Pryor on 3/5/17.
//
#include "configure.h"
#include <iostream>
namespace PRISM {
	entry_func execute_plan;
//	ms_output_func<double> buildMultisliceOutput;
	void configure(Metadata<PRISM_FLOAT_PRECISION>& meta) {
		if (meta.algorithm == 0) {
			std::cout << "Execution plan: PRISM w/ single FP configuration" << std::endl;
			execute_plan = PRISM_entry;
//		execute_plan = PRISM_entry<float>;
		} else{
			std::cout << "Execution plan: Multislice w/ single FP configuration" << std::endl;
//			execute_plan = Multislice_entry;
//			buildMultisliceOutput = buildMultisliceOutput_cpuOnly<double>;
		}
	}

//	template void configure(Metadata<double>&);
//	template void configure(Metadata<float>&);
}