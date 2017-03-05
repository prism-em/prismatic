//
// Created by AJ Pryor on 3/5/17.
//
#include "configure.h"
#include <iostream>
namespace PRISM {
	entry_func execute_plan;
	template <>
	void configure(Metadata<double>& meta) {
		if (meta.algorithm == 0) {
			std::cout << "Execution plan: PRISM w/ single FP configuration" << std::endl;
			execute_plan = PRISM_entry<double>;
//		execute_plan = PRISM_entry<float>;
		} else{
			std::cout << "Execution plan: Multislice w/ single FP configuration" << std::endl;
			execute_plan = Multislice_entry<double>;
		}
	}

//	template void configure(Metadata<double>&);
//	template void configure(Metadata<float>&);
}