//
// Created by AJ Pryor on 3/5/17.
//
#include "configure.h"
#include <iostream>
namespace PRISM {
	entry_func execute_plan;
	void configure() {
		std::cout << "Execution plan: PRISM w/ single FP configuration" << std::endl;
		execute_plan = PRISM_entry<double>;
	}
}