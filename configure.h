//
// Created by AJ Pryor on 3/5/17.
//

// This file is for configuring the behavior of PRISM such as setting
// algorithm, CPU/GPU configuration to use, output formatting, etc.

#ifndef PRISM_CONFIGURE_H
#define PRISM_CONFIGURE_H
#include "meta.h"
#include "PRISM_entry.h"

namespace PRISM {
	void configure();

	using entry_func = int (*)(Metadata<double>&);
//	using entry_func = int (*)(Metadata<float>&);
	extern entry_func execute_plan;
}

#endif //PRISM_CONFIGURE_H
