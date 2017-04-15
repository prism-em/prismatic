//
// Created by AJ Pryor on 3/5/17.
//

// This file is for configuring the behavior of PRISM such as setting
// algorithm, CPU/GPU configuration to use, output formatting, etc.

#ifndef PRISM_CONFIGURE_H
#define PRISM_CONFIGURE_H

#include "defines.h"
#include <complex>
#include "meta.h"
#include "ArrayND.h"
#include "params.h"

#ifdef PRISM_ENABLE_GPU
#define CUDA_API_PER_THREAD_DEFAULT_STREAM
#include <cuda_runtime.h>
#endif //PRISM_ENABLE_GPU

#ifdef PRISM_BUILDING_GPU
#include "prism_progressbar.h"
#endif

namespace PRISM {

	using entry_func     = Parameters<PRISM_FLOAT_PRECISION>  (*)(Metadata<PRISM_FLOAT_PRECISION>&);
	using ms_output_func = void (*)(Parameters<PRISM_FLOAT_PRECISION>&);

	using prism_output_func = void (*)(Parameters<PRISM_FLOAT_PRECISION>&);

	using format_output_func = void (*)( Parameters<PRISM_FLOAT_PRECISION>&,
	                                     Array2D< std::complex<PRISM_FLOAT_PRECISION> >&,
	                                     const Array2D<PRISM_FLOAT_PRECISION>&,
	                                     const size_t,
	                                     const size_t);
	using fill_Scompact_func = void(*)(Parameters<PRISM_FLOAT_PRECISION> &pars);

	extern entry_func execute_plan;
	extern ms_output_func buildMultisliceOutput;
	extern prism_output_func buildPRISMOutput;
	extern format_output_func formatOutput_CPU;
	extern fill_Scompact_func fill_Scompact;
#ifdef PRISM_ENABLE_GPU
	using format_output_func_GPU = void (*)(Parameters<PRISM_FLOAT_PRECISION>&,
	                                        PRISM_FLOAT_PRECISION *,
	                                        const PRISM_FLOAT_PRECISION *,
	                                        PRISM_FLOAT_PRECISION *,
	                                        PRISM_FLOAT_PRECISION*,
	                                        const size_t,
	                                        const size_t,
	                                        const size_t&,
	                                        const size_t&,
	                                        const cudaStream_t&,
	                                        const long&);
	extern format_output_func_GPU formatOutput_GPU;

	template <class T>
	StreamingMode transferMethodAutoChooser(PRISM::Metadata<T>& meta);
#endif //PRISM_ENABLE_GPU
	void configure(Metadata<PRISM_FLOAT_PRECISION>&);
}

#endif //PRISM_CONFIGURE_H
