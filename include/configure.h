// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// Prismatic is distributed under the GNU General Public License (GPL)
// If you use Prismatic, we kindly ask that you cite the following papers:

// 1. Ophus, C.: A fast image simulation algorithm for scanning
//    transmission electron microscopy. Advanced Structural and
//    Chemical Imaging 3(1), 13 (2017)

// 2. Pryor, Jr., A., Ophus, C., and Miao, J.: A Streaming Multi-GPU
//    Implementation of Image Simulation Algorithms for Scanning
//	  Transmission Electron Microscopy. arXiv:1706.08563 (2017)

// This file is for configuring the behavior of PRISM such as setting
// algorithm, CPU/GPU configuration to use, output formatting, etc.

#ifndef PRISM_CONFIGURE_H
#define PRISM_CONFIGURE_H

#include "defines.h"
#include <complex>
#include "meta.h"
#include "ArrayND.h"
#include "params.h"

#ifdef PRISMATIC_ENABLE_GPU
//#define CUDA_API_PER_THREAD_DEFAULT_STREAM
#include <cuda_runtime.h>
#endif //PRISMATIC_ENABLE_GPU

#ifdef PRISMATIC_BUILDING_GUI
#include "prism_progressbar.h"
#endif

namespace Prismatic {
	using entry_func     = Parameters<PRISMATIC_FLOAT_PRECISION>  (*)(Metadata<PRISMATIC_FLOAT_PRECISION>&);
	using ms_output_func = void (*)(Parameters<PRISMATIC_FLOAT_PRECISION>&);

	using prism_output_func = void (*)(Parameters<PRISMATIC_FLOAT_PRECISION>&);

	using format_output_func = void (*)( Parameters<PRISMATIC_FLOAT_PRECISION>&,
	                                     Array2D< std::complex<PRISMATIC_FLOAT_PRECISION> >&,
	                                     const Array2D<PRISMATIC_FLOAT_PRECISION>&,
	                                     const size_t,
	                                     const size_t);
	using fill_Scompact_func = void(*)(Parameters<PRISMATIC_FLOAT_PRECISION> &pars);

	extern entry_func execute_plan;
	extern ms_output_func buildMultisliceOutput;
	extern prism_output_func buildPRISMOutput;
	extern format_output_func formatOutput_CPU;
	extern fill_Scompact_func fill_Scompact;
#ifdef PRISMATIC_ENABLE_GPU
	using format_output_func_GPU = void (*)(Parameters<PRISMATIC_FLOAT_PRECISION>&,
	                                        PRISMATIC_FLOAT_PRECISION *,
	                                        const PRISMATIC_FLOAT_PRECISION *,
	                                        PRISMATIC_FLOAT_PRECISION *,
	                                        PRISMATIC_FLOAT_PRECISION*,
	                                        const size_t,
	                                        const size_t,
	                                        const size_t&,
	                                        const size_t&,
	                                        const cudaStream_t&,
	                                        const long&);
	extern format_output_func_GPU formatOutput_GPU;

	template <class T>
	StreamingMode transferMethodAutoChooser(Prismatic::Metadata<T>& meta);
#endif //PRISMATIC_ENABLE_GPU
	void configure(Metadata<PRISMATIC_FLOAT_PRECISION>&);
}

#endif //PRISM_CONFIGURE_H
