//
// Created by AJ Pryor on 3/5/17.
//

// This file is for configuring the behavior of PRISM such as setting
// algorithm, CPU/GPU configuration to use, output formatting, etc.

#ifndef PRISM_CONFIGURE_H
#define PRISM_CONFIGURE_H

//#define PRISM_ENABLE_DOUBLE_PRECISION
#ifdef PRISM_ENABLE_DOUBLE_PRECISION
typedef double PRISM_FLOAT_PRECISION;
#define PRISM_FFTW_PLAN fftw_plan
#define PRISM_FFTW_PLAN_DFT_2D fftw_plan_dft_2d
#define PRISM_FFTW_EXECUTE fftw_execute
#define PRISM_FFTW_DESTROY_PLAN fftw_destroy_plan
#define PRISM_FFTW_COMPLEX fftw_complex
#else
typedef float PRISM_FLOAT_PRECISION;
#define PRISM_FFTW_PLAN fftwf_plan
#define PRISM_FFTW_PLAN_DFT_2D fftwf_plan_dft_2d
#define PRISM_FFTW_EXECUTE fftwf_execute
#define PRISM_FFTW_DESTROY_PLAN fftwf_destroy_plan
#define PRISM_FFTW_COMPLEX fftwf_complex
#endif //PRISM_ENABLE_DOUBLE_PRECISION

#include <complex>
#include "meta.h"
#include "ArrayND.h"
#include "params.h"
//#include "PRISM_entry.h"
//#include "Multislice_entry.h"
//#include "Multislice.h"
namespace PRISM {
	using entry_func     = int  (*)(Metadata<PRISM_FLOAT_PRECISION>&);
	using ms_output_func = void (*)(Parameters<PRISM_FLOAT_PRECISION>&,
                                    Array3D<std::complex<PRISM_FLOAT_PRECISION> >&,
                                    Array2D<std::complex<PRISM_FLOAT_PRECISION> >&,
                                    Array2D<PRISM_FLOAT_PRECISION>&);

	extern entry_func execute_plan;
	extern ms_output_func buildMultisliceOutput;

	template <class T>
	using Array1D = PRISM::ArrayND<1, std::vector<T> >;
	template <class T>
	using Array2D = PRISM::ArrayND<2, std::vector<T> >;
	template <class T>
	using Array3D = PRISM::ArrayND<3, std::vector<T> >;
	template <class T>
	using Array4D = PRISM::ArrayND<4, std::vector<T> >;

	void configure(Metadata<PRISM_FLOAT_PRECISION>&);





}

#endif //PRISM_CONFIGURE_H
