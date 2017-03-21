#ifndef PRISM_MULTISLICE_CUH
#define PRISM_MULTISLICE_CUH
#define CUDA_API_PER_THREAD_DEFAULT_STREAM
#include <cuda_runtime.h>
#include "params.h"
#include "configure.h"
#include <complex>
namespace PRISM {

    void buildMultisliceOutput_GPU_singlexfer(Parameters <PRISM_FLOAT_PRECISION> &pars,
                                              Array3D <std::complex<PRISM_FLOAT_PRECISION>> &trans,
                                              Array2D <std::complex<PRISM_FLOAT_PRECISION>> &PsiProbeInit,
                                              Array2D <PRISM_FLOAT_PRECISION> &alphaInd);

	void buildMultisliceOutput_GPU_streaming(Parameters <PRISM_FLOAT_PRECISION> &pars,
	                                         Array3D <std::complex<PRISM_FLOAT_PRECISION>> &trans,
	                                         Array2D <std::complex<PRISM_FLOAT_PRECISION>> &PsiProbeInit,
	                                         Array2D <PRISM_FLOAT_PRECISION> &alphaInd);

	void getMultisliceProbe_GPU_singlexfer(Parameters<PRISM_FLOAT_PRECISION>& pars,
	                                       PRISM_CUDA_COMPLEX_FLOAT* trans_d,
	                                       PRISM_CUDA_COMPLEX_FLOAT* PsiProbeInit_d,
	                                       PRISM_CUDA_COMPLEX_FLOAT* psi_ds,
	                                       PRISM_FLOAT_PRECISION* output_ph,
	                                       PRISM_FLOAT_PRECISION* psi_intensity_ds,
	                                       PRISM_FLOAT_PRECISION* integratedOutput_ds,
	                                       const PRISM_FLOAT_PRECISION* qya_d,
	                                       const PRISM_FLOAT_PRECISION* qxa_d,
	                                       const PRISM_CUDA_COMPLEX_FLOAT* prop_d,
	                                       const size_t& ay,
	                                       const size_t& ax,
	                                       const size_t dimj,
	                                       const size_t dimi,
	                                       const PRISM_FLOAT_PRECISION* alphaInd_d,
	                                       const cufftHandle& plan,
	                                       cudaStream_t& stream);

	void getMultisliceProbe_GPU_streaming(Parameters<PRISM_FLOAT_PRECISION>& pars,
	                                      PRISM_CUDA_COMPLEX_FLOAT* trans_d,
	                                      const std::complex<PRISM_FLOAT_PRECISION>* trans_ph,
	                                      PRISM_CUDA_COMPLEX_FLOAT* PsiProbeInit_d,
	                                      PRISM_CUDA_COMPLEX_FLOAT* psi_ds,
	                                      PRISM_FLOAT_PRECISION* output_ph,
	                                      PRISM_FLOAT_PRECISION* psi_intensity_ds,
	                                      PRISM_FLOAT_PRECISION* integratedOutput_ds,
	                                      const PRISM_FLOAT_PRECISION* qya_d,
	                                      const PRISM_FLOAT_PRECISION* qxa_d,
	                                      const PRISM_CUDA_COMPLEX_FLOAT* prop_d,
	                                      const size_t& ay,
	                                      const size_t& ax,
	                                      const size_t dimj,
	                                      const size_t dimi,
	                                      const PRISM_FLOAT_PRECISION* alphaInd_d,
	                                      const cufftHandle& plan,
	                                      cudaStream_t& stream);
}
#endif //PRISM_MULTISLICE_CUH
