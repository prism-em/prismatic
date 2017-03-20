#ifndef PRISM_PRISM02_CUH
#define PRISM_PRISM02_CUH
#include "params.h"
#include "configure.h"
#include <complex>
namespace PRISM{
	void fill_Scompact_GPU_singlexfer(Parameters<PRISM_FLOAT_PRECISION> &pars);
	void fill_Scompact_GPU_streaming(Parameters <PRISM_FLOAT_PRECISION> &pars);
	void propagatePlaneWave_GPU_singlexfer(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                            PRISM_CUDA_COMPLEX_FLOAT* trans_d,
	                            PRISM_CUDA_COMPLEX_FLOAT* psi_d,
	                            PRISM_CUDA_COMPLEX_FLOAT* psi_small_d,
	                            std::complex<PRISM_FLOAT_PRECISION>* Scompact_slice_ph,
	                            const size_t* qyInd_d,
	                            const size_t* qxInd_d,
	                            const PRISM_CUDA_COMPLEX_FLOAT* prop_d,
	                            const size_t* beamsIndex,
	                            const size_t& beamNumber,
	                            const cufftHandle& plan,
	                            const cufftHandle& plan_small,
	                            cudaStream_t& stream);

	void propagatePlaneWave_GPU_streaming(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                                      PRISM_CUDA_COMPLEX_FLOAT* trans_d,
	                                      const std::complex<PRISM_FLOAT_PRECISION> *trans_ph,
	                                      PRISM_CUDA_COMPLEX_FLOAT* psi_d,
	                                      PRISM_CUDA_COMPLEX_FLOAT* psi_small_d,
	                                      std::complex<PRISM_FLOAT_PRECISION>* Scompact_slice_ph,
	                                      const size_t* qyInd_d,
	                                      const size_t* qxInd_d,
	                                      const PRISM_CUDA_COMPLEX_FLOAT* prop_d,
	                                      const size_t* beamsIndex,
	                                      const size_t& beamNumber,
	                                      const cufftHandle& plan,
	                                      const cufftHandle& plan_small,
	                                      cudaStream_t& stream);
}
#endif //PRISM_PRISM02_CUH