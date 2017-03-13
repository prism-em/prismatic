#ifndef PRISM_MULTISLICE_CUH
#define PRISM_MULTISLICE_CUH
#define CUDA_API_PER_THREAD_DEFAULT_STREAM
#include <cuda_runtime.h>
#include "params.h"
#include "configure.h"
#include <complex>
namespace PRISM {

    void buildMultisliceOutput_GPU(Parameters <PRISM_FLOAT_PRECISION> &pars,
                                   Array3D <std::complex<PRISM_FLOAT_PRECISION>> &trans,
                                   Array2D <std::complex<PRISM_FLOAT_PRECISION>> &PsiProbeInit,
                                   Array2D <PRISM_FLOAT_PRECISION> &alphaInd);
	void formatOutput_GPU_integrate(Parameters<PRISM_FLOAT_PRECISION>& pars,
                                    PRISM_FLOAT_PRECISION * psi,
                                    const PRISM_FLOAT_PRECISION * alphaInd,
                                    PRISM_FLOAT_PRECISION * stack_ph,
                                    PRISM_FLOAT_PRECISION * integratedOutput,
                                    const size_t& ay,
                                    const size_t& ax,
	                                const size_t&,
	                                const size_t&,
	                                cudaStream_t& stream);
}
#endif //PRISM_MULTISLICE_CUH
