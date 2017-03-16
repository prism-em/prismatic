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

}
#endif //PRISM_MULTISLICE_CUH
