#ifndef PRISM_PARAMS_CUH
#define PRISM_PARAMS_CUH
#include "defines.h"
#include "configure.h"
namespace PRISM {
	template<class T>
	class CudaParameters {
	public:
		PRISM_CUDA_COMPLEX_FLOAT **permuted_Scompact_d;
		PRISM_CUDA_COMPLEX_FLOAT **PsiProbeInit_d;
		PRISM_FLOAT_PRECISION **qxaReduce_d;
		PRISM_FLOAT_PRECISION **qyaReduce_d;
		size_t **yBeams_d;
		size_t **xBeams_d;
		PRISM_FLOAT_PRECISION **alphaInd_d;
		PRISM_CUDA_COMPLEX_FLOAT **psi_ds;
		PRISM_CUDA_COMPLEX_FLOAT **phaseCoeffs_ds;
		PRISM_FLOAT_PRECISION **psi_intensity_ds;
		long **y_ds;
		long **x_ds;
		PRISM_FLOAT_PRECISION **integratedOutput_ds;
		cufftHandle *cufft_plans;
		cudaStream_t *streams;

		// pinned memory buffers on host
		PRISM_FLOAT_PRECISION               **output_ph;
		std::complex<PRISM_FLOAT_PRECISION> *permuted_Scompact_ph;
		std::complex<PRISM_FLOAT_PRECISION> *PsiProbeInit_ph;
		PRISM_FLOAT_PRECISION               *qxaReduce_ph;
		PRISM_FLOAT_PRECISION               *qyaReduce_ph;
		PRISM_FLOAT_PRECISION               *alphaInd_ph;
		size_t                              *xBeams_ph;
		size_t                              *yBeams_ph;

	};
}
#endif //PRISM_PARAMS_CUH