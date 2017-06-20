#ifndef PRISMATIC_PARAMS_CUH
#define PRISMATIC_PARAMS_CUH
#include "defines.h"
#include "configure.h"
namespace Prismatic {
	template<class T>
	class CudaParameters {
	public:
		PRISMATIC_CUDA_COMPLEX_FLOAT **permuted_Scompact_d;
		PRISMATIC_CUDA_COMPLEX_FLOAT **PsiProbeInit_d;
		PRISMATIC_FLOAT_PRECISION **qxaReduce_d;
		PRISMATIC_FLOAT_PRECISION **qyaReduce_d;
		size_t **yBeams_d;
		size_t **xBeams_d;
		PRISMATIC_FLOAT_PRECISION **alphaInd_d;
		PRISMATIC_CUDA_COMPLEX_FLOAT **psi_ds;
		PRISMATIC_CUDA_COMPLEX_FLOAT **psi_small_ds;
		PRISMATIC_CUDA_COMPLEX_FLOAT **phaseCoeffs_ds;
		PRISMATIC_FLOAT_PRECISION **psi_intensity_ds;
		long **y_ds;
		long **x_ds;
		PRISMATIC_FLOAT_PRECISION **integratedOutput_ds;
		cufftHandle *cufft_plans, *cufft_plans_small;
		cudaStream_t *streams;
		std::complex<PRISMATIC_FLOAT_PRECISION>** Scompact_slice_ph;
		std::complex<PRISMATIC_FLOAT_PRECISION> *trans_ph;
		std::complex<PRISMATIC_FLOAT_PRECISION> *prop_ph;
		size_t **qxInd_ph;
		size_t **qyInd_ph;
		size_t **beamsIndex_ph;
        PRISMATIC_CUDA_COMPLEX_FLOAT **trans_d;
        PRISMATIC_CUDA_COMPLEX_FLOAT **prop_d;
        size_t **qxInd_d;
        size_t **qyInd_d;
        size_t **beamsIndex_d;

		// pinned memory buffers on host
		PRISMATIC_FLOAT_PRECISION               **output_ph;
		std::complex<PRISMATIC_FLOAT_PRECISION> *permuted_Scompact_ph;
		std::complex<PRISMATIC_FLOAT_PRECISION> *PsiProbeInit_ph;
		PRISMATIC_FLOAT_PRECISION               *qxaReduce_ph;
		PRISMATIC_FLOAT_PRECISION               **qxa_ph;
		PRISMATIC_FLOAT_PRECISION               **qxa_d;
		PRISMATIC_FLOAT_PRECISION               *qyaReduce_ph;
		PRISMATIC_FLOAT_PRECISION               **qya_ph;
		PRISMATIC_FLOAT_PRECISION               **qya_d;
		PRISMATIC_FLOAT_PRECISION               *alphaInd_ph;
		size_t                              *xBeams_ph;
		size_t                              *yBeams_ph;

	};
}
#endif //PRISMATIC_PARAMS_CUH