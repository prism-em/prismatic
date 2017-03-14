#include "utility.cuh"
#define PI 3.14159265359
// define some constants
__device__ __constant__ PRISM_FLOAT_PRECISION pi       = PI;
__device__ __constant__ PRISM_CUDA_COMPLEX_FLOAT i     = {0, 1};
__device__ __constant__ PRISM_CUDA_COMPLEX_FLOAT pi_cx = {PI, 0};
__device__ __constant__ PRISM_CUDA_COMPLEX_FLOAT minus_2pii = {0, -2*PI};

// computes exp(real(a) + i * imag(a))
__device__ __forceinline__ cuDoubleComplex exp_cx(const cuDoubleComplex a){
	double e = exp(a.x);
	double s,c;
	sincos(a.y, &s, &c);
	return make_cuDoubleComplex(e*c, e*s);
}
__device__ __forceinline__ cuFloatComplex exp_cx(const cuFloatComplex a){
	float e = expf(a.x);
	float s,c;
	sincosf(a.y, &s, &c);
	return make_cuFloatComplex(e*c, e*s);
}

// creates initial probe using existing GPU memory rather than streaming each probe
__global__ void initializePsi_oneNonzero(cuFloatComplex *psi_d, const size_t N, const size_t beamLoc){
	int idx = threadIdx.x + blockDim.x*blockIdx.x;
	if (idx < N) {
		psi_d[idx] = (idx == beamLoc) ? make_cuFloatComplex(1,0):make_cuFloatComplex(0,0);
	}
}

__global__ void initializePsi_oneNonzero(cuDoubleComplex *psi_d, const size_t N, const size_t beamLoc){
	int idx = threadIdx.x + blockDim.x*blockIdx.x;
	if (idx < N) {
		psi_d[idx] = (idx == beamLoc) ? make_cuDoubleComplex(1,0):make_cuDoubleComplex(0,0);
	}
}

// multiply two complex arrays
__global__ void multiply_inplace(PRISM_CUDA_COMPLEX_FLOAT* arr,
                                 const PRISM_CUDA_COMPLEX_FLOAT* other,
                                 const size_t N){
	int idx = threadIdx.x + blockDim.x*blockIdx.x;
	if (idx < N) {
		PRISM_CUDA_COMPLEX_FLOAT a = arr[idx];
		PRISM_CUDA_COMPLEX_FLOAT o = other[idx];
		arr[idx].x = a.x * o.x - a.y * o.y;
		arr[idx].y = a.x * o.y + a.y * o.x;
	}
}

// divide two complex arrays
__global__ void divide_inplace(PRISM_CUDA_COMPLEX_FLOAT* arr,
                               const PRISM_FLOAT_PRECISION val,
                               const size_t N){
	int idx = threadIdx.x + blockDim.x*blockIdx.x;
	if (idx < N) {
		arr[idx].x /= val;
		arr[idx].y /= val;
	}
}
// set all array values to val
__global__ void setAll(PRISM_FLOAT_PRECISION *data, PRISM_FLOAT_PRECISION val, size_t N) {
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	if (idx<N) {
		data[idx] = val;
	}
}

// creates initial probe using existing GPU memory rather than streaming each probe
__global__ void initializePsi(PRISM_CUDA_COMPLEX_FLOAT *psi_d,
                              const PRISM_CUDA_COMPLEX_FLOAT* PsiProbeInit_d,
                              const PRISM_FLOAT_PRECISION* qya_d,
                              const PRISM_FLOAT_PRECISION* qxa_d,
                              const size_t N,
                              const PRISM_FLOAT_PRECISION yp,
                              const PRISM_FLOAT_PRECISION xp){
	int idx = threadIdx.x + blockDim.x*blockIdx.x;
	if (idx < N) {
		PRISM_CUDA_COMPLEX_FLOAT arg;
		arg = (PRISM_CUDA_COMPLEX_FLOAT)make_cuFloatComplex(qxa_d[idx]*xp + qya_d[idx]*yp, 0);
		psi_d[idx] = cuCmulf(PsiProbeInit_d[idx], exp_cx(cuCmulf(minus_2pii,arg)));
	}
}
// compute modulus squared of other and store in arr
__global__ void abs_squared(PRISM_FLOAT_PRECISION* arr,
                            const PRISM_CUDA_COMPLEX_FLOAT* other,
                            const size_t N){
	int idx = threadIdx.x + blockDim.x*blockIdx.x;
	if (idx < N) {
		PRISM_FLOAT_PRECISION re = other[idx].x;
		PRISM_FLOAT_PRECISION im = other[idx].y;
		arr[idx] = re*re + im*im;
	}
}
__global__ void array_subset(const PRISM_CUDA_COMPLEX_FLOAT* psi_d,
                             PRISM_CUDA_COMPLEX_FLOAT* psi_small_d,
                             const size_t* qyInd_d,
                             const size_t* qxInd_d,
                             const size_t dimi,
                             const size_t dimj_small,
                             const size_t dimi_small,
                             const size_t N){
	int idx = threadIdx.x + blockDim.x*blockIdx.x;
	if (idx < N) {
		int y = idx / dimj_small;
		int x = idx % dimi_small;
		int idxBig = qyInd_d[y]*dimi +  qxInd_d[x];
		psi_small_d[idx] = psi_d[idxBig];

//			arr[idx].x /= val;
//			arr[idx].y /= val;
	}
}