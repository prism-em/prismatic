#include "utility.cuh"
#include <iostream>

#define PI 3.14159265359
// define some constants
__device__ __constant__ float pi_f       = PI;
__device__ __constant__ cuFloatComplex i_f     = {0, 1};
__device__ __constant__ cuFloatComplex pi_cx_f = {PI, 0};
__device__ __constant__ cuFloatComplex minus_2pii_f = {0, -2*PI};
__device__ __constant__ double pi       = PI;
__device__ __constant__ cuDoubleComplex i     = {0, 1};
__device__ __constant__ cuDoubleComplex pi_cx = {PI, 0};
__device__ __constant__ cuDoubleComplex minus_2pii = {0, -2*PI};

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
__global__ void multiply_inplace(cuDoubleComplex* arr,
                                 const cuDoubleComplex* other,
                                 const size_t N){
	int idx = threadIdx.x + blockDim.x*blockIdx.x;
	if (idx < N) {
		cuDoubleComplex a = arr[idx];
		cuDoubleComplex o = other[idx];
		arr[idx].x = a.x * o.x - a.y * o.y;
		arr[idx].y = a.x * o.y + a.y * o.x;
	}
}

// multiply two complex arrays
__global__ void multiply_inplace(cuFloatComplex* arr,
                                 const cuFloatComplex* other,
                                 const size_t N){
	int idx = threadIdx.x + blockDim.x*blockIdx.x;
	if (idx < N) {
		cuFloatComplex a = arr[idx];
		cuFloatComplex o = other[idx];
		arr[idx].x = a.x * o.x - a.y * o.y;
		arr[idx].y = a.x * o.y + a.y * o.x;
	}
}

// multiply two complex arrays
__global__ void multiply_cx(cuDoubleComplex* arr,
                             const cuDoubleComplex* other,
                             const size_t N){
	int idx = threadIdx.x + blockDim.x*blockIdx.x;
	if (idx < N) {
//		cuDoubleComplex a = arr[idx];
//		cuDoubleComplex o = other[idx];
		arr[idx] = cuCmul(arr[idx], other[idx]);
	}
}

// multiply two complex arrays
__global__ void multiply_cx(cuFloatComplex* arr,
                            const cuFloatComplex* other,
                            const size_t N){
	int idx = threadIdx.x + blockDim.x*blockIdx.x;
	if (idx < N) {
//		cuFloatComplex a = arr[idx];
//		cuFloatComplex o = other[idx];
		arr[idx] = cuCmulf(arr[idx], other[idx]);
	}
}

//// divide two complex arrays
//__global__ void divide_inplace(PRISM_CUDA_COMPLEX_FLOAT* arr,
//                               const PRISM_FLOAT_PRECISION val,
//                               const size_t N){
//	int idx = threadIdx.x + blockDim.x*blockIdx.x;
//	if (idx < N) {
//		arr[idx].x /= val;
//		arr[idx].y /= val;
//	}
//}

__global__ void divide_inplace(cuDoubleComplex* arr,
                               const cuDoubleComplex val,
                               const size_t N){
	int idx = threadIdx.x + blockDim.x*blockIdx.x;
	if (idx < N) {
		arr[idx] = cuCdiv(arr[idx], val);
	}
}

__global__ void divide_inplace(cuFloatComplex* arr,
                               const cuFloatComplex val,
                               const size_t N){
	int idx = threadIdx.x + blockDim.x*blockIdx.x;
	if (idx < N) {
		arr[idx] = cuCdivf(arr[idx], val);
	}
}

// set all array values to val
__global__ void setAll(double *data, double val, size_t N) {
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	if (idx<N) {
		data[idx] = val;
	}
}

// set all array values to val
__global__ void setAll(float *data, float val, size_t N) {
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	if (idx<N) {
		data[idx] = val;
	}
}

// creates initial probe using existing GPU memory rather than streaming each probe
__global__ void initializePsi(cuDoubleComplex *psi_d,
                              const cuDoubleComplex* PsiProbeInit_d,
                              const double* qya_d,
                              const double* qxa_d,
                              const size_t N,
                              const double yp,
                              const double xp){
	int idx = threadIdx.x + blockDim.x*blockIdx.x;
	if (idx < N) {
		cuDoubleComplex arg;
		arg = make_cuDoubleComplex(qxa_d[idx]*xp + qya_d[idx]*yp, 0);
		psi_d[idx] = cuCmul(PsiProbeInit_d[idx], exp_cx(cuCmul(minus_2pii,arg)));
	}
}

// creates initial probe using existing GPU memory rather than streaming each probe
__global__ void initializePsi(cuFloatComplex *psi_d,
                              const cuFloatComplex* PsiProbeInit_d,
                              const float* qya_d,
                              const float* qxa_d,
                              const size_t N,
                              const float yp,
                              const float xp){
	int idx = threadIdx.x + blockDim.x*blockIdx.x;
	if (idx < N) {
		cuFloatComplex arg;
		arg = make_cuFloatComplex(qxa_d[idx]*xp + qya_d[idx]*yp, 0);
		psi_d[idx] = cuCmulf(PsiProbeInit_d[idx], exp_cx(cuCmulf(minus_2pii_f,arg)));
	}
}


// compute modulus squared of other and store in arr
__global__ void abs_squared(double* arr,
                            const cuDoubleComplex* other,
                            const size_t N){
	int idx = threadIdx.x + blockDim.x*blockIdx.x;
	if (idx < N) {
		double re = other[idx].x;
		double im = other[idx].y;
		arr[idx] = re*re + im*im;
	}
}

// compute modulus squared of other and store in arr
__global__ void abs_squared(float* arr,
                            const cuFloatComplex* other,
                            const size_t N){
	int idx = threadIdx.x + blockDim.x*blockIdx.x;
	if (idx < N) {
		float re = other[idx].x;
		float im = other[idx].y;
		arr[idx] = re*re + im*im;
	}
}

__global__ void array_subset(const cuDoubleComplex* psi_d,
                             cuDoubleComplex* psi_small_d,
                             const size_t* qyInd_d,
                             const size_t* qxInd_d,
                             const size_t dimi,
                             const size_t dimj_small,
                             const size_t dimi_small){
	int idx = threadIdx.x + blockDim.x*blockIdx.x;
	if (idx < dimj_small*dimi_small) {
		int y = idx / (int)dimi_small;
		int x = idx % (int)dimi_small;
		int idxBig = qyInd_d[y] * dimi + qxInd_d[x];
		psi_small_d[idx] = psi_d[idxBig];
//		psi_small_d[idx] = make_cuFloatComplex(idx,idxBig);
	}
}
__global__ void array_subset(const cuFloatComplex* psi_d,
                             cuFloatComplex* psi_small_d,
                             const size_t* qyInd_d,
                             const size_t* qxInd_d,
                             const size_t dimi,
                             const size_t dimj_small,
                             const size_t dimi_small){
	int idx = threadIdx.x + blockDim.x*blockIdx.x;
	if (idx < dimj_small*dimi_small) {
		int y = idx / (int)dimi_small;
		int x = idx % (int)dimi_small;
		int idxBig = qyInd_d[y] * dimi + qxInd_d[x];
		psi_small_d[idx] = psi_d[idxBig];
//		psi_small_d[idx] = make_cuFloatComplex(idx,idxBig);
	}
}