#include "Multislice.cuh"
#include "Multislice.h"
#include "cuComplex.h"
#include "cufft.h"




#include <iostream>
//#include "cuda.h"

#include "fftw3.h"
#include "utility.h"
#include "../../../../../../usr/local/cuda/include/driver_types.h"

#ifdef PRISM_ENABLE_DOUBLE_PRECISION
typedef cuDoubleComplex PRISM_CUDA_COMPLEX_FLOAT;
#else
typedef cuFloatComplex PRISM_CUDA_COMPLEX_FLOAT;
#endif //PRISM_ENABLE_DOUBLE_PRECISION
#define NX 64
#define NY 64
#define NZ 128
#define PI 3.14159265359
#define BLOCK_SIZE1D 1024

// helpful function for checking CUDA errors.
// Source: http://stackoverflow.com/questions/14038589/what-is-the-canonical-way-to-check-for-errors-using-the-cuda-runtime-api
#define cudaErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true){
	if (code != cudaSuccess)
	{
		fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}

#define cufftErrchk(ans) { gpuAssert_cufft((ans), __FILE__, __LINE__); }
inline void gpuAssert_cufft(int code, const char *file, int line, bool abort=true){
	if (code != CUFFT_SUCCESS)
	{
		fprintf(stderr,"GPUassert: %s %d\n", file, line);
		if (abort) exit(code);
	}
}

namespace PRISM{

	__device__ __constant__ PRISM_FLOAT_PRECISION pi       = PI;
	__device__ __constant__ PRISM_CUDA_COMPLEX_FLOAT i     = {0, 1};
	__device__ __constant__ PRISM_CUDA_COMPLEX_FLOAT pi_cx = {PI, 0};
	__device__ __constant__ PRISM_CUDA_COMPLEX_FLOAT minus_2pii = {0, -2*PI};
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
		arg = make_cuFloatComplex(qxa_d[idx]*xp + qya_d[idx]*yp, 0);
		psi_d[idx] = cuCmulf(PsiProbeInit_d[idx], exp_cx(cuCmulf(minus_2pii,arg)));
	}
}
	__global__ void kernel_multiply_inplace(PRISM_CUDA_COMPLEX_FLOAT* arr,
	                                        const PRISM_CUDA_COMPLEX_FLOAT* other,
	                                        const size_t N){
		int idx = threadIdx.x + blockDim.x*blockIdx.x;
		if (idx < N) {
			PRISM_CUDA_COMPLEX_FLOAT a = arr[idx];
			PRISM_CUDA_COMPLEX_FLOAT o = other[idx];

			arr[idx].x = a.x*o.x - a.y*o.y;
			arr[idx].y = a.x*o.y + a.y*o.x;
		}
	}
//	__global__ void kernel_multiply_inplace(PRISM_CUDA_COMPLEX_FLOAT* arr,
//	                                        const PRISM_CUDA_COMPLEX_FLOAT* other,
//	                                        const size_t N){
//		int idx = threadIdx.x + blockDim.x*blockIdx.x;
//		if (idx < N) {
//			arr[idx] = cuCmulf(arr[idx], other[idx]);
//			PRISM_FLOAT_PRECISION ax = arr[idx].x;
//			PRISM_FLOAT_PRECISION ay = arr[idx].y;
//			PRISM_FLOAT_PRECISION ox = other[idx].x;
//			PRISM_FLOAT_PRECISION oy = other[idx].y;
//
//			arr[idx].x = ax*ox - ay*oy;
//			arr[idx].y = ax*oy + ay*ox;
//		}
//	}


	__global__ void kernel_divide_inplace(PRISM_CUDA_COMPLEX_FLOAT* arr,
	                                        const PRISM_FLOAT_PRECISION val,
	                                        const size_t N){
		int idx = threadIdx.x + blockDim.x*blockIdx.x;
		if (idx < N) {
			arr[idx].x /= val;
			arr[idx].y /= val;
		}
	}

//	__global__ void kernel_divide_inplace(PRISM_CUDA_COMPLEX_FLOAT* arr,
//	                                      const PRISM_CUDA_COMPLEX_FLOAT* other,
//	                                      const size_t N){
//		int idx = threadIdx.x + blockDim.x*blockIdx.x;
//		if (idx < N) {
//			arr[idx].x *= other[idx].x;
//			arr[idx].y *= other[idx].y;
//		}
//	}

	__global__ void abs_squared(PRISM_FLOAT_PRECISION* arr,
	                            const PRISM_CUDA_COMPLEX_FLOAT* other,
	                            const size_t N){
		int idx = threadIdx.x + blockDim.x*blockIdx.x;
		if (idx < N) {
//			arr[idx] = cuCabsf(other[idx]);
			PRISM_FLOAT_PRECISION re = other[idx].x;
			PRISM_FLOAT_PRECISION im = other[idx].y;
			arr[idx] = re*re + im*im;
		}
	}

	__global__ void setAll(PRISM_FLOAT_PRECISION *data, PRISM_FLOAT_PRECISION val, size_t N) {
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	if (idx<N) {
		data[idx] = val;
	}
	}
	__global__ void integrateDetector(const PRISM_FLOAT_PRECISION* psiabssquared_d,
	                                  const PRISM_FLOAT_PRECISION* alphaInd_d,
	                                  PRISM_FLOAT_PRECISION* integratedOutput,
	                                  const size_t N,
	                                  const size_t num_integration_bins) {
		int idx = threadIdx.x + blockDim.x * blockIdx.x;
		if (idx < N) {
			size_t alpha = (size_t)alphaInd_d[idx];
			if (alpha <= num_integration_bins)
				atomicAdd(&integratedOutput[alpha-1], psiabssquared_d[idx]);
//				atomicAdd(&integratedOutput[alpha], 1);
		}
	}
  void formatOutput_gpu_integrate(Parameters<PRISM_FLOAT_PRECISION>& pars,
                                  PRISM_FLOAT_PRECISION * psiabssquared_d,
                                  const PRISM_FLOAT_PRECISION * alphaInd_d,
//                                  PRISM_FLOAT_PRECISION* integratedOutput_d,
                                  const size_t& ay,
                                  const size_t& ax,
                                  const size_t& dimj,
                                  const size_t& dimi,
                                  cudaStream_t& stream){
	  PRISM_FLOAT_PRECISION * integratedOutput_d;
	  size_t num_integration_bins = pars.detectorAngles.size();
	  cudaErrchk(cudaMalloc((void **) &integratedOutput_d,
	                        num_integration_bins * sizeof(PRISM_FLOAT_PRECISION)));
//	  cudaErrchk(cudaMemset(integratedOutput, 1, pars.detectorAngles.size() * sizeof(PRISM_FLOAT_PRECISION)));
	  setAll<<< (num_integration_bins - 1)/BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream>>>(integratedOutput_d, 0, num_integration_bins);
	  integrateDetector<<< (dimj*dimi - 1)/BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream>>>(psiabssquared_d, alphaInd_d, integratedOutput_d, dimj*dimi, num_integration_bins);

	  cudaErrchk(cudaMemcpy(&pars.stack.at(ay,ax,0,0),integratedOutput_d,
	                        num_integration_bins * sizeof(PRISM_FLOAT_PRECISION),
	                        cudaMemcpyDeviceToHost));

	  cudaErrchk(cudaFree(integratedOutput_d));
//	  int a;

}

__host__ void getMultisliceProbe_gpu(Parameters<PRISM_FLOAT_PRECISION>& pars,
									 PRISM_CUDA_COMPLEX_FLOAT* trans_d,
									 PRISM_CUDA_COMPLEX_FLOAT* PsiProbeInit_d,
//									 PRISM_FLOAT_PRECISION* integratedOutput_d,
									 const PRISM_FLOAT_PRECISION* qya_d,
									 const PRISM_FLOAT_PRECISION* qxa_d,
									 const PRISM_CUDA_COMPLEX_FLOAT* prop_d,
									 const size_t& ay,
									 const size_t& ax,
									 const size_t dimj,
									 const size_t dimi,
									 const PRISM_FLOAT_PRECISION* alphaInd_d,
									 cudaStream_t& stream){
//		cout << "stream = " << stream << endl;
//		cout << "ax = " << ax <<endl;
//		cout << "ay = " << ay <<endl;
		// create cuFFT plan
		cufftHandle plan;
		cufftErrchk(cufftPlan2d(&plan, dimi, dimj, CUFFT_C2C));

//		cudaStream_t stream2;
//		cudaErrchk(cudaStreamCreate(&stream2));
		// set the stream
//		cufftErrchk(cufftSetStream(plan, stream2));
		cufftErrchk(cufftSetStream(plan, stream));

		// initialize psi
		PRISM_FLOAT_PRECISION yp = pars.yp[ay];
		PRISM_FLOAT_PRECISION xp = pars.xp[ax];
		PRISM_CUDA_COMPLEX_FLOAT *psi_d;
		PRISM_FLOAT_PRECISION *psi_abssquared_d;
		const size_t N = dimj*dimi;
		cudaErrchk(cudaMalloc((void**)&psi_d, dimj*dimi*sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
		cudaErrchk(cudaMalloc((void**)&psi_abssquared_d, dimj*dimi*sizeof(PRISM_FLOAT_PRECISION)));
//		cout << " xp = " << xp << endl;
//		cout << " yp = " << yp << endl;


		//initializePsi<<<(N-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_d, PsiProbeInit_d, qya_d, qxa_d, dimj*dimi, yp, xp);
		initializePsi<<<(N-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_d, PsiProbeInit_d, qya_d, qxa_d, dimj*dimi, yp, xp);

//		initializePsi<<<(N-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D>>>(psi_d, PsiProbeInit_d, qya_d, qxa_d, dimj*dimi, yp, xp);
//
//
//		for (auto planeNum = 0; planeNum < 1; ++planeNum) {
//			cufftExecC2C(plan, psi_d, psi_d, CUFFT_INVERSE);
//			kernel_multiply_inplace<<<(N-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream2>>>(psi_d, &trans_d[planeNum*N], N);
//			cufftExecC2C(plan, psi_d, psi_d, CUFFT_FORWARD);
//			kernel_multiply_inplace<<<(N-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream2>>>(psi_d, prop_d, N);
//			kernel_divide_inplace<<<(N-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream2>>>(psi_d, N, N);
//		}


		for (auto planeNum = 0; planeNum < pars.numPlanes; ++planeNum) {
			cufftErrchk(cufftExecC2C(plan, &psi_d[0], &psi_d[0], CUFFT_INVERSE));
			kernel_multiply_inplace<<<(N-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_d, &trans_d[planeNum*N], N);
			cufftErrchk(cufftExecC2C(plan, &psi_d[0], &psi_d[0], CUFFT_FORWARD));
			kernel_multiply_inplace<<<(N-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_d, prop_d, N);
			kernel_divide_inplace<<<(N-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_d, N, N);
		}

		abs_squared<<<(N-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_abssquared_d, psi_d, N);
		//cudaMemcpyAsync(output, psi_abssquared_d,N*sizeof(PRISM_FLOAT_PRECISION), cudaMemcpyDeviceToHost,stream2);
//	cudaErrchk(cudaMemcpy(output, psi_abssquared_d,N*sizeof(PRISM_FLOAT_PRECISION), cudaMemcpyDeviceToHost));


	formatOutput_gpu_integrate(pars, psi_abssquared_d, alphaInd_d, ay, ax, dimj, dimi,stream);

//		if(ay==0 && ax==0) {
//			for (auto j = 0; j < 25; ++j)cout << "output[j] = " << output[j] << endl;
////		Array2D<PRISM_FLOAT_PRECISION> = zeros_ND<2, PRISM_FLOAT_PRECISION>({{dimj, dimi}});
//		}
//
//	std::complex<PRISM_FLOAT_PRECISION> answer;
//
////	cudaMemcpy(&answer, psi_d+(ax),1*sizeof(PRISM_CUDA_COMPLEX_FLOAT),cudaMemcpyDeviceToHost);
//////	cout << " answer = " << answer << endl;
//		if (ax==0 && ay==0) {
//			for (auto j = 0; j < 10; ++j) {
//				cudaErrchk(cudaMemcpy(&answer, psi_d + j, 1 * sizeof(PRISM_CUDA_COMPLEX_FLOAT), cudaMemcpyDeviceToHost));
//				cout << " psi_d = " << answer << endl; //<< "xp = " << xp << "yp = " << yp << endl;
//			}
//		}
//
//		PRISM_FLOAT_PRECISION answer2;
//
////	cudaMemcpy(&answer, psi_d+(ax),1*sizeof(PRISM_CUDA_COMPLEX_FLOAT),cudaMemcpyDeviceToHost);
//////	cout << " answer = " << answer << endl;
//		if (ax==0 && ay==0) {
//			for (auto j = 0; j < 10; ++j) {
//				cudaErrchk(cudaMemcpy(&answer2, psi_abssquared_d + j, 1 * sizeof(PRISM_FLOAT_PRECISION), cudaMemcpyDeviceToHost));
//				cout << " psi_abssquared_d = " << answer2 << endl << "xp = " << xp << "yp = " << yp << endl;
//			}
//		}
//
//		std::complex<PRISM_FLOAT_PRECISION> answer_propd;
//
////	cudaMemcpy(&answer, psi_d+(ax),1*sizeof(PRISM_CUDA_COMPLEX_FLOAT),cudaMemcpyDeviceToHost);
//////	cout << " answer = " << answer << endl;
//		if (ax==0 && ay==0) {
//			for (auto j = 0; j < 10; ++j) {
//				cudaErrchk(cudaMemcpy(&answer_propd, prop_d + j, 1 * sizeof(std::complex<PRISM_FLOAT_PRECISION>), cudaMemcpyDeviceToHost));
//				cout << " prop_d = " << answer_propd << endl;
//			}
//		}
//
//		std::complex<PRISM_FLOAT_PRECISION> answer_trans_d;
//
////	cudaMemcpy(&answer, psi_d+(ax),1*sizeof(PRISM_CUDA_COMPLEX_FLOAT),cudaMemcpyDeviceToHost);
//////	cout << " answer = " << answer << endl;
//		if (ax==0 && ay==0) {
//			for (auto j = 0; j < 10; ++j) {
//				cudaErrchk(cudaMemcpy(&answer_trans_d, trans_d + j, 1 * sizeof(std::complex<PRISM_FLOAT_PRECISION>), cudaMemcpyDeviceToHost));
//				cout << " answer_trans_d = " << answer_trans_d << endl;
//			}
//		}
//
//

/*
psi(:) = PsiProbeInit ...
		.* exp(-2i*pi ...
		*(emdSTEM.qxa*emdSTEM.MULTIxp(a0) ...
		+ emdSTEM.qya*emdSTEM.MULTIyp(a1)));


copy PsiProbeInit,qxa, qya to GPU cuComplex types. Allocate a pinned output that is same type as stack
run kernel to accomplish above line. If necessary, write a cuComplex mult(const cuComplex& a,const cuComplex& b) function and a cuComplex my_exp(const cuComplex& a,const cuComplex& b) that uses exp(a + bi) = exp(a) *cos(a) + exp(a) * sin(a)*i

create two cufft plans, write a complex inplace multiplication kernel, and apply the following loop
for a2 = 1:emdSTEM.numPlanes
		psi = fft2(ifft2(psi).*trans(:,:,a2)).*emdSTEM.prop;
	end

take final FFT and abs square

async stream the result to the pinned stack

after done copy the pinned stack to original

*/
	cufftErrchk(cufftDestroy(plan));
	cudaErrchk(cudaFree(psi_d));
	cudaErrchk(cudaFree(psi_abssquared_d));
//	cudaErrchk(cudaStreamDestroy(stream2));
}
    __host__ void buildMultisliceOutput_gpu(Parameters <PRISM_FLOAT_PRECISION> &pars,
                                            Array3D <std::complex<PRISM_FLOAT_PRECISION>> &trans,
                                            Array2D <std::complex<PRISM_FLOAT_PRECISION>> &PsiProbeInit,
                                            Array2D <PRISM_FLOAT_PRECISION> &alphaInd) {
		cout << "debug pars.prop" << endl;
		for (auto i =0; i < 10; ++i)cout << pars.prop[i] << endl;
		// populate the Multislice output stack dividing the work between GPUs and CPU cores.
		// this version assumes the full trans array fits into DRAM on each GPU

		using namespace std;
		cout << "DEBUG PsiProbeInit.at(0,1) = " << PsiProbeInit.at(0,1) << endl;
		cout << "DEBUG PsiProbeInit[1] = " << PsiProbeInit[1] << endl;
		cout << "DEBUG PsiProbeInit[2] = " << PsiProbeInit[2] << endl;

		// create CUDA streams
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;
		cudaStream_t streams[total_num_streams];
		cout <<"total_num_streams = " << total_num_streams<< endl;
		for (auto j = 0; j < total_num_streams; ++j){
			cudaSetDevice(j % pars.meta.NUM_GPUS);
			cudaErrchk(cudaStreamCreate(&streams[j]));
		}
		cout << "first stream = " << streams[0] << endl;
        cout << "Test GPU function from CUDA host" << endl;
		const PRISM_FLOAT_PRECISION cpu_stop = std::floor(pars.meta.cpu_gpu_ratio*pars.yp.size());
		vector<thread> workers_gpu;
		vector<thread> workers_cpu;
		workers_gpu.reserve(total_num_streams); // prevents multiple reallocations
		workers_cpu.reserve(pars.meta.NUM_THREADS); // prevents multiple reallocations

		// pointers to GPU memory (one copy per GPU)
		PRISM_CUDA_COMPLEX_FLOAT *PsiProbeInit_d[pars.meta.NUM_GPUS];
		PRISM_CUDA_COMPLEX_FLOAT *trans_d[pars.meta.NUM_GPUS];
		PRISM_CUDA_COMPLEX_FLOAT *prop_d[pars.meta.NUM_GPUS];
		PRISM_FLOAT_PRECISION    *qxa_d[pars.meta.NUM_GPUS];
		PRISM_FLOAT_PRECISION    *qya_d[pars.meta.NUM_GPUS];
	    PRISM_FLOAT_PRECISION    *alphaInd_d[pars.meta.NUM_GPUS];
//		PRISM_FLOAT_PRECISION    *integratedOutput_d[pars.meta.NUM_GPUS];

		// allocate memory on each GPU
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaMalloc((void **) &PsiProbeInit_d[g],     PsiProbeInit.size()        * sizeof(PsiProbeInit[0])));
			cudaErrchk(cudaMalloc((void **) &trans_d[g],            trans.size()               * sizeof(trans[0])));
			cudaErrchk(cudaMalloc((void **) &qxa_d[g],              pars.qxa.size()            * sizeof(pars.qxa[0])));
			cudaErrchk(cudaMalloc((void **) &qya_d[g],              pars.qya.size()            * sizeof(pars.qya[0])));
			cudaErrchk(cudaMalloc((void **) &prop_d[g],             pars.prop.size()           * sizeof(pars.prop[0])));
			cudaErrchk(cudaMalloc((void **) &alphaInd_d[g],         alphaInd.size()            * sizeof(alphaInd[0])));
//			cudaErrchk(cudaMalloc((void **) &integratedOutput_d[g], pars.detectorAngles.size() * sizeof(alphaInd[0])));
		}

		// copy memory to each GPU (this can be made asynchronous if necessary by copying to pinned memory first)
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaMemcpy(PsiProbeInit_d[g], &PsiProbeInit[0],
			                      PsiProbeInit.size() * sizeof(PsiProbeInit[0]), cudaMemcpyHostToDevice));
			cudaErrchk(cudaMemcpy(trans_d[g], &trans[0],
			                      trans.size() * sizeof(trans[0]), cudaMemcpyHostToDevice));
			cudaErrchk(cudaMemcpy(qxa_d[g], &pars.qxa[0],
			                      pars.qxa.size() * sizeof(pars.qxa[0]), cudaMemcpyHostToDevice));
			cudaErrchk(cudaMemcpy(qya_d[g], &pars.qya[0],
			                      pars.qya.size() * sizeof(pars.qya[0]), cudaMemcpyHostToDevice));
			cudaErrchk(cudaMemcpy(prop_d[g], &pars.prop[0],
			                      pars.prop.size() * sizeof(pars.prop[0]), cudaMemcpyHostToDevice));
			cudaErrchk(cudaMemcpy(alphaInd_d[g], &alphaInd[0],
			                      alphaInd.size() * sizeof(alphaInd[0]), cudaMemcpyHostToDevice));
		}
//
//		std::complex<float> b;
//		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
//			cudaSetDevice(g);
//			cudaMemcpy(&b, PsiProbeInit_d[g]+1, 1 * sizeof(PsiProbeInit[0]), cudaMemcpyDeviceToHost);
//			cout << "this b = " << b << endl;
//		}

		size_t psi_size = PsiProbeInit.size();

		auto WORK_CHUNK_SIZE_GPU = std::floor( (((pars.yp.size() - cpu_stop - 1)) / total_num_streams) + 1); //TODO: divide work more generally than just splitting up by yp. If input isn't square this might not do a good job
		cout << "WORK_CHUNK_SIZE_GPU = " << WORK_CHUNK_SIZE_GPU << endl;
		auto start = cpu_stop;// gpu work starts where the cpu work will stop
		auto stop = start + WORK_CHUNK_SIZE_GPU;
		int stream_count = 0;
		while (start < pars.yp.size()) {
			int gpu_num = stream_count % pars.meta.NUM_GPUS; // determine which gpu handles this job
			cudaStream_t& current_stream = streams[stream_count];
			cout << " STREAM = " << streams[stream_count] << endl;
			cout << "Launching thread to compute all x-probe positions for y-probes "
				 << start << "-" << std::min((size_t)stop,pars.yp.size()) << " on stream #" << stream_count << " of GPU #" << gpu_num << '\n';

			// get pointers to the pre-copied arrays, making sure to get those on the current GPU
			PRISM_CUDA_COMPLEX_FLOAT *current_PsiProbeInit_d = PsiProbeInit_d[gpu_num];
			PRISM_CUDA_COMPLEX_FLOAT *current_trans_d = trans_d[gpu_num];
			PRISM_CUDA_COMPLEX_FLOAT *current_prop_d  = prop_d[gpu_num];
			PRISM_FLOAT_PRECISION *current_qxa_d      = qxa_d[gpu_num];
			PRISM_FLOAT_PRECISION *current_qya_d      = qya_d[gpu_num];
			PRISM_FLOAT_PRECISION *current_alphaInd_d = alphaInd_d[gpu_num];
//			PRISM_FLOAT_PRECISION *current_integratedOutput_d = integratedOutput_d[gpu_num];

			// launch a new thread
			// emplace_back is better whenever constructing a new object
			workers_gpu.emplace_back(thread([&pars, current_trans_d, current_PsiProbeInit_d, &alphaInd, current_alphaInd_d,
					                                start, stop, gpu_num, current_qya_d, current_qxa_d,
					                                current_prop_d, &current_stream, &psi_size, &PsiProbeInit]() {

				// page-locked (pinned) memory for async streaming of the result back
//				PRISM_FLOAT_PRECISION *pinned_output;

				// figure out how much pinned memory to allocate in this job
				size_t num_probe_positions = pars.xp.size() * (min((size_t) stop, pars.yp.size()) - start);
//				size_t pinned_output_size = psi_size * num_probe_positions;

				// allocate pinned memory on host
//				cout << "pinned_output_size * sizeof(PRISM_FLOAT_PRECISION)) = " << pinned_output_size * sizeof(PRISM_FLOAT_PRECISION) << endl;
//				cout << "num_probe_positions = " << num_probe_positions<< endl;
//				cout << "psi_size = " << psi_size << endl;
//				cout << "PsiProbeInit.get_dimj() = " << PsiProbeInit.get_dimj() <<endl;
//				cout << "PsiProbeInit.get_dimi() = " << PsiProbeInit.get_dimi() <<endl;
//				cudaErrchk(cudaMallocHost((void **) &pinned_output, pinned_output_size * sizeof(PRISM_FLOAT_PRECISION)));

				// set the GPU context
				cudaErrchk(cudaSetDevice(gpu_num)); // set current gpu
//				PRISM_FLOAT_PRECISION *pinned_output_begin = pinned_output; // pointer to the beginning of corresponding output layer in the 3D array
				for (auto ay = start; ay < std::min((size_t) stop, pars.yp.size()); ++ay) {
					for (auto ax = 0; ax < pars.xp.size(); ++ax) {
//				for (auto ay = start; ay < 25; ++ay) {
//					for (auto ax = 0; ax < 25; ++ax) {
						getMultisliceProbe_gpu(pars, current_trans_d, current_PsiProbeInit_d, current_qya_d,
						                       current_qxa_d,
						                       current_prop_d, ay, ax, PsiProbeInit.get_dimj(), PsiProbeInit.get_dimi(),
						                       current_alphaInd_d, current_stream);
//						pinned_output_begin += psi_size; // advance the start point of the output
					}
				}

//				cudaErrchk(cudaDeviceSynchronize());
//				{
//					for (auto j = 0; j < 25; ++j)cout << "pinned_output[j] = " << pinned_output[j] << endl;
//				}
//				PRISM_FLOAT_PRECISION c = 0;
//				for (auto i = 0; i < pinned_output_size; ++i) {
//					c += pinned_output[i];
//				}
//				cout << "c = " << c << endl;
//				PRISM_FLOAT_PRECISION *counts = pinned_output;

//				Array2D<PRISM_FLOAT_PRECISION> db = zeros_ND<2, PRISM_FLOAT_PRECISION>({{PsiProbeInit.get_dimj(), PsiProbeInit.get_dimi()}});
//				auto db_ptr = db.begin();
//				for (auto ay = start; ay < 25; ++ay) {
//					for (auto ax = 0; ax < 25; ++ax) {
//				for (auto ay = start; ay < std::min((size_t) stop, pars.yp.size()); ++ay) {
//					for (auto ax = 0; ax < pars.xp.size(); ++ax) {
//						auto idx = alphaInd.begin();
//						while (idx != alphaInd.end()) {
//							if (ay==0 && ax==0)*db_ptr = *counts;
//							if (*idx <= pars.Ndet) {
//								cout << "count = " << *counts << endl;
//								pars.stack.at(ay, ax, (*idx) - 1, 0) += (*counts);

//								cout << "ax = " << ax << endl;
//								cout << "ay = " << ay << endl;
//								cout << "(*idx) - 1 = " << (*idx) - 1 << endl;
//								cout << "pars.stack.at(ay, ax, (*idx) - 1, 0) = "
//								     << pars.stack.at(ay, ax, (*idx) - 1, 0)
//								     << endl;
//							}
//							db_ptr++;
//							++idx;
//							++counts;
//						}
//					}
//				}
//				db.toMRC_f("db_intOutput.mrc");
//			cout << " DEBUG 1 " << endl;
//				if (start == 0 ) {
//					Array2D <PRISM_FLOAT_PRECISION> prism_image;
//					prism_image = zeros_ND<2, PRISM_FLOAT_PRECISION>({{pars.stack.get_diml(), pars.stack.get_dimk()}});
//					for (auto y = 0; y < pars.stack.get_diml(); ++y) {
//						for (auto x = 0; x < pars.stack.get_dimk(); ++x) {
//							for (auto b = 13; b < 18; ++b) {
//								prism_image.at(y, x) += pars.stack.at(y, x, b, 0);
////								cout << "prism_image.at(y, x) = " << prism_image.at(y, x) << endl;
//							}
//						}
//					}
////					prism_image.toMRC_f("TEST.mrc");
//					cout <<" debug written" <<endl;
//				}
//				for (auto y = 0; y < PsiProbeInit.get_dimj(),; ++y) {
//					for (auto x = 0; x < PsiProbeInit.get_dimi(),; ++x) {
//						auto idx = alphaInd.at(y,x);
//						if ( alphaInd.at(y,x) <= pars.Ndet){
//							pars.stack.at(ay,ax,(*idx)-1, 0) += *counts * pars.scale;
//						}
//					}
//				}
//				auto idx = alphaInd.begin();
//				for (auto counts = intOutput.begin(); counts != intOutput.end(); ++counts){
//					if (*idx <= pars.Ndet){
//						pars.stack.at(ay,ax,(*idx)-1, 0) += *counts * pars.scale;
//					}
//					++idx;
//				};

//				auto stack_ptr = &pars.stack[start*];
//				cudaErrchk(cudaFreeHost(pinned_output));
			}));

			stream_count++;
			start += WORK_CHUNK_SIZE_GPU;
			if (start >= pars.yp.size())break;
			stop += WORK_CHUNK_SIZE_GPU;
		}


		// now launch CPU work
		auto WORK_CHUNK_SIZE_CPU = std::floor(((cpu_stop - 1) / pars.meta.NUM_THREADS) + 1); //TODO: divide work more generally than just splitting up by yp. If input isn't square this might not do a good job
		cout << "WORK_CHUNK_SIZE_CPU = " << WORK_CHUNK_SIZE_CPU << endl;
                start = 0;// cpu work starts at beginning
                stop = start + WORK_CHUNK_SIZE_CPU;
                while (start < cpu_stop) {
                        cout << "Launching thread to compute all x-probe positions for y-probes "
                                 << start << "-" << std::min(stop,cpu_stop) << " on CPU\n";
                        // emplace_back is better whenever constructing a new object
                        workers_cpu.emplace_back(thread([&pars, &trans,
                                                                                                &alphaInd, &PsiProbeInit,
                                                                                                start, cpu_stop,stop]() {
                                for (auto ay = start; ay < std::min(stop, cpu_stop); ++ay) {
                                        for (auto ax = 0; ax < pars.xp.size(); ++ax) {
                                                getMultisliceProbe_cpu(pars, trans, PsiProbeInit, ay, ax, alphaInd);
                                        }
                                }
                        }));
                        start += WORK_CHUNK_SIZE_CPU;
                        if (start >= cpu_stop)break;
                        stop += WORK_CHUNK_SIZE_CPU;
                }

		// synchronize threads
		cout << "waiting on threads" << endl;
		for (auto& t:workers_gpu)t.join();
		for (auto& t:workers_cpu)t.join();
		cout << "threads done, cleaning up" << endl;
		// synchronize GPUs and cleanup data
		for (auto j = 0; j < pars.meta.NUM_GPUS; ++j){
			cudaErrchk(cudaSetDevice(j));
			cudaErrchk(cudaDeviceSynchronize());
			cudaErrchk(cudaFree(PsiProbeInit_d[j]));
			cudaErrchk(cudaFree(trans_d[j]));
			cudaErrchk(cudaFree(qxa_d[j]));
			cudaErrchk(cudaFree(qya_d[j]));
			cudaErrchk(cudaFree(prop_d[j]));
			cudaErrchk(cudaFree(alphaInd_d[j]));
//			cudaErrchk(cudaFree(integratedOutput_d[j]));
		}
		// destroy CUDA streams
		for (auto j = 0; j < total_num_streams; ++j)cudaErrchk(cudaStreamDestroy(streams[j]));
		for (auto j = 0; j < pars.meta.NUM_GPUS; ++j) {
			cudaErrchk(cudaSetDevice(j));
			cudaErrchk(cudaDeviceReset());
		}
	}
}
