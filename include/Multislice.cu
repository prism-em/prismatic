#include "Multislice.cuh"
#include "Multislice.h"
#include "cuComplex.h"
#include "cufft.h"
#include <iostream>
//#include "cuda.h"

#include "fftw3.h"
#include "utility.h"
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
		float e = exp(a.x);
		float s,c;
		sincosf(a.y, &s, &c);
		return make_cuFloatComplex(e*c, e*s);
	}
	__global__ void initializePsi(PRISM_CUDA_COMPLEX_FLOAT *psi_d,
	                              const PRISM_CUDA_COMPLEX_FLOAT* PsiProbeInit_d,
	                              const PRISM_FLOAT_PRECISION* qya_d,
	                              const PRISM_FLOAT_PRECISION* qxa_d,
	                              const size_t dimj,
	                              const size_t dimi,
	                              const PRISM_FLOAT_PRECISION yp,
	                              const PRISM_FLOAT_PRECISION xp){
	int idx = threadIdx.x + blockDim.x*blockIdx.x;
	if (idx < (dimj * dimi)) {
		PRISM_CUDA_COMPLEX_FLOAT arg;
		arg = make_cuFloatComplex(qxa_d[idx]*xp + qya_d[idx]*yp, 0);
		psi_d[idx] = cuCmulf(PsiProbeInit_d[idx], exp_cx(cuCmulf(minus_2pii,arg)));
	}
};
__host__ void getMultisliceProbe_gpu(Parameters<PRISM_FLOAT_PRECISION>& pars,
									 PRISM_CUDA_COMPLEX_FLOAT* trans_d,
									 PRISM_CUDA_COMPLEX_FLOAT* PsiProbeInit_d,
									 const PRISM_FLOAT_PRECISION* qya_d,
									 const PRISM_FLOAT_PRECISION* qxa_d,
									 const size_t& ay,
									 const size_t& ax,
									 const size_t dimj,
									 const size_t dimi,
									 Array2D<PRISM_FLOAT_PRECISION> &alphaInd,
									 const cudaStream_t& stream,
									 std::complex<PRISM_FLOAT_PRECISION>* const output){


	PRISM_FLOAT_PRECISION yp = pars.yp[ay];
	PRISM_FLOAT_PRECISION xp = pars.xp[ax];
//		std::cout << "XP = " << xp << std::endl;
	PRISM_CUDA_COMPLEX_FLOAT *psi_d;
	cudaMalloc((void**)&psi_d, dimj*dimi*sizeof(PRISM_CUDA_COMPLEX_FLOAT));
	initializePsi<<<(dimj*dimi-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D>>>(psi_d, PsiProbeInit_d, qya_d, qxa_d, dimj, dimi, yp, xp);

	std::complex<PRISM_FLOAT_PRECISION> answer;

	cudaMemcpy(&answer, psi_d+(ax),1*sizeof(PRISM_CUDA_COMPLEX_FLOAT),cudaMemcpyDeviceToHost);
//	cout << " answer = " << answer << endl;
		if (ax==0 && ay==0) {
			for (auto j = 0; j < 10; ++j) {
				cudaMemcpy(&answer, psi_d + j, 1 * sizeof(PRISM_CUDA_COMPLEX_FLOAT), cudaMemcpyDeviceToHost);
				cout << " answer = " << answer << endl << "xp = " << xp << "yp = " << yp << endl;
			}
		}
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
	cufftHandle plan;
	cufftComplex *data1, *data2;
	cudaMalloc((void**)&data1, sizeof(cufftComplex)*NX*NY*NZ);
	cudaMalloc((void**)&data2, sizeof(cufftComplex)*NX*NY*NZ);
	/* Create a 3D FFT plan. */
	cufftPlan3d(&plan, NX, NY, NZ, CUFFT_C2C);
	cufftDestroy(plan);
	cudaFree(psi_d);
}
    __host__ void buildMultisliceOutput_gpu(Parameters <PRISM_FLOAT_PRECISION> &pars,
                                            Array3D <std::complex<PRISM_FLOAT_PRECISION>> &trans,
                                            Array2D <std::complex<PRISM_FLOAT_PRECISION>> &PsiProbeInit,
                                            Array2D <PRISM_FLOAT_PRECISION> &alphaInd) {

		// populate the Multislice output stack dividing the work between GPUs and CPU cores.
		// this version assumes the full trans array fits into DRAM on each GPU

		using namespace std;
		cout << "DEBUG PsiProbeInit.at(0,1) = " << PsiProbeInit.at(0,1) << endl;
		cout << "DEBUG PsiProbeInit[1] = " << PsiProbeInit[1] << endl;
		cout << "DEBUG PsiProbeInit[2] = " << PsiProbeInit[2] << endl;
		// create CUDA streams
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;
		cudaStream_t streams[total_num_streams];
		for (auto j = 0; j < total_num_streams; ++j)cudaStreamCreate(&streams[j]);

        cout << "Test GPU function from CUDA host" << endl;
		const PRISM_FLOAT_PRECISION cpu_stop = std::floor(pars.meta.cpu_gpu_ratio*pars.yp.size());
		vector<thread> workers_gpu;
		vector<thread> workers_cpu;
		workers_gpu.reserve(total_num_streams); // prevents multiple reallocations
		workers_cpu.reserve(pars.meta.NUM_THREADS); // prevents multiple reallocations

		// pointers to GPU memory (one copy per GPU)
		PRISM_CUDA_COMPLEX_FLOAT *PsiProbeInit_d[pars.meta.NUM_GPUS];
		PRISM_CUDA_COMPLEX_FLOAT *trans_d[pars.meta.NUM_GPUS];
		PRISM_FLOAT_PRECISION    *qxa_d[pars.meta.NUM_GPUS];
		PRISM_FLOAT_PRECISION    *qya_d[pars.meta.NUM_GPUS];

		// allocate memory on each GPU
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			cudaSetDevice(g);
			cudaMalloc((void **) &PsiProbeInit_d[g], PsiProbeInit.size() * sizeof(PsiProbeInit[0]));
			cudaMalloc((void **) &trans_d[g], trans.size() * sizeof(trans[0]));
			cudaMalloc((void **) &qxa_d[g], pars.qxa.size() * sizeof(pars.qxa[0]));
			cudaMalloc((void **) &qya_d[g], pars.qya.size() * sizeof(pars.qya[0]));
		}

		// copy memory to each GPU (this can be made asynchronous if necessary by copying to pinned memory first)
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			cudaSetDevice(g);
			cudaMemcpy(PsiProbeInit_d[g], &PsiProbeInit[0], PsiProbeInit.size() * sizeof(PsiProbeInit[0]), cudaMemcpyHostToDevice);
			cudaMemcpy(trans_d[g], &trans[0], trans.size() * sizeof(trans[0]), cudaMemcpyHostToDevice);
			cudaMemcpy(qxa_d[g], &pars.qxa[0], pars.qxa.size() * sizeof(pars.qxa[0]), cudaMemcpyHostToDevice);
			cudaMemcpy(qya_d[g], &pars.qya[0], pars.qya.size() * sizeof(pars.qya[0]), cudaMemcpyHostToDevice);
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
			int gpu_num = stream_count++ % pars.meta.NUM_GPUS; // determine which gpu handles this job
			cudaStream_t& current_stream = streams[stream_count];
			cout << "Launching thread to compute all x-probe positions for y-probes "
				 << start << "-" << std::min((size_t)stop,pars.yp.size()) << " on stream #" << stream_count << " of GPU #" << gpu_num << '\n';

			// get pointers to the pre-copied arrays, making sure to get those on the current GPU
			PRISM_CUDA_COMPLEX_FLOAT *current_PsiProbeInit_d = PsiProbeInit_d[gpu_num];
			PRISM_CUDA_COMPLEX_FLOAT *current_trans_d = trans_d[gpu_num];
			PRISM_FLOAT_PRECISION *current_qxa_d = qxa_d[gpu_num];
			PRISM_FLOAT_PRECISION *current_qya_d = qya_d[gpu_num];

			// launch a new thread
			// emplace_back is better whenever constructing a new object
			workers_gpu.emplace_back(thread([&pars, current_trans_d, current_PsiProbeInit_d, &alphaInd,
					                                start, stop, gpu_num, current_qya_d, current_qxa_d,
					                                &current_stream, &psi_size, &PsiProbeInit]() {

				// page-locked (pinned) memory for async streaming of the result back
				std::complex<PRISM_FLOAT_PRECISION>* pinned_stack;

				// figure out how much pinned memory to allocate in this job
				size_t num_probe_positions = pars.xp.size() * (min((size_t)stop, pars.yp.size()) - start);
				size_t pinned_output_size = psi_size * num_probe_positions;

				// allocate pinned memory on host
				cudaMallocHost((void**)&pinned_stack, pinned_output_size);

				// set the GPU context
				cudaSetDevice(gpu_num); // set current gpu
				std::complex<PRISM_FLOAT_PRECISION>* pinned_stack_begin = pinned_stack; // pointer to the beginning of corresponding output layer in the 3D array
				for (auto ay = start; ay < std::min((size_t) stop, pars.yp.size()); ++ay) {
					for (auto ax = 0; ax < pars.xp.size(); ++ax) {
						getMultisliceProbe_gpu(pars, current_trans_d, current_PsiProbeInit_d,current_qya_d, current_qxa_d,
						                       ay, ax, PsiProbeInit.get_dimj(), PsiProbeInit.get_dimi(),
						                       alphaInd, current_stream, pinned_stack_begin);
						pinned_stack_begin += psi_size; // advance the start point of the output
					}
				}
				cudaFreeHost(pinned_stack);
			}));


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
                        workers_gpu.emplace_back(thread([&pars, &trans,
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
		for (auto& t:workers_gpu)t.join();
		for (auto& t:workers_cpu)t.join();
		// synchronize GPUs and cleanup data
		for (auto j = 0; j < pars.meta.NUM_GPUS; ++j){
			cudaSetDevice(j);
			cudaDeviceSynchronize();
			cudaFree(PsiProbeInit_d[j]);
			cudaFree(trans_d[j]);
			cudaFree(qxa_d[j]);
			cudaFree(qya_d[j]);
		}
		// destroy CUDA streams
		for (auto j = 0; j < total_num_streams; ++j)cudaStreamDestroy(streams[j]);

	}
}
