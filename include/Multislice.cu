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
			arr[idx].x = arr[idx].x*other[idx].x - arr[idx].y*other[idx].y;
			arr[idx].y = arr[idx].x*other[idx].y + arr[idx].y*other[idx].x;
		}
	}


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


__host__ void getMultisliceProbe_gpu(Parameters<PRISM_FLOAT_PRECISION>& pars,
									 PRISM_CUDA_COMPLEX_FLOAT* trans_d,
									 PRISM_CUDA_COMPLEX_FLOAT* PsiProbeInit_d,
									 const PRISM_FLOAT_PRECISION* qya_d,
									 const PRISM_FLOAT_PRECISION* qxa_d,
									 const PRISM_CUDA_COMPLEX_FLOAT* prop_d,
									 const size_t& ay,
									 const size_t& ax,
									 const size_t dimj,
									 const size_t dimi,
									 Array2D<PRISM_FLOAT_PRECISION> &alphaInd,
									 const cudaStream_t& stream,
									 PRISM_FLOAT_PRECISION* const output){
		// create cuFFT plan
		cufftHandle plan;
		cufftPlan2d(&plan, dimi, dimj, CUFFT_C2C);

		cudaStream_t stream2;
		cudaStreamCreate(&stream2);
		// set the stream
		cufftSetStream(plan, stream2);

		// initialize psi
		PRISM_FLOAT_PRECISION yp = pars.yp[ay];
		PRISM_FLOAT_PRECISION xp = pars.xp[ax];
		PRISM_CUDA_COMPLEX_FLOAT *psi_d;
		PRISM_FLOAT_PRECISION *psi_abssquared_d;
		const size_t N = dimj*dimi;
		cudaMalloc((void**)&psi_d, dimj*dimi*sizeof(PRISM_CUDA_COMPLEX_FLOAT));
		cudaMalloc((void**)&psi_abssquared_d, dimj*dimi*sizeof(PRISM_FLOAT_PRECISION));
//		cout << " xp = " << xp << endl;
//		cout << " yp = " << yp << endl;


		//initializePsi<<<(N-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_d, PsiProbeInit_d, qya_d, qxa_d, dimj*dimi, yp, xp);
		initializePsi<<<(N-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream2>>>(psi_d, PsiProbeInit_d, qya_d, qxa_d, dimj*dimi, yp, xp);

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
			cufftExecC2C(plan, &psi_d[0], &psi_d[0], CUFFT_INVERSE);
			kernel_multiply_inplace<<<(N-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream2>>>(psi_d, &trans_d[planeNum*N], N);
			cufftExecC2C(plan, &psi_d[0], &psi_d[0], CUFFT_FORWARD);
			kernel_multiply_inplace<<<(N-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream2>>>(psi_d, prop_d, N);
			kernel_divide_inplace<<<(N-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream2>>>(psi_d, N, N);
		}

		abs_squared<<<(N-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream2>>>(psi_abssquared_d, psi_d, N);
		//cudaMemcpyAsync(output, psi_abssquared_d,N*sizeof(PRISM_FLOAT_PRECISION), cudaMemcpyDeviceToHost,stream2);
		cudaMemcpy(output, psi_abssquared_d,N*sizeof(PRISM_FLOAT_PRECISION), cudaMemcpyDeviceToHost);
		if(ay==0 && ax==0) {
			for (auto j = 0; j < 25; ++j)cout << "output[j] = " << output[j] << endl;
//		Array2D<PRISM_FLOAT_PRECISION> = zeros_ND<2, PRISM_FLOAT_PRECISION>({{dimj, dimi}});
		}

	std::complex<PRISM_FLOAT_PRECISION> answer;

//	cudaMemcpy(&answer, psi_d+(ax),1*sizeof(PRISM_CUDA_COMPLEX_FLOAT),cudaMemcpyDeviceToHost);
////	cout << " answer = " << answer << endl;
		if (ax==0 && ay==0) {
			for (auto j = 0; j < 10; ++j) {
				cudaMemcpy(&answer, psi_d + j, 1 * sizeof(PRISM_CUDA_COMPLEX_FLOAT), cudaMemcpyDeviceToHost);
				cout << " psi_d = " << answer << endl; //<< "xp = " << xp << "yp = " << yp << endl;
			}
		}

		PRISM_FLOAT_PRECISION answer2;

//	cudaMemcpy(&answer, psi_d+(ax),1*sizeof(PRISM_CUDA_COMPLEX_FLOAT),cudaMemcpyDeviceToHost);
////	cout << " answer = " << answer << endl;
		if (ax==0 && ay==0) {
			for (auto j = 0; j < 10; ++j) {
				cudaMemcpy(&answer2, psi_abssquared_d + j, 1 * sizeof(PRISM_FLOAT_PRECISION), cudaMemcpyDeviceToHost);
				cout << " psi_abssquared_d = " << answer2 << endl << "xp = " << xp << "yp = " << yp << endl;
			}
		}

		std::complex<PRISM_FLOAT_PRECISION> answer_propd;

//	cudaMemcpy(&answer, psi_d+(ax),1*sizeof(PRISM_CUDA_COMPLEX_FLOAT),cudaMemcpyDeviceToHost);
////	cout << " answer = " << answer << endl;
		if (ax==0 && ay==0) {
			for (auto j = 0; j < 10; ++j) {
				cudaMemcpy(&answer_propd, prop_d + j, 1 * sizeof(std::complex<PRISM_FLOAT_PRECISION>), cudaMemcpyDeviceToHost);
				cout << " prop_d = " << answer_propd << endl;
			}
		}

		std::complex<PRISM_FLOAT_PRECISION> answer_trans_d;

//	cudaMemcpy(&answer, psi_d+(ax),1*sizeof(PRISM_CUDA_COMPLEX_FLOAT),cudaMemcpyDeviceToHost);
////	cout << " answer = " << answer << endl;
		if (ax==0 && ay==0) {
			for (auto j = 0; j < 10; ++j) {
				cudaMemcpy(&answer_trans_d, trans_d + j, 1 * sizeof(std::complex<PRISM_FLOAT_PRECISION>), cudaMemcpyDeviceToHost);
				cout << " answer_trans_d = " << answer_trans_d << endl;
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
	cufftDestroy(plan);
	cudaFree(psi_d);
	cudaStreamDestroy(stream2);
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
		PRISM_CUDA_COMPLEX_FLOAT *prop_d[pars.meta.NUM_GPUS];
		PRISM_FLOAT_PRECISION    *qxa_d[pars.meta.NUM_GPUS];
		PRISM_FLOAT_PRECISION    *qya_d[pars.meta.NUM_GPUS];

		// allocate memory on each GPU
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			cudaSetDevice(g);
			cudaMalloc((void **) &PsiProbeInit_d[g], PsiProbeInit.size() * sizeof(PsiProbeInit[0]));
			cudaMalloc((void **) &trans_d[g], trans.size() * sizeof(trans[0]));
			cudaMalloc((void **) &qxa_d[g], pars.qxa.size() * sizeof(pars.qxa[0]));
			cudaMalloc((void **) &qya_d[g], pars.qya.size() * sizeof(pars.qya[0]));
			cudaMalloc((void **) &prop_d[g], pars.prop.size() * sizeof(pars.prop[0]));
		}

		// copy memory to each GPU (this can be made asynchronous if necessary by copying to pinned memory first)
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			cudaSetDevice(g);
			cudaMemcpy(PsiProbeInit_d[g], &PsiProbeInit[0], PsiProbeInit.size() * sizeof(PsiProbeInit[0]), cudaMemcpyHostToDevice);
			cudaMemcpy(trans_d[g], &trans[0], trans.size() * sizeof(trans[0]), cudaMemcpyHostToDevice);
			cudaMemcpy(qxa_d[g], &pars.qxa[0], pars.qxa.size() * sizeof(pars.qxa[0]), cudaMemcpyHostToDevice);
			cudaMemcpy(qya_d[g], &pars.qya[0], pars.qya.size() * sizeof(pars.qya[0]), cudaMemcpyHostToDevice);
			cudaMemcpy(prop_d[g], &pars.prop[0], pars.prop.size() * sizeof(pars.prop[0]), cudaMemcpyHostToDevice);
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
			PRISM_CUDA_COMPLEX_FLOAT *current_prop_d  = prop_d[gpu_num];
			PRISM_FLOAT_PRECISION *current_qxa_d      = qxa_d[gpu_num];
			PRISM_FLOAT_PRECISION *current_qya_d      = qya_d[gpu_num];

			// launch a new thread
			// emplace_back is better whenever constructing a new object
			workers_gpu.emplace_back(thread([&pars, current_trans_d, current_PsiProbeInit_d, &alphaInd,
					                                start, stop, gpu_num, current_qya_d, current_qxa_d,
					                                current_prop_d, &current_stream, &psi_size, &PsiProbeInit]() {

				// page-locked (pinned) memory for async streaming of the result back
				PRISM_FLOAT_PRECISION *pinned_output;

				// figure out how much pinned memory to allocate in this job
				size_t num_probe_positions = pars.xp.size() * (min((size_t) stop, pars.yp.size()) - start);
				size_t pinned_output_size = psi_size * num_probe_positions;

				// allocate pinned memory on host
				cudaMallocHost((void **) &pinned_output, pinned_output_size * sizeof(PRISM_FLOAT_PRECISION));

				// set the GPU context
				cudaSetDevice(gpu_num); // set current gpu
				PRISM_FLOAT_PRECISION *pinned_output_begin = pinned_output; // pointer to the beginning of corresponding output layer in the 3D array
				for (auto ay = start; ay < std::min((size_t) stop, pars.yp.size()); ++ay) {
					for (auto ax = 0; ax < pars.xp.size(); ++ax) {
//				for (auto ay = start; ay < 1; ++ay) {
//					for (auto ax = 0; ax < 1; ++ax) {
						getMultisliceProbe_gpu(pars, current_trans_d, current_PsiProbeInit_d, current_qya_d,
						                       current_qxa_d,
						                       current_prop_d, ay, ax, PsiProbeInit.get_dimj(), PsiProbeInit.get_dimi(),
						                       alphaInd, current_stream, pinned_output_begin);
						pinned_output_begin += psi_size; // advance the start point of the output
					}
				}

				cudaDeviceSynchronize();
//				{
//					for (auto j = 0; j < 25; ++j)cout << "pinned_output[j] = " << pinned_output[j] << endl;
//				}
				PRISM_FLOAT_PRECISION c = 0;
				for (auto i = 0; i < pinned_output_size; ++i) {
					c += pinned_output[i];
				}
				cout << "c = " << c << endl;
				PRISM_FLOAT_PRECISION *counts = pinned_output;

				Array2D<PRISM_FLOAT_PRECISION> db = zeros_ND<2, PRISM_FLOAT_PRECISION>({{PsiProbeInit.get_dimj(), PsiProbeInit.get_dimi()}});
				auto db_ptr = db.begin();
				for (auto ay = start; ay < std::min((size_t) stop, pars.yp.size()); ++ay) {
					for (auto ax = 0; ax < pars.xp.size(); ++ax) {
						auto idx = alphaInd.begin();
						while (idx != alphaInd.end()) {
							if (ay==0 && ax==0)*db_ptr = *counts;
							if (*idx <= pars.Ndet) {
//								cout << "count = " << *counts << endl;
								pars.stack.at(ay, ax, (*idx) - 1, 0) += (*counts);

//								cout << "ax = " << ax << endl;
//								cout << "ay = " << ay << endl;
//								cout << "(*idx) - 1 = " << (*idx) - 1 << endl;
//								cout << "pars.stack.at(ay, ax, (*idx) - 1, 0) = "
//								     << pars.stack.at(ay, ax, (*idx) - 1, 0)
//								     << endl;
							}
							db_ptr++;
							++idx;
							++counts;
						}
					}
				}
				db.toMRC_f("db_intOutput.mrc");
			cout << " DEBUG 1 " << endl;
				if (start == 0 ) {
					Array2D <PRISM_FLOAT_PRECISION> prism_image;
					prism_image = zeros_ND<2, PRISM_FLOAT_PRECISION>({{pars.stack.get_diml(), pars.stack.get_dimk()}});
					for (auto y = 0; y < pars.stack.get_diml(); ++y) {
						for (auto x = 0; x < pars.stack.get_dimk(); ++x) {
							for (auto b = 13; b < 18; ++b) {
								prism_image.at(y, x) += pars.stack.at(y, x, b, 0);
//								cout << "prism_image.at(y, x) = " << prism_image.at(y, x) << endl;
							}
						}
					}
					prism_image.toMRC_f("TEST.mrc");
					cout <<" debug written" <<endl;
				}
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
				cudaFreeHost(pinned_output);
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
			cudaFree(prop_d[j]);
		}
		// destroy CUDA streams
		for (auto j = 0; j < total_num_streams; ++j)cudaStreamDestroy(streams[j]);
//		cudaDeviceReset();
	}
}
