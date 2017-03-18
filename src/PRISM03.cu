// For variable naming, the suffixes are "_d" for "device" (1 copy per GPU), "_ds" for "device stream (1 copy per stream), "_ph" for "pinned host"


#include <iostream>
#include <complex>
#include <thread>
#include <vector>
#include "getWorkID.h"
#include "PRISM03.cuh"
#include "PRISM03.h"
#include "configure.h"
#include "ArrayND.h"
#include "params.h"
#include "utility.cuh"

namespace PRISM {
	// define some constants
	__device__ __constant__ float pi_f                  = PI;
	__device__ __constant__ cuFloatComplex i_f          = {0, 1};
	__device__ __constant__ cuFloatComplex pi_cx_f      = {PI, 0};
	__device__ __constant__ cuFloatComplex minus_2pii_f = {0, -2*PI};
	__device__ __constant__ double pi                   = PI;
	__device__ __constant__ cuDoubleComplex i           = {0, 1};
	__device__ __constant__ cuDoubleComplex pi_cx       = {PI, 0};
	__device__ __constant__ cuDoubleComplex minus_2pii  = {0, -2*PI};

template <size_t BlockSizeX>
__global__ void scaleReduceS(const PRISM_CUDA_COMPLEX_FLOAT *permuted_Scompact_d,
                             const PRISM_CUDA_COMPLEX_FLOAT *phaseCoeffs_ds,
                             PRISM_CUDA_COMPLEX_FLOAT *psi_ds,
                             const long *z_ds,
                             const long* y_ds,
                             const size_t numberBeams,
                             const size_t dimk_S,
                             const size_t dimj_S,
                             const size_t dimj_psi,
                             const size_t dimi_psi) {
	// This code is heavily modeled after Mark Harris's presentation on optimized parallel reduction
	// http://developer.download.nvidia.com/compute/cuda/1.1-Beta/x86_website/projects/reduction/doc/reduction.pdf

	// for the permuted Scompact matrix, the x direction runs along the number of beams, leaving y and z to represent the
//	 2D array of reduced values in psi
	volatile __shared__ cuFloatComplex scaled_values[BlockSizeX];
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	int y = blockIdx.y;
	int z = blockIdx.z;

	scaled_values[idx] = make_cuFloatComplex(0,0); // guarantee the memory is initialized to 0 so we can accumulate without bounds checking
	__syncthreads();
////	 read in first values
	if (idx < numberBeams) {
		scaled_values[idx] = cuCmulf(permuted_Scompact_d[z_ds[z]*numberBeams*dimj_S + y_ds[y]*numberBeams + idx],
		                             phaseCoeffs_ds[idx]);
		__syncthreads();
	}

	// step through global memory accumulating until values have been reduced to BlockSizeX elements in shared memory
	size_t offset = BlockSizeX;
	while (offset < numberBeams){
		if (idx + offset < numberBeams){
			scaled_values[idx] = cuCaddf(scaled_values[idx],
										 cuCmulf(  permuted_Scompact_d[z_ds[z]*numberBeams*dimj_S + y_ds[y]*numberBeams + idx + offset],
									              phaseCoeffs_ds[idx + offset]));
		}
		offset += BlockSizeX;
		__syncthreads(); // maybe can skip syncing every iteration and just do once at end
	}

	// At this point we have exactly BlockSizeX elements to reduce from shared memory which we will add by recursively
	// dividing the array in half

	// Take advantage of templates. Because BlockSizeX is passed at compile time, all of these comparisons are also
	// evaluated at compile time
	if (BlockSizeX >= 1024){
		if (idx < 512){
			scaled_values[idx] = cuCaddf(scaled_values[idx], scaled_values[idx + 512]);
		}
		__syncthreads();
	}

	if (BlockSizeX >= 512){
		if (idx < 256){
			scaled_values[idx] = cuCaddf(scaled_values[idx], scaled_values[idx + 256]);
		}
		__syncthreads();
	}

	if (BlockSizeX >= 256){
		if (idx < 128){
			scaled_values[idx] = cuCaddf(scaled_values[idx], scaled_values[idx + 128]);
		}
		__syncthreads();
	}

	if (BlockSizeX >= 128){
		if (idx < 64){
			scaled_values[idx] = cuCaddf(scaled_values[idx], scaled_values[idx + 64]);
		}
		__syncthreads();
	}

	// Now we have reduced to 64 elements, and thus the remaining work need only be done by a single warp of 32 threads.
	// This means we can completely unroll the remaining loops with no synchronization
//
//	if (idx < 32){
//		if (256 >= 64)scaled_values[idx] = cuCaddf(scaled_values[idx], scaled_values[idx + 32]);
//		if (256 >= 32)scaled_values[idx] = cuCaddf(scaled_values[idx], scaled_values[idx + 16]);
//		if (256 >= 16)scaled_values[idx] = cuCaddf(scaled_values[idx], scaled_values[idx + 8]);
//		if (256 >= 8) scaled_values[idx] = cuCaddf(scaled_values[idx], scaled_values[idx + 4]);
//		if (256 >= 4) scaled_values[idx] = cuCaddf(scaled_values[idx], scaled_values[idx + 2]);
//		if (256 >= 2) scaled_values[idx] = cuCaddf(scaled_values[idx], scaled_values[idx + 1]);

//	if (idx < 32){
//		if (BlockSizeX >= 64){
//			scaled_values[idx].x +=  scaled_values[idx + 32].x;
//			scaled_values[idx].y +=  scaled_values[idx + 32].y;
//		}
//		if (BlockSizeX >= 32){
//			scaled_values[idx].x +=  scaled_values[idx + 16].x;
//			scaled_values[idx].y +=  scaled_values[idx + 16].y;
//		}
//		if (BlockSizeX >= 16){
//			scaled_values[idx].x +=  scaled_values[idx + 8].x;
//			scaled_values[idx].y +=  scaled_values[idx + 8].y;
//		}
//		if (BlockSizeX >= 8){
//			scaled_values[idx].x +=  scaled_values[idx + 4].x;
//			scaled_values[idx].y +=  scaled_values[idx + 4].y;
//		}
//		if (BlockSizeX >= 4){
//			scaled_values[idx].x +=  scaled_values[idx + 2].x;
//			scaled_values[idx].y +=  scaled_values[idx + 2].y;
//		}
//		if (BlockSizeX >= 2){
//			scaled_values[idx].x +=  scaled_values[idx + 1].x;
//			scaled_values[idx].y +=  scaled_values[idx + 1].y;
//		}
//	}

//	 TODO: figure out why I need syncthreads here... you shouldn't hjave to synchronize at warp level
	if (idx < 32){
		if (BlockSizeX >= 64)scaled_values[idx] = cuCaddf(scaled_values[idx], scaled_values[idx + 32]);		__syncthreads();
		if (BlockSizeX >= 32)scaled_values[idx] = cuCaddf(scaled_values[idx], scaled_values[idx + 16]);		__syncthreads();
		if (BlockSizeX >= 16)scaled_values[idx] = cuCaddf(scaled_values[idx], scaled_values[idx + 8]);		__syncthreads();
		if (BlockSizeX >= 8) scaled_values[idx] = cuCaddf(scaled_values[idx], scaled_values[idx + 4]);		__syncthreads();
		if (BlockSizeX >= 4) scaled_values[idx] = cuCaddf(scaled_values[idx], scaled_values[idx + 2]);		__syncthreads();
		if (BlockSizeX >= 2) scaled_values[idx] = cuCaddf(scaled_values[idx], scaled_values[idx + 1]);		__syncthreads();
	}
//	if (idx < 32){
//		if (BlockSizeX >= 64)scaled_values[idx] = cuCaddf(scaled_values[idx], scaled_values[idx + 32]);
//		if (BlockSizeX >= 32)scaled_values[idx] = cuCaddf(scaled_values[idx], scaled_values[idx + 16]);
//		if (BlockSizeX >= 16)scaled_values[idx] = cuCaddf(scaled_values[idx], scaled_values[idx + 8]);
//		if (BlockSizeX >= 8) scaled_values[idx] = cuCaddf(scaled_values[idx], scaled_values[idx + 4]);
//		if (BlockSizeX >= 4) scaled_values[idx] = cuCaddf(scaled_values[idx], scaled_values[idx + 2]);
//		if (BlockSizeX >= 2) scaled_values[idx] = cuCaddf(scaled_values[idx], scaled_values[idx + 1]);
//	}


	// write out the result
	if (idx == 0)psi_ds[z*dimi_psi + y] = scaled_values[0];
	//	if (idx == 0)psi_ds[z*dimi_psi + y] = make_cuFloatComplex(2,3);
}
	using namespace std;
//	__global__ void computePhaseCoeffs(cuFloatComplex* phaseCoeffs,
//	                                   const cuFloatComplex *PsiProbeInit_d,
//	                                   const float * qyaReduce_d,
//	                                   const float * qxaReduce_d,
//	                                   const size_t *yBeams_d,
//	                                   const size_t *xBeams_d,
//	                                   const float yp,
//	                                   const float xp,
//	                                   const float yTiltShift,
//	                                   const float xTiltShift,
//	                                   const size_t dimi,
//	                                   const size_t numBeams){
//		int idx = threadIdx.x + blockDim.x*blockIdx.x;
//		if (idx < numBeams) {
//			size_t yB = yBeams_d[idx];
//			size_t xB = xBeams_d[idx];
//			cuFloatComplex xp_cx = make_cuFloatComplex(xp, 0);
//			cuFloatComplex yp_cx = make_cuFloatComplex(yp, 0);
//			cuFloatComplex xTiltShift_cx = make_cuFloatComplex(xTiltShift, 0);
//			cuFloatComplex yTiltShift_cx = make_cuFloatComplex(yTiltShift, 0);
//			cuFloatComplex qya = make_cuFloatComplex(qyaReduce_d[yB * dimi + xB], 0);
//			cuFloatComplex qxa = make_cuFloatComplex(qxaReduce_d[yB * dimi + xB], 0);
//			cuFloatComplex arg1 = cuCmulf(qxa, cuCaddf(xp_cx, xTiltShift_cx));
//			cuFloatComplex arg2 = cuCmulf(qya, cuCaddf(yp_cx, yTiltShift_cx));
//			cuFloatComplex arg = cuCaddf(arg1, arg2);
//			cuFloatComplex phase_shift = exp_cx(cuCmulf(minus_2pii_f, arg));
////		phaseCoeffs[idx] = phase_shift;
//			phaseCoeffs[idx] = cuCmulf(phase_shift, PsiProbeInit_d[yB * dimi + xB]);
////		phaseCoeffs[idx] = make_cuFloatComplex(1,2);
//		}
//	}
//
//	__global__ void scaleReduceS(const PRISM_CUDA_COMPLEX_FLOAT *permuted_Scompact_d,
//	                             const PRISM_CUDA_COMPLEX_FLOAT *phaseCoeffs_ds,
//	                             PRISM_CUDA_COMPLEX_FLOAT *psi_ds,
//	                             const size_t numberBeams,
//	                             const size_t dimk_S,
//	                             const size_t dimj_S,
//	                             const size_t dimj_psi,
//	                             const size_t dimi_psi) {
//		// for the permuted Scompact matrix, the x direction runs along the number of beams, leaving y and z to represent the
//		// 2D array of reduced values in psi
//		extern __shared__ cuFloatComplex scaled_values[];
//		int idx = threadIdx.x + blockDim.x * blockIdx.x;
//
//		if (idx < numberBeams) {
//			int y = blockIdx.y;
//			int z = blockIdx.z;
//			scaled_values[idx] = cuCmulf(permuted_Scompact_d[z*numberBeams*dimj_S + y*numberBeams + idx], phaseCoeffs_ds[idx]);
//			__syncthreads();
//			if (idx == 0){
//				psi_ds[z*dimi_psi + y] = scaled_values[0];
//			}
//		}
//
//	}
	void buildPRISMOutput_GPU(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                          const PRISM_FLOAT_PRECISION xTiltShift,
	                          const PRISM_FLOAT_PRECISION yTiltShift,
	                          const Array2D<PRISM_FLOAT_PRECISION> &alphaInd,
	                          const Array2D<std::complex<PRISM_FLOAT_PRECISION> > &PsiProbeInit) {

		// create CUDA streams
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;
		cudaStream_t streams[total_num_streams];
		cufftHandle cufft_plan[total_num_streams];

//		cout << "pars.imageSizeReduce[0] = " << pars.imageSizeReduce[0] << endl;
//		cout << "pars.imageSizeReduce[1] = " << pars.imageSizeReduce[1] << endl;
		for (auto j = 0; j < total_num_streams; ++j) {
			cudaSetDevice(j % pars.meta.NUM_GPUS);
			cudaErrchk(cudaStreamCreate(&streams[j]));
			cufftErrchk(cufftPlan2d(&cufft_plan[j], pars.imageSizeReduce[1], pars.imageSizeReduce[0], PRISM_CUFFT_PLAN_TYPE));
			cufftErrchk(cufftSetStream(cufft_plan[j], streams[j]));
		}

		// pointers to pinned host memory for async transfers
		PRISM_FLOAT_PRECISION               *output_ph[total_num_streams]; // one output array per stream
		std::complex<PRISM_FLOAT_PRECISION> *permuted_Scompact_ph; // see below for explanation of why this is permuted
		std::complex<PRISM_FLOAT_PRECISION> *PsiProbeInit_ph;
//		PRISM_FLOAT_PRECISION               *xVec_ph;
//		PRISM_FLOAT_PRECISION               *yVec_ph;
		PRISM_FLOAT_PRECISION               *qxaReduce_ph;
		PRISM_FLOAT_PRECISION               *qyaReduce_ph;
		PRISM_FLOAT_PRECISION               *alphaInd_ph;
		size_t                              *xBeams_ph;
		size_t                              *yBeams_ph;


		// allocate pinned memory
		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaMallocHost((void **) &output_ph[s],
			                          pars.stack.get_dimj() * pars.stack.get_dimi() *
			                          sizeof(PRISM_FLOAT_PRECISION)));
		}
		cudaErrchk(cudaMallocHost((void **) &permuted_Scompact_ph, pars.Scompact.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
		cudaErrchk(cudaMallocHost((void **) &PsiProbeInit_ph, PsiProbeInit.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
//		cudaErrchk(cudaMallocHost((void **) &xVec_ph, pars.xVec.size() * sizeof(PRISM_FLOAT_PRECISION)));
//		cudaErrchk(cudaMallocHost((void **) &yVec_ph, pars.xVec.size() * sizeof(PRISM_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **) &qxaReduce_ph, pars.qxaReduce.size() * sizeof(PRISM_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **) &qyaReduce_ph, pars.qyaReduce.size() * sizeof(PRISM_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **) &alphaInd_ph,  alphaInd.size()       * sizeof(PRISM_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **) &xBeams_ph, pars.xyBeams.get_dimj()  * sizeof(size_t)));
		cudaErrchk(cudaMallocHost((void **) &yBeams_ph, pars.xyBeams.get_dimj()  * sizeof(size_t)));

//		cout << "NUMBEAMS = " << pars.numberBeams << endl;

		// copy host memory to pinned
		for (auto s = 0; s < total_num_streams; ++s) {
			memset(output_ph[s], 0, pars.stack.get_dimj() * pars.stack.get_dimi() *
			                        sizeof(PRISM_FLOAT_PRECISION));
		}
//		for (auto ii=0; ii < pars.detectorAngles.size(); ++ii){
//			std::cout << "first output_ph[" << ii <<"] = " << output_ph[0][ii] << std::endl;
//		}
		// the GPU computational model operates on Scompact in a different order than the CPU, and it is
		// more optimal to permute the dimensions so that the consecutive elements represent different
		// beams on the GPU as opposed to consecutive x-probe positions on the CPU
		{
			auto S_ptr = permuted_Scompact_ph;
			for (auto jj = 0; jj < pars.Scompact.get_dimj(); ++jj){
				for (auto ii = 0; ii < pars.Scompact.get_dimi(); ++ii){
					for (auto kk = 0; kk < pars.Scompact.get_dimk(); ++kk){
						*S_ptr++ = pars.Scompact.at(kk, jj, ii);
					}
				}
			}
		}
		memcpy(PsiProbeInit_ph, &(*PsiProbeInit.begin()), PsiProbeInit.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>));
//		memcpy(xVec_ph, &(*pars.xVec.begin()), pars.xVec.size() * sizeof(PRISM_FLOAT_PRECISION));
//		memcpy(yVec_ph, &(*pars.yVec.begin()), pars.yVec.size() * sizeof(PRISM_FLOAT_PRECISION));
		memcpy(alphaInd_ph,  &(*alphaInd.begin()),       alphaInd.size()       * sizeof(PRISM_FLOAT_PRECISION));
		memcpy(qxaReduce_ph, &(*pars.qxaReduce.begin()), pars.qxaReduce.size() * sizeof(PRISM_FLOAT_PRECISION));
		memcpy(qyaReduce_ph, &(*pars.qyaReduce.begin()), pars.qyaReduce.size() * sizeof(PRISM_FLOAT_PRECISION));

		{
			auto x_ptr = xBeams_ph;
			auto y_ptr = yBeams_ph;
			for (auto jj = 0; jj < pars.xyBeams.get_dimj(); ++jj){
				*y_ptr++ = pars.xyBeams.at(jj,0);
				*x_ptr++ = pars.xyBeams.at(jj,1);
			}
		}


//		for (auto i = 0; i < 10; ++i){
//			cout << "PsiProbeInit[" << i << "] = " << PsiProbeInit[i] << endl;
//			cout << "PsiProbeInit_ph[" << i << "] = " << PsiProbeInit_ph[i] << endl;
//		}
//		for (auto i = 0; i < 10; ++i){
//			cout << "permuted_Scompact_ph[" << i << "] = " << permuted_Scompact_ph[i] << endl;
//			cout << "pars.Scompact.at("  << i << ",0,0) = " << pars.Scompact.at(i,0,0) << endl;
//		}

		// pointers to read-only GPU memory (one copy per GPU)
		PRISM_CUDA_COMPLEX_FLOAT *permuted_Scompact_d[pars.meta.NUM_GPUS];
		PRISM_CUDA_COMPLEX_FLOAT *PsiProbeInit_d[pars.meta.NUM_GPUS];
		PRISM_FLOAT_PRECISION    *qxaReduce_d[pars.meta.NUM_GPUS];
		PRISM_FLOAT_PRECISION    *qyaReduce_d[pars.meta.NUM_GPUS];
		PRISM_FLOAT_PRECISION    *alphaInd_d[pars.meta.NUM_GPUS];
		size_t                   *yBeams_d[pars.meta.NUM_GPUS];
		size_t                   *xBeams_d[pars.meta.NUM_GPUS];
//
		// pointers to read/write GPU memory (one per stream)
		PRISM_CUDA_COMPLEX_FLOAT *psi_ds[total_num_streams];
		PRISM_CUDA_COMPLEX_FLOAT *phaseCoeffs_ds[total_num_streams];
		PRISM_FLOAT_PRECISION    *psi_intensity_ds[total_num_streams];
		PRISM_FLOAT_PRECISION    *integratedOutput_ds[total_num_streams];
		long                     *y_ds[total_num_streams];
		long                     *x_ds[total_num_streams];
//
		// allocate memory on each GPU
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaMalloc((void **) &permuted_Scompact_d[g], pars.Scompact.size() * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMalloc((void **) &PsiProbeInit_d[g],      PsiProbeInit.size()  * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMalloc((void **) &qxaReduce_d[g], pars.qxaReduce.size() * sizeof(PRISM_FLOAT_PRECISION)));
			cudaErrchk(cudaMalloc((void **) &qyaReduce_d[g], pars.qyaReduce.size() * sizeof(PRISM_FLOAT_PRECISION)));
			cudaErrchk(cudaMalloc((void **) &alphaInd_d[g],  alphaInd.size()       * sizeof(PRISM_FLOAT_PRECISION)));
			cudaErrchk(cudaMalloc((void **) &yBeams_d[g], pars.xyBeams.get_dimj()  * sizeof(size_t)));
			cudaErrchk(cudaMalloc((void **) &xBeams_d[g], pars.xyBeams.get_dimj()  * sizeof(size_t)));
		}

		// allocate memory per stream and 0 it
		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaSetDevice(s % pars.meta.NUM_GPUS));
			cudaErrchk(cudaMalloc((void **) &psi_ds[s],
			                      pars.imageSizeReduce[0] * pars.imageSizeReduce[1] * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMalloc((void **) &phaseCoeffs_ds[s],
			                      pars.numberBeams * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMalloc((void **) &psi_intensity_ds[s],
			                      pars.imageSizeReduce[0] * pars.imageSizeReduce[1] * sizeof(PRISM_FLOAT_PRECISION)));
			cudaErrchk(cudaMalloc((void **) &y_ds[s],
			                      pars.xyBeams.get_dimj() * sizeof(long)));
			cudaErrchk(cudaMalloc((void **) &x_ds[s],
			                      pars.xyBeams.get_dimj() * sizeof(long)));
			cudaErrchk(cudaMalloc((void **) &integratedOutput_ds[s],
			                      pars.detectorAngles.size() * sizeof(PRISM_FLOAT_PRECISION)));
			cudaErrchk(cudaMemset(psi_ds[s], 0,
			                      pars.imageSizeReduce[0] * pars.imageSizeReduce[1] * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMemset(phaseCoeffs_ds[s], 0,
			                      pars.numberBeams * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMemset(psi_intensity_ds[s], 0,
			                      pars.imageSizeReduce[0] * pars.imageSizeReduce[1] * sizeof(PRISM_FLOAT_PRECISION)));
			cudaErrchk(cudaMemset(y_ds[s], 0,
			                      pars.xyBeams.get_dimj() * sizeof(long)));
			cudaErrchk(cudaMemset(x_ds[s], 0,
			                      pars.xyBeams.get_dimj() * sizeof(long)));
			cudaErrchk(cudaMemset(integratedOutput_ds[s], 0,
			                      pars.detectorAngles.size() * sizeof(PRISM_FLOAT_PRECISION)));
		}

		// Copy memory to each GPU asynchronously from the pinned host memory spaces.
		// The streams are laid out so that consecutive streams represent different GPUs. If we
		// have more than one stream per GPU, then we want to interleave as much as possible
		int stream_id = 0;
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			stream_id = g;
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaMemcpyAsync(permuted_Scompact_d[g], &permuted_Scompact_ph[0],
			                           pars.Scompact.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>),
			                           cudaMemcpyHostToDevice, streams[stream_id]));

			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cout << "stream_id = " << stream_id << endl;
			cudaErrchk(cudaMemcpyAsync(PsiProbeInit_d[g], &PsiProbeInit_ph[0],
			                           PsiProbeInit.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>),
			                           cudaMemcpyHostToDevice, streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cout << "stream_id = " << stream_id << endl;
			cudaErrchk(cudaMemcpyAsync(qxaReduce_d[g], &qxaReduce_ph[0],
			                           pars.qxaReduce.size() * sizeof(PRISM_FLOAT_PRECISION), cudaMemcpyHostToDevice, streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cout << "stream_id = " << stream_id << endl;
			cudaErrchk(cudaMemcpyAsync(qyaReduce_d[g], &qyaReduce_ph[0],
			                           pars.qyaReduce.size() * sizeof(PRISM_FLOAT_PRECISION), cudaMemcpyHostToDevice, streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(alphaInd_d[g], &alphaInd_ph[0],
			                           alphaInd.size() * sizeof(alphaInd[0]), cudaMemcpyHostToDevice, streams[stream_id]));
			cout << "stream_id = " << stream_id << endl;
			cudaErrchk(cudaMemcpyAsync(yBeams_d[g], &yBeams_ph[0],
			                           pars.xyBeams.get_dimj() * sizeof(size_t), cudaMemcpyHostToDevice,
			                           streams[stream_id]));
			cudaErrchk(cudaMemcpyAsync(xBeams_d[g], &xBeams_ph[0],
			                           pars.xyBeams.get_dimj() * sizeof(size_t), cudaMemcpyHostToDevice,
			                           streams[stream_id]));
		}

//		for (auto i : alphaInd)cout << "alphaInd = " << i << endl;
		// make sure transfers are complete
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			cudaSetDevice(g);
			cudaDeviceSynchronize();
		}


		// launch threads to compute results for batches of xp, yp
		// I do this by dividing the xp points among threads, and each computes
		// all of the relevant yp for each of its xp. This seems an okay strategy
		// as long as the number of xp and yp are similar. If that is not the case
		// this may need to be adapted
		vector<thread> workers_GPU;
		workers_GPU.reserve(pars.meta.NUM_THREADS); // prevents multiple reallocations
		int stream_count = 0;
		setWorkStartStop(0, pars.xp.size() * pars.yp.size());
//		setWorkStartStop(0, 1);
		for (auto t = 0; t < total_num_streams; ++t) {

			int GPU_num = stream_count % pars.meta.NUM_GPUS; // determine which GPU handles this job
			cudaSetDevice(GPU_num);
			cudaStream_t &current_stream = streams[stream_count];
			cout << "Launching GPU worker on stream #" << stream_count << " of GPU #" << GPU_num << '\n';


//			PRISM_CUDA_COMPLEX_FLOAT *psi_ds[total_num_streams];
//			PRISM_CUDA_COMPLEX_FLOAT *phaseCoeffs_ds[total_num_streams];
//			PRISM_FLOAT_PRECISION    *psi_intensity_ds[total_num_streams];
//			size_t                   *y_ds[total_num_streams];
//			size_t                   *x_ds[total_num_streams];
			// get pointers to the pre-copied arrays, making sure to get those on the current GPU
			PRISM_CUDA_COMPLEX_FLOAT *current_permuted_Scompact_d = permuted_Scompact_d[GPU_num];
			PRISM_CUDA_COMPLEX_FLOAT *current_PsiProbeInit_d      = PsiProbeInit_d[GPU_num];
			PRISM_FLOAT_PRECISION *current_qxaReduce_d            = qxaReduce_d[GPU_num];
			PRISM_FLOAT_PRECISION *current_qyaReduce_d            = qyaReduce_d[GPU_num];
			size_t *current_yBeams_d                              = yBeams_d[GPU_num];
			size_t *current_xBeams_d                              = xBeams_d[GPU_num];
			PRISM_FLOAT_PRECISION *current_alphaInd_d             = alphaInd_d[GPU_num];

			// get pointers to per-stream arrays
			PRISM_CUDA_COMPLEX_FLOAT *current_psi_ds              = psi_ds[stream_count];
			PRISM_CUDA_COMPLEX_FLOAT *current_phaseCoeffs_ds      = phaseCoeffs_ds[stream_count];
			PRISM_FLOAT_PRECISION *current_psi_intensity_ds       = psi_intensity_ds[stream_count];
			PRISM_FLOAT_PRECISION *current_integratedOutput_ds    = integratedOutput_ds[stream_count];
			long *current_y_ds                                    = y_ds[stream_count];
			long *current_x_ds                                    = x_ds[stream_count];
			cufftHandle &current_cufft_plan                       = cufft_plan[stream_count];

			// get pointer to output pinned memory
			PRISM_FLOAT_PRECISION *current_output_ph              = output_ph[stream_count];


			// emplace_back is better whenever constructing a new object
			workers_GPU.emplace_back(thread([&pars, GPU_num, stream_count, &yTiltShift, &xTiltShift, current_permuted_Scompact_d,
					                                current_alphaInd_d, current_PsiProbeInit_d, current_qxaReduce_d, current_qyaReduce_d,
					                                current_yBeams_d, current_xBeams_d, current_psi_ds, current_phaseCoeffs_ds,
					                                current_psi_intensity_ds, current_y_ds, current_x_ds, current_integratedOutput_ds,
					                                current_output_ph, &current_cufft_plan, &current_stream]() {
				size_t Nstart, Nstop, ay, ax;
				while (getWorkID(pars, Nstart, Nstop)) { // synchronously get work assignment
					while (Nstart != Nstop) {
						ay = Nstart / pars.xp.size();
						ax = Nstart % pars.xp.size();
						buildSignal_GPU(pars, ay, ax, yTiltShift, xTiltShift, current_permuted_Scompact_d,
						                current_PsiProbeInit_d, current_qxaReduce_d, current_qyaReduce_d,
						                current_yBeams_d, current_xBeams_d, current_alphaInd_d, current_psi_ds,
						                current_phaseCoeffs_ds, current_psi_intensity_ds, current_y_ds,
						                current_x_ds, current_output_ph, current_integratedOutput_ds, current_cufft_plan, current_stream );



//						buildSignal_CPU(pars, ay, ax, yTiltShift, xTiltShift, alphaInd, PsiProbeInit);
						++Nstart;
					}
				}
			}));
			++stream_count;
		}



	////// The CUDA implementation of this reduction is so much faster that it is better to not do any CPU work here
//		// Now launch CPU work
//		if (pars.meta.also_do_CPU_work) {
//			vector<thread> workers_CPU;
//			workers_CPU.reserve(pars.meta.NUM_THREADS); // prevents multiple reallocations
//			for (auto t = 0; t < pars.meta.NUM_THREADS; ++t) {
//				cout << "Launching CPU worker thread #" << t << " to compute partial PRISM result\n";
//				// emplace_back is better whenever constructing a new object
//				workers_CPU.emplace_back(thread([&pars, &xTiltShift, &yTiltShift,
//						                            &alphaInd, &PsiProbeInit, t]() {
//					size_t Nstart, Nstop, ay, ax, early_CPU_stop;
//					Nstop = 0;
//					//early_CPU_stop = pars.xp.size() * pars.yp.size() * (1-pars.meta.cpu_gpu_ratio);
//					early_CPU_stop = pars.xp.size() * pars.yp.size() * (1-pars.meta.cpu_gpu_ratio);
////					early_CPU_stop = 0;
//					while (getWorkID(pars, Nstart, Nstop)) { // synchronously get work assignment
//						while (Nstart != Nstop) {
//							ay = Nstart / pars.xp.size();
//							ax = Nstart % pars.xp.size();
//							buildSignal_CPU(pars, ay, ax, yTiltShift, xTiltShift, alphaInd, PsiProbeInit);
//							++Nstart;
//						}
//						if (Nstop >= early_CPU_stop) break;
//					}
//					cout << "CPU worker #" << t << " finished\n";
//				}));
//			}
//			cout << "Waiting for CPU threads...\n";
//			for (auto& t:workers_CPU)t.join();
//		}


		// synchronize
		cout << "Waiting for GPU threads...\n";
		for (auto &t:workers_GPU)t.join();

		// free pinned memory
		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaFreeHost(output_ph[s]));
		}
		cudaErrchk(cudaFreeHost(permuted_Scompact_ph));
		cudaErrchk(cudaFreeHost(PsiProbeInit_ph));
//		cudaErrchk(cudaFreeHost(xVec_ph));
//		cudaErrchk(cudaFreeHost(yVec_ph));
		cudaErrchk(cudaFreeHost(qxaReduce_ph));
		cudaErrchk(cudaFreeHost(qyaReduce_ph));
		cudaErrchk(cudaFreeHost(xBeams_ph));
		cudaErrchk(cudaFreeHost(yBeams_ph));

		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaFree(permuted_Scompact_d[g]));
			cudaErrchk(cudaFree(PsiProbeInit_d[g]));
			cudaErrchk(cudaFree(qxaReduce_d[g]));
			cudaErrchk(cudaFree(qyaReduce_d[g]));
			cudaErrchk(cudaFree(yBeams_d[g]));
			cudaErrchk(cudaFree(xBeams_d[g]));
			cudaErrchk(cudaFree(alphaInd_d[g]));
		}

		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaSetDevice(s % pars.meta.NUM_GPUS));
			cudaErrchk(cudaFree(psi_ds[s]));
			cudaErrchk(cudaFree(phaseCoeffs_ds[s]));
			cudaErrchk(cudaFree(psi_intensity_ds[s]));
			cudaErrchk(cudaFree(y_ds[s]));
			cudaErrchk(cudaFree(x_ds[s]));
			cudaErrchk(cudaFree(integratedOutput_ds[s]));
			cufftErrchk(cufftDestroy(cufft_plan[s]));
		}


	}

	void buildSignal_GPU(Parameters<PRISM_FLOAT_PRECISION>&  pars,
	                     const size_t& ay,
	                     const size_t& ax,
	                     const PRISM_FLOAT_PRECISION& yTiltShift,
	                     const PRISM_FLOAT_PRECISION& xTiltShift,
	                     const PRISM_CUDA_COMPLEX_FLOAT *permuted_Scompact_d,
	                     const PRISM_CUDA_COMPLEX_FLOAT *PsiProbeInit_d,
	                     const PRISM_FLOAT_PRECISION *qxaReduce_d,
	                     const PRISM_FLOAT_PRECISION *qyaReduce_d,
	                     const size_t *yBeams_d,
	                     const size_t *xBeams_d,
	                     const PRISM_FLOAT_PRECISION *alphaInd_d,
	                     PRISM_CUDA_COMPLEX_FLOAT *psi_ds,
	                     PRISM_CUDA_COMPLEX_FLOAT *phaseCoeffs_ds,
	                     PRISM_FLOAT_PRECISION *psi_intensity_ds,
	                     long  *y_ds,
	                     long  *x_ds,
	                     PRISM_FLOAT_PRECISION *output_ph,
	                     PRISM_FLOAT_PRECISION *integratedOutput_ds,
	                     const cufftHandle &cufft_plan,
	                     const cudaStream_t& stream){

//		cout << "ax = " << ax << endl;
//		cout << "ay = " << ay << endl;
		const PRISM_FLOAT_PRECISION yp = pars.yp[ay];
		const PRISM_FLOAT_PRECISION xp = pars.xp[ax];
		const size_t psi_size = pars.imageSizeReduce[0] * pars.imageSizeReduce[1];
//		cout << "xp = " << xp << endl;
//		cout << "yp = " << yp << endl;
//		cout << "pars.pixelSize[0] = " << pars.pixelSize[0]<< endl;
		shiftIndices <<<(pars.imageSizeReduce[0] - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream>>> (
			y_ds, std::round(yp / pars.pixelSizeOutput[0]),pars.imageSize[0], pars.imageSizeReduce[0]);

		shiftIndices <<<(pars.imageSizeReduce[1] - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream>>> (
			x_ds, std::round(xp / pars.pixelSizeOutput[1]), pars.imageSize[1], pars.imageSizeReduce[1]);

		computePhaseCoeffs <<<(pars.numberBeams - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream>>>(
                   phaseCoeffs_ds, PsiProbeInit_d, qyaReduce_d, qxaReduce_d,
		           yBeams_d, xBeams_d, yp, xp, yTiltShift, xTiltShift, pars.imageSizeReduce[1], pars.numberBeams);
		dim3 grid(1, pars.imageSizeReduce[0], pars.imageSizeReduce[1]);
//		cout << "pars.Scompact.get_dimk() = " << pars.Scompact.get_dimk() << endl;
//		cout << "pars.Scompact.get_dimj() = " << pars.Scompact.get_dimj() << endl;
//		cout << "pars.Scompact.get_dimi() = " << pars.Scompact.get_dimi() << endl;
//		dim3 block(32, 1, 1); // should make multiple of 32
//		scaleReduceS<32> <<< grid, block, 0, stream>>>(
//				permuted_Scompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
//						pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]);
		dim3 block(512, 1, 1); // should make multiple of 32
		scaleReduceS<512> <<< grid, block, 0, stream>>>(
				    permuted_Scompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
					pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]);
//		{
////			cudaDeviceSynchronize();
//			complex<PRISM_FLOAT_PRECISION > ans;
//			for (auto i = 0; i < pars.imageSizeReduce[0]; ++i){
//				cudaErrchk(cudaMemcpy(&ans, (psi_ds + i), sizeof(complex<PRISM_FLOAT_PRECISION>),cudaMemcpyDeviceToHost));
//				cout << "psi_ds[" << i << "] = " << ans << endl;
//			}
//		}

		cufftErrchk(PRISM_CUFFT_EXECUTE(cufft_plan, &psi_ds[0], &psi_ds[0], CUFFT_FORWARD));
		abs_squared <<< (psi_size - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >>> (psi_intensity_ds, psi_ds, psi_size);
//		cout << "pars.alphaInd"
		formatOutput_GPU_integrate(pars, psi_intensity_ds, alphaInd_d, output_ph,
		                           integratedOutput_ds, ay, ax, pars.imageSizeReduce[0],
		                           pars.imageSizeReduce[1], stream, pars.scale);


//		setAll<<< (num_integration_bins - 1)/BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream>>>(output_ph, 0, num_integration_bins);

		// Copy result. For the integration case the 4th dim of stack is 1, so the offset strides need only consider k and j
//		cudaErrchk(cudaMemcpyAsync(&stack_ph[ay*pars.stack.get_dimk()*pars.stack.get_dimj()+ ax*pars.stack.get_dimj()],integratedOutput_ds,
//		                           num_integration_bins * sizeof(PRISM_FLOAT_PRECISION),
//		                           cudaMemcpyDeviceToHost, stream));

		// wait for the copy to complete and then copy on the host. Other host threads exist doing work so this wait isn't costing anything
//		const size_t stack_start_offset = ay*pars.stack.get_dimk()*pars.stack.get_dimj()+ ax*pars.stack.get_dimj();
//		memcpy(&pars.stack[stack_start_offset], &stack_ph[stack_start_offset], num_integration_bins * sizeof(PRISM_FLOAT_PRECISION));



//		{
//			size_t ans;
//			for (auto i = 0; i < pars.imageSizeReduce[0]; ++i){
//				cudaErrchk(cudaMemcpy(&ans, (yBeams_d + i), sizeof(size_t),cudaMemcpyDeviceToHost));
//				cout << "yBeams_d[" << i << "] = " << ans << endl;
//			}
//			for (auto i = 0; i < pars.imageSizeReduce[0]; ++i){
//				cudaErrchk(cudaMemcpy(&ans, (xBeams_d + i), sizeof(size_t),cudaMemcpyDeviceToHost));
//				cout << "xBeams_d[" << i << "] = " << ans << endl;
//			}
//		}
//
//		{
//			long ans;
//			for (auto i = 0; i < pars.imageSizeReduce[0]; ++i){
//				cudaErrchk(cudaMemcpy(&ans, (y_ds + i), sizeof(long),cudaMemcpyDeviceToHost));
//				cout << "y_ds[" << i << "] = " << ans << endl;
//			}
//			for (auto i = 0; i < pars.imageSizeReduce[0]; ++i){
//				cudaErrchk(cudaMemcpy(&ans, (x_ds + i), sizeof(long),cudaMemcpyDeviceToHost));
//				cout << "x_ds[" << i << "] = " << ans << endl;
//			}
//		}



//		{
////			cudaDeviceSynchronize();
//			complex<PRISM_FLOAT_PRECISION > ans;
//			for (auto i = 0; i < pars.imageSizeReduce[0]; ++i){
//				cudaErrchk(cudaMemcpy(&ans, (phaseCoeffs_ds + i), sizeof(complex<PRISM_FLOAT_PRECISION>),cudaMemcpyDeviceToHost));
//				cout << "phaseCoeffs_ds[" << i << "] = " << ans << endl;
//			}
//		}
//

//		cout <<"test\n";
	}
}