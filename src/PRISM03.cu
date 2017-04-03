// For variable naming, the suffixes are "_d" for "device" (1 copy per GPU), "_ds" for "device stream (1 copy per stream), "_ph" for "pinned host"


#include <iostream>
#include <complex>
#include <thread>
#include <vector>
#include "WorkDispatcher.h"
#include "PRISM03.cuh"
#include "PRISM03.h"
#include "configure.h"
#include "ArrayND.h"
#include "params.h"
#include "utility.cuh"

namespace PRISM {
    extern std::mutex fftw_plan_lock;

	// define some constants
	__device__ __constant__ float pi_f                  = PI;
	__device__ __constant__ cuFloatComplex i_f          = {0, 1};
	__device__ __constant__ cuFloatComplex pi_cx_f      = {PI, 0};
	__device__ __constant__ cuFloatComplex minus_2pii_f = {0, -2*PI};
	__device__ __constant__ double pi                   = PI;
	__device__ __constant__ cuDoubleComplex i           = {0, 1};
	__device__ __constant__ cuDoubleComplex pi_cx       = {PI, 0};
	__device__ __constant__ cuDoubleComplex minus_2pii  = {0, -2*PI};



	//There are a number of template specializations that follow. They exploit the fact that scaleReduceS provides one
	// element of shared memory per thread. Therefore, for every reduction operation other than the one involving BlockSizeX/2 threads
	// the bounds checking operation can be skipped -- the latter threads just operate on meaningless values
	template <size_t BlockSizeX>
	__device__  void warpReduce_cx(volatile cuFloatComplex* sdata, int idx){
		// When 32 or fewer threads remain, there is only a single warp remaining and no need to synchronize; however,
		// the volatile keyword is necessary otherwise the compiler will optimize these operations into registers
		// and the result will be incorrect
		if (BlockSizeX >= 64){
			sdata[idx].x += sdata[idx + 32].x;
			sdata[idx].y += sdata[idx + 32].y;
		}
		if (BlockSizeX >= 32){
			sdata[idx].x += sdata[idx + 16].x;
			sdata[idx].y += sdata[idx + 16].y;
		}
		if (BlockSizeX >= 16){
			sdata[idx].x += sdata[idx + 8].x;
			sdata[idx].y += sdata[idx + 8].y;
		}
		if (BlockSizeX >= 8){
			sdata[idx].x += sdata[idx + 4].x;
			sdata[idx].y += sdata[idx + 4].y;
		}
		if (BlockSizeX >= 4){
			sdata[idx].x += sdata[idx + 2].x;
			sdata[idx].y += sdata[idx + 2].y;
		}
		if (BlockSizeX >= 2){
			sdata[idx].x += sdata[idx + 1].x;
			sdata[idx].y += sdata[idx + 1].y;
		}

	}

	template <>
	__device__  void warpReduce_cx<32>(volatile cuFloatComplex* sdata, int idx){
		// When 32 or fewer threads remain, there is only a single warp remaining and no need to synchronize; however,
		// the volatile keyword is necessary otherwise the compiler will optimize these operations into registers
		// and the result will be incorrect
			if (idx < 16) {
				sdata[idx].x += sdata[idx + 16].x;
				sdata[idx].y += sdata[idx + 16].y;
			}
			sdata[idx].x += sdata[idx + 8].x;
			sdata[idx].y += sdata[idx + 8].y;
			sdata[idx].x += sdata[idx + 4].x;
			sdata[idx].y += sdata[idx + 4].y;
			sdata[idx].x += sdata[idx + 2].x;
			sdata[idx].y += sdata[idx + 2].y;
			sdata[idx].x += sdata[idx + 1].x;
			sdata[idx].y += sdata[idx + 1].y;
	}

	template <>
	__device__  void warpReduce_cx<16>(volatile cuFloatComplex* sdata, int idx){
		// When 32 or fewer threads remain, there is only a single warp remaining and no need to synchronize; however,
		// the volatile keyword is necessary otherwise the compiler will optimize these operations into registers
		// and the result will be incorrect
		if (idx < 8) {
			sdata[idx].x += sdata[idx + 8].x;
			sdata[idx].y += sdata[idx + 8].y;
		}
		sdata[idx].x += sdata[idx + 4].x;
		sdata[idx].y += sdata[idx + 4].y;
		sdata[idx].x += sdata[idx + 2].x;
		sdata[idx].y += sdata[idx + 2].y;
		sdata[idx].x += sdata[idx + 1].x;
		sdata[idx].y += sdata[idx + 1].y;
	}

	template <>
	__device__  void warpReduce_cx<8>(volatile cuFloatComplex* sdata, int idx){
		// When 32 or fewer threads remain, there is only a single warp remaining and no need to synchronize; however,
		// the volatile keyword is necessary otherwise the compiler will optimize these operations into registers
		// and the result will be incorrect
		if (idx < 4) {
			sdata[idx].x += sdata[idx + 4].x;
			sdata[idx].y += sdata[idx + 4].y;
		}
		sdata[idx].x += sdata[idx + 2].x;
		sdata[idx].y += sdata[idx + 2].y;
		sdata[idx].x += sdata[idx + 1].x;
		sdata[idx].y += sdata[idx + 1].y;
	}


	template <>
	__device__  void warpReduce_cx<4>(volatile cuFloatComplex* sdata, int idx){
		// When 32 or fewer threads remain, there is only a single warp remaining and no need to synchronize; however,
		// the volatile keyword is necessary otherwise the compiler will optimize these operations into registers
		// and the result will be incorrect
		if (idx < 2) {
			sdata[idx].x += sdata[idx + 2].x;
			sdata[idx].y += sdata[idx + 2].y;
		}
		sdata[idx].x += sdata[idx + 1].x;
		sdata[idx].y += sdata[idx + 1].y;
	}

	template <>
	__device__  void warpReduce_cx<2>(volatile cuFloatComplex* sdata, int idx){
		// When 32 or fewer threads remain, there is only a single warp remaining and no need to synchronize; however,
		// the volatile keyword is necessary otherwise the compiler will optimize these operations into registers
		// and the result will be incorrect
		if (idx < 1) {
			sdata[idx].x += sdata[idx + 1].x;
			sdata[idx].y += sdata[idx + 1].y;
		}
	}



	template <size_t BlockSizeX>
	__device__  void warpReduce_cx(volatile cuDoubleComplex* sdata, int idx){
		// When 32 or fewer threads remain, there is only a single warp remaining and no need to synchronize; however,
		// the volatile keyword is necessary otherwise the compiler will optimize these operations into registers
		// and the result will be incorrect
		if (BlockSizeX >= 64){
			sdata[idx].x += sdata[idx + 32].x;
			sdata[idx].y += sdata[idx + 32].y;
		}
		if (BlockSizeX >= 32){
			sdata[idx].x += sdata[idx + 16].x;
			sdata[idx].y += sdata[idx + 16].y;
		}
		if (BlockSizeX >= 16){
			sdata[idx].x += sdata[idx + 8].x;
			sdata[idx].y += sdata[idx + 8].y;
		}
		if (BlockSizeX >= 8){
			sdata[idx].x += sdata[idx + 4].x;
			sdata[idx].y += sdata[idx + 4].y;
		}
		if (BlockSizeX >= 4){
			sdata[idx].x += sdata[idx + 2].x;
			sdata[idx].y += sdata[idx + 2].y;
		}
		if (BlockSizeX >= 2){
			sdata[idx].x += sdata[idx + 1].x;
			sdata[idx].y += sdata[idx + 1].y;
		}

	}

	template <>
	__device__  void warpReduce_cx<32>(volatile cuDoubleComplex* sdata, int idx){
		// When 32 or fewer threads remain, there is only a single warp remaining and no need to synchronize; however,
		// the volatile keyword is necessary otherwise the compiler will optimize these operations into registers
		// and the result will be incorrect
		if (idx < 16) {
			sdata[idx].x += sdata[idx + 16].x;
			sdata[idx].y += sdata[idx + 16].y;
		}
		sdata[idx].x += sdata[idx + 8].x;
		sdata[idx].y += sdata[idx + 8].y;
		sdata[idx].x += sdata[idx + 4].x;
		sdata[idx].y += sdata[idx + 4].y;
		sdata[idx].x += sdata[idx + 2].x;
		sdata[idx].y += sdata[idx + 2].y;
		sdata[idx].x += sdata[idx + 1].x;
		sdata[idx].y += sdata[idx + 1].y;
	}

	template <>
	__device__  void warpReduce_cx<16>(volatile cuDoubleComplex* sdata, int idx){
		// When 32 or fewer threads remain, there is only a single warp remaining and no need to synchronize; however,
		// the volatile keyword is necessary otherwise the compiler will optimize these operations into registers
		// and the result will be incorrect
		if (idx < 8) {
			sdata[idx].x += sdata[idx + 8].x;
			sdata[idx].y += sdata[idx + 8].y;
		}
		sdata[idx].x += sdata[idx + 4].x;
		sdata[idx].y += sdata[idx + 4].y;
		sdata[idx].x += sdata[idx + 2].x;
		sdata[idx].y += sdata[idx + 2].y;
		sdata[idx].x += sdata[idx + 1].x;
		sdata[idx].y += sdata[idx + 1].y;
	}

	template <>
	__device__  void warpReduce_cx<8>(volatile cuDoubleComplex* sdata, int idx){
		// When 32 or fewer threads remain, there is only a single warp remaining and no need to synchronize; however,
		// the volatile keyword is necessary otherwise the compiler will optimize these operations into registers
		// and the result will be incorrect
		if (idx < 4) {
			sdata[idx].x += sdata[idx + 4].x;
			sdata[idx].y += sdata[idx + 4].y;
		}
		sdata[idx].x += sdata[idx + 2].x;
		sdata[idx].y += sdata[idx + 2].y;
		sdata[idx].x += sdata[idx + 1].x;
		sdata[idx].y += sdata[idx + 1].y;
	}


	template <>
	__device__  void warpReduce_cx<4>(volatile cuDoubleComplex* sdata, int idx){
		// When 32 or fewer threads remain, there is only a single warp remaining and no need to synchronize; however,
		// the volatile keyword is necessary otherwise the compiler will optimize these operations into registers
		// and the result will be incorrect
		if (idx < 2) {
			sdata[idx].x += sdata[idx + 2].x;
			sdata[idx].y += sdata[idx + 2].y;
		}
		sdata[idx].x += sdata[idx + 1].x;
		sdata[idx].y += sdata[idx + 1].y;
	}

	template <>
	__device__  void warpReduce_cx<2>(volatile cuDoubleComplex* sdata, int idx){
		// When 32 or fewer threads remain, there is only a single warp remaining and no need to synchronize; however,
		// the volatile keyword is necessary otherwise the compiler will optimize these operations into registers
		// and the result will be incorrect
		if (idx < 1) {
			sdata[idx].x += sdata[idx + 1].x;
			sdata[idx].y += sdata[idx + 1].y;
		}
	}



	template <size_t BlockSizeX>
__global__ void scaleReduceS(const cuFloatComplex *permuted_Scompact_d,
                             const cuFloatComplex *phaseCoeffs_ds,
                             cuFloatComplex *psi_ds,
                             const long *z_ds,
                             const long* y_ds,
                             const size_t numberBeams,
                             const size_t dimk_S,
                             const size_t dimj_S,
                             const size_t dimj_psi,
                             const size_t dimi_psi) {
		// This code is heavily modeled after Mark Harris's presentation on optimized parallel reduction
		// http://developer.download.nvidia.com/compute/cuda/1.1-Beta/x86_website/projects/reduction/doc/reduction.pdf

		// shared memory
		__shared__ cuFloatComplex scaled_values[BlockSizeX]; // for holding the values to reduce
		extern __shared__ cuFloatComplex coeff_cache []; // cache the coefficients to prevent repeated global reads

		// for the permuted Scompact matrix, the x direction runs along the number of beams, leaving y and z to represent the
		//	2D array of reduced values in psi
		int idx = threadIdx.x + blockDim.x * blockIdx.x;
		int y = blockIdx.y;
		int z = blockIdx.z;

		// determine grid size for stepping through the array
		int gridSizeY = gridDim.y * blockDim.y;
		int gridSizeZ = gridDim.z * blockDim.z;

		// guarantee the shared memory is initialized to 0 so we can accumulate without bounds checking
		scaled_values[idx] = make_cuFloatComplex(0,0);
		__syncthreads();

		// read the coefficients into shared memory once
		size_t offset_phase_idx = 0;
		while (offset_phase_idx < numberBeams){
			if (idx + offset_phase_idx < numberBeams){
				coeff_cache[idx + offset_phase_idx] = phaseCoeffs_ds[idx + offset_phase_idx];
			}
		offset_phase_idx += BlockSizeX;
		}
		__syncthreads();

		// each block processes several reductions strided by the grid size
		int y_saved = y;
		while (z < dimj_psi){
			y=y_saved; // reset y
			while(y < dimi_psi){

	//	 read in first values
		if (idx < numberBeams) {
			scaled_values[idx] = cuCmulf(permuted_Scompact_d[z_ds[z]*numberBeams*dimj_S + y_ds[y]*numberBeams + idx],
			                             coeff_cache[idx]);
			__syncthreads();
		}

		// step through global memory accumulating until values have been reduced to BlockSizeX elements in shared memory
		size_t offset = BlockSizeX;
		while (offset < numberBeams){
			if (idx + offset < numberBeams){
				scaled_values[idx] = cuCaddf(scaled_values[idx],
		                                 cuCmulf( permuted_Scompact_d[z_ds[z]*numberBeams*dimj_S + y_ds[y]*numberBeams + idx + offset],
					                              coeff_cache[idx + offset]));
			}
			offset += BlockSizeX;
			__syncthreads();
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

		// use a special optimization for the last reductions
		if (idx < 32)warpReduce_cx<BlockSizeX>(scaled_values, idx);

		// write out the result
		if (idx == 0)psi_ds[z*dimi_psi + y] = scaled_values[0];

		// increment
		y+=gridSizeY;
		__syncthreads();
	}
		z+=gridSizeZ;
		__syncthreads();
	}
}

	template <size_t BlockSizeX>
	// double precision version, see float version above for comments
	__global__ void scaleReduceS(const cuDoubleComplex *permuted_Scompact_d,
	                             const cuDoubleComplex *phaseCoeffs_ds,
	                             cuDoubleComplex *psi_ds,
	                             const long *z_ds,
	                             const long* y_ds,
	                             const size_t numberBeams,
	                             const size_t dimk_S,
	                             const size_t dimj_S,
	                             const size_t dimj_psi,
	                             const size_t dimi_psi) {
		__shared__ cuDoubleComplex scaled_values[BlockSizeX];
		extern __shared__ cuDoubleComplex coeff_cache_double[];
		int idx = threadIdx.x + blockDim.x * blockIdx.x;
		int y = blockIdx.y;
		int z = blockIdx.z;
		int gridSizeY = gridDim.y * blockDim.y;
		int gridSizeZ = gridDim.z * blockDim.z;
		scaled_values[idx] = make_cuDoubleComplex(0,0);
		__syncthreads();
		size_t offset_phase_idx = 0;
		while (offset_phase_idx < numberBeams){
			if (idx + offset_phase_idx < numberBeams){
				coeff_cache_double[idx + offset_phase_idx] = phaseCoeffs_ds[idx + offset_phase_idx];
			}
			offset_phase_idx += BlockSizeX;
		}
		__syncthreads();
		int y_saved = y;
		while (z < dimj_psi){
			y=y_saved;
			while(y < dimi_psi){
				if (idx < numberBeams) {
					scaled_values[idx] = cuCmul(permuted_Scompact_d[z_ds[z]*numberBeams*dimj_S + y_ds[y]*numberBeams + idx],
					                             coeff_cache_double[idx]);
					__syncthreads();
				}
				size_t offset = BlockSizeX;
				while (offset < numberBeams){
					if (idx + offset < numberBeams){
						scaled_values[idx] = cuCadd(scaled_values[idx],
						                             cuCmul( permuted_Scompact_d[z_ds[z]*numberBeams*dimj_S + y_ds[y]*numberBeams + idx + offset],
						                                      coeff_cache_double[idx + offset]));
					}
					offset += BlockSizeX;
					__syncthreads();
				}
				if (BlockSizeX >= 1024){
					if (idx < 512){
						scaled_values[idx] = cuCadd(scaled_values[idx], scaled_values[idx + 512]);
					}
					__syncthreads();
				}
				if (BlockSizeX >= 512){
					if (idx < 256){
						scaled_values[idx] = cuCadd(scaled_values[idx], scaled_values[idx + 256]);
					}
					__syncthreads();
				}
				if (BlockSizeX >= 256){
					if (idx < 128){
						scaled_values[idx] = cuCadd(scaled_values[idx], scaled_values[idx + 128]);
					}
					__syncthreads();
				}
				if (BlockSizeX >= 128){
					if (idx < 64){
						scaled_values[idx] = cuCadd(scaled_values[idx], scaled_values[idx + 64]);
					}
					__syncthreads();
				}
				if (idx < 32)warpReduce_cx<BlockSizeX>(scaled_values, idx);
				if (idx == 0)psi_ds[z*dimi_psi + y] = scaled_values[0];
				y+=gridSizeY;
				__syncthreads();
			}
			z+=gridSizeZ;
			__syncthreads();
		}
	}


	using namespace std;
	void buildPRISMOutput_GPU_singlexfer(Parameters<PRISM_FLOAT_PRECISION> &pars){
		// construct the PRISM output array using GPUs

		// set device flags
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g){
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaSetDeviceFlags(cudaDeviceBlockingSync));
		}

		// create CUDA streams and cuFFT plans
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;
		cudaStream_t *streams   = new cudaStream_t[total_num_streams];
		cufftHandle *cufft_plan = new cufftHandle[total_num_streams];

		for (auto j = 0; j < total_num_streams; ++j) {
			cudaSetDevice(j % pars.meta.NUM_GPUS);
			cudaErrchk(cudaStreamCreate(&streams[j]));
			cufftErrchk(cufftPlan2d(&cufft_plan[j], pars.imageSizeReduce[0], pars.imageSizeReduce[1], PRISM_CUFFT_PLAN_TYPE));
			cufftErrchk(cufftSetStream(cufft_plan[j], streams[j]));
		}


		// pointers to pinned host memory for async transfers
        PRISM_FLOAT_PRECISION               **output_ph = new PRISM_FLOAT_PRECISION*[total_num_streams]; // one output array per stream
//        PRISM_FLOAT_PRECISION               *output_ph[total_num_streams]; // one output array per stream
		std::complex<PRISM_FLOAT_PRECISION> *permuted_Scompact_ph; // see below for explanation of why this is permuted
		std::complex<PRISM_FLOAT_PRECISION> *PsiProbeInit_ph;
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
		cudaErrchk(cudaMallocHost((void **) &PsiProbeInit_ph, pars.psiProbeInit.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
		cudaErrchk(cudaMallocHost((void **) &qxaReduce_ph, pars.qxaReduce.size() * sizeof(PRISM_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **) &qyaReduce_ph, pars.qyaReduce.size() * sizeof(PRISM_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **) &alphaInd_ph,  pars.alphaInd.size()       * sizeof(PRISM_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **) &xBeams_ph, pars.xyBeams.get_dimj()  * sizeof(size_t)));
		cudaErrchk(cudaMallocHost((void **) &yBeams_ph, pars.xyBeams.get_dimj()  * sizeof(size_t)));


		// copy host memory to pinned
		for (auto s = 0; s < total_num_streams; ++s) {
			memset(output_ph[s], 0, pars.stack.get_dimj() * pars.stack.get_dimi() *
			                        sizeof(PRISM_FLOAT_PRECISION));
		}


		// Copy to pinned memory
		// The GPU computational model operates on Scompact in a different order than the CPU, and it is
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

		memcpy(PsiProbeInit_ph, &(*pars.psiProbeInit.begin()), pars.psiProbeInit.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>));
		memcpy(alphaInd_ph,  &(*pars.alphaInd.begin()),       pars.alphaInd.size()       * sizeof(PRISM_FLOAT_PRECISION));
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



		// pointers to read-only GPU memory (one copy per GPU)
		PRISM_CUDA_COMPLEX_FLOAT **permuted_Scompact_d = new PRISM_CUDA_COMPLEX_FLOAT*[pars.meta.NUM_GPUS];
		PRISM_CUDA_COMPLEX_FLOAT **PsiProbeInit_d 		= new PRISM_CUDA_COMPLEX_FLOAT*[pars.meta.NUM_GPUS];
		PRISM_FLOAT_PRECISION    **qxaReduce_d    		= new PRISM_FLOAT_PRECISION*[pars.meta.NUM_GPUS];
		PRISM_FLOAT_PRECISION    **qyaReduce_d    		= new PRISM_FLOAT_PRECISION*[pars.meta.NUM_GPUS];
		PRISM_FLOAT_PRECISION    **alphaInd_d 	  		= new PRISM_FLOAT_PRECISION*[pars.meta.NUM_GPUS];
		size_t                   **yBeams_d 	  		= new size_t*[pars.meta.NUM_GPUS];
		size_t                   **xBeams_d 	 		 = new size_t*[pars.meta.NUM_GPUS];

		// pointers to read/write GPU memory (one per stream)

		PRISM_CUDA_COMPLEX_FLOAT **psi_ds 				= new PRISM_CUDA_COMPLEX_FLOAT*[total_num_streams];
		PRISM_CUDA_COMPLEX_FLOAT **phaseCoeffs_ds 		= new PRISM_CUDA_COMPLEX_FLOAT*[total_num_streams];
		PRISM_FLOAT_PRECISION    **psi_intensity_ds 	= new PRISM_FLOAT_PRECISION*[total_num_streams];
		PRISM_FLOAT_PRECISION    **integratedOutput_ds  = new PRISM_FLOAT_PRECISION*[total_num_streams];
		long                     **y_ds = new long*[total_num_streams];
		long                     **x_ds = new long*[total_num_streams];


		// allocate memory on each GPU
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaMalloc((void **) &permuted_Scompact_d[g], pars.Scompact.size() * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMalloc((void **) &PsiProbeInit_d[g],      pars.psiProbeInit.size()  * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMalloc((void **) &qxaReduce_d[g], pars.qxaReduce.size() * sizeof(PRISM_FLOAT_PRECISION)));
			cudaErrchk(cudaMalloc((void **) &qyaReduce_d[g], pars.qyaReduce.size() * sizeof(PRISM_FLOAT_PRECISION)));
			cudaErrchk(cudaMalloc((void **) &alphaInd_d[g],  pars.alphaInd.size()       * sizeof(PRISM_FLOAT_PRECISION)));
			cudaErrchk(cudaMalloc((void **) &yBeams_d[g], pars.xyBeams.get_dimj()  * sizeof(size_t)));
			cudaErrchk(cudaMalloc((void **) &xBeams_d[g], pars.xyBeams.get_dimj()  * sizeof(size_t)));
		}

		// allocate memory per stream and zero it
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
		// have more than one stream per GPU, then we try to interleave as much as possible
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
			                           pars.psiProbeInit.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>),
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
			                           pars.alphaInd.size() * sizeof(pars.alphaInd[0]), cudaMemcpyHostToDevice, streams[stream_id]));
			cout << "stream_id = " << stream_id << endl;
			cudaErrchk(cudaMemcpyAsync(yBeams_d[g], &yBeams_ph[0],
			                           pars.xyBeams.get_dimj() * sizeof(size_t), cudaMemcpyHostToDevice,
			                           streams[stream_id]));
			cudaErrchk(cudaMemcpyAsync(xBeams_d[g], &xBeams_ph[0],
			                           pars.xyBeams.get_dimj() * sizeof(size_t), cudaMemcpyHostToDevice,
			                           streams[stream_id]));
		}

		// wait for transfers to complete
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaDeviceSynchronize());
		}


		// launch threads that will consume work provided by getWorkID
		vector<thread> workers_GPU;
		workers_GPU.reserve(total_num_streams); // prevents multiple reallocations
		int stream_count = 0;
//		setWorkStartStop(0, pars.xp.size() * pars.yp.size(), 1);
		WorkDispatcher dispatcher(0, pars.xp.size() * pars.yp.size(), 1);

//		setWorkStartStop(0, 1, 1);
		for (auto t = 0; t < total_num_streams; ++t) {

			int GPU_num = stream_count % pars.meta.NUM_GPUS; // determine which GPU handles this job

			cudaStream_t &current_stream = streams[stream_count];
			cout << "Launching GPU worker on stream #" << stream_count << " of GPU #" << GPU_num << '\n';

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

			// push_back is better whenever constructing a new object
			workers_GPU.push_back(thread([&pars, GPU_num, stream_count, current_permuted_Scompact_d, &dispatcher,													 current_alphaInd_d, current_PsiProbeInit_d, current_qxaReduce_d, current_qyaReduce_d,
					                                current_yBeams_d, current_xBeams_d, current_psi_ds, current_phaseCoeffs_ds,
					                                current_psi_intensity_ds, current_y_ds, current_x_ds, current_integratedOutput_ds,
					                                current_output_ph, &current_cufft_plan, &current_stream]() {
				cudaErrchk(cudaSetDevice(GPU_num));
				size_t Nstart, Nstop, ay, ax;
                Nstart=Nstop=0;
//				while (getWorkID(pars, Nstart, Nstop)) { // synchronously get work assignment
				while (dispatcher.getWork(Nstart, Nstop)) { // synchronously get work assignment
					while (Nstart != Nstop) {
						ay = Nstart / pars.xp.size();
						ax = Nstart % pars.xp.size();
						buildSignal_GPU_singlexfer(pars, ay, ax, current_permuted_Scompact_d,
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


		// Now launch CPU work
		if (pars.meta.also_do_CPU_work) {
			PRISM_FFTW_INIT_THREADS();
			PRISM_FFTW_PLAN_WITH_NTHREADS(pars.meta.NUM_THREADS);
			vector<thread> workers_CPU;
			workers_CPU.reserve(pars.meta.NUM_THREADS); // prevents multiple reallocations
//			for (auto t = 0; t < pars.meta.NUM_THREADS; ++t) {
			for (auto t = 0; t < 1; ++t) {
				cout << "Launching CPU worker thread #" << t << " to compute partial PRISM result\n";
				// push_back is better whenever constructing a new object
				workers_CPU.push_back(thread([&pars, &dispatcher, t]() {
					size_t Nstart, Nstop, ay, ax, early_CPU_stop;
                    Nstart=Nstop=0;
//					early_CPU_stop = pars.xp.size() * pars.yp.size();
					early_CPU_stop = pars.xp.size() * pars.yp.size() - (1./pars.meta.cpu_gpu_ratio);
//					while (getWorkID(pars, Nstart, Nstop)) { // synchronously get work assignment
					if(dispatcher.getWork(Nstart, Nstop)) { // synchronously get work assignment
						Array2D<std::complex<PRISM_FLOAT_PRECISION> > psi = PRISM::zeros_ND<2, std::complex<PRISM_FLOAT_PRECISION> > (
								{{pars.imageSizeReduce[0], pars.imageSizeReduce[1]}});
						unique_lock<mutex> gatekeeper(fftw_plan_lock);


						PRISM_FFTW_PLAN plan = PRISM_FFTW_PLAN_DFT_2D(psi.get_dimj(), psi.get_dimi(),
																	  reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi[0]),
																	  reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi[0]),
																	  FFTW_FORWARD, FFTW_MEASURE);
						gatekeeper.unlock();
						do {
							while (Nstart != Nstop) {
								ay = Nstart / pars.xp.size();
								ax = Nstart % pars.xp.size();
								buildSignal_CPU(pars, ay, ax, plan, psi);
								++Nstart;
							}
if (Nstop >= early_CPU_stop) break;
						} while(dispatcher.getWork(Nstart, Nstop));
						gatekeeper.lock();
						PRISM_FFTW_DESTROY_PLAN(plan);
						gatekeeper.unlock();
					}
				}));
			}
			cout << "Waiting for CPU threads...\n";
			for (auto& t:workers_CPU)t.join();
			PRISM_FFTW_CLEANUP_THREADS();
		}


		// synchronize
		cout << "Waiting for GPU threads...\n";
		for (auto &t:workers_GPU)t.join();

		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g){
			cudaSetDevice(g);
			cudaDeviceSynchronize();
		}

		// free pinned memory
		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaFreeHost(output_ph[s]));
		}
		cudaErrchk(cudaFreeHost(permuted_Scompact_ph));
		cudaErrchk(cudaFreeHost(PsiProbeInit_ph));
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

		for (auto j = 0; j < pars.meta.NUM_GPUS; ++j) {
			cudaErrchk(cudaSetDevice(j));
			cudaErrchk(cudaDeviceReset());
		}

		delete[] streams;
		delete[] cufft_plan;
		delete[] output_ph;
		delete[] PsiProbeInit_d;
		delete[] qxaReduce_d;
		delete[] qyaReduce_d;
		delete[] alphaInd_d;
		delete[] yBeams_d;
		delete[] xBeams_d;
		delete[] permuted_Scompact_d;
		delete[] psi_ds;
		delete[] phaseCoeffs_ds;
		delete[] psi_intensity_ds;
		delete[] integratedOutput_ds;
		delete[] y_ds;
		delete[] x_ds;
	}

	void buildPRISMOutput_GPU_streaming(Parameters<PRISM_FLOAT_PRECISION> &pars){
		// construct the PRISM output array using GPUs

		// set device flags
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g){
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaSetDeviceFlags(cudaDeviceBlockingSync));
		}

		// create CUDA streams and cuFFT plans
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;
        cudaStream_t *streams   = new cudaStream_t[total_num_streams];
        cufftHandle *cufft_plan = new cufftHandle[total_num_streams];
//		cudaStream_t streams[total_num_streams];
//		cufftHandle cufft_plan[total_num_streams];

		for (auto j = 0; j < total_num_streams; ++j) {
			cudaSetDevice(j % pars.meta.NUM_GPUS);
			cudaErrchk(cudaStreamCreate(&streams[j]));
			cufftErrchk(cufftPlan2d(&cufft_plan[j], pars.imageSizeReduce[0], pars.imageSizeReduce[1], PRISM_CUFFT_PLAN_TYPE));
			cufftErrchk(cufftSetStream(cufft_plan[j], streams[j]));
		}


		// pointers to pinned host memory for async transfers
//		PRISM_FLOAT_PRECISION               *output_ph[total_num_streams]; // one output array per stream
		PRISM_FLOAT_PRECISION               **output_ph = new PRISM_FLOAT_PRECISION*[total_num_streams];
		std::complex<PRISM_FLOAT_PRECISION> *permuted_Scompact_ph; // see below for explanation of why this is permuted
		std::complex<PRISM_FLOAT_PRECISION> *PsiProbeInit_ph;
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
		cudaErrchk(cudaMallocHost((void **) &PsiProbeInit_ph, pars.psiProbeInit.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
		cudaErrchk(cudaMallocHost((void **) &qxaReduce_ph, pars.qxaReduce.size() * sizeof(PRISM_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **) &qyaReduce_ph, pars.qyaReduce.size() * sizeof(PRISM_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **) &alphaInd_ph,  pars.alphaInd.size()       * sizeof(PRISM_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **) &xBeams_ph, pars.xyBeams.get_dimj()  * sizeof(size_t)));
		cudaErrchk(cudaMallocHost((void **) &yBeams_ph, pars.xyBeams.get_dimj()  * sizeof(size_t)));


		// copy host memory to pinned
		for (auto s = 0; s < total_num_streams; ++s) {
			memset(output_ph[s], 0, pars.stack.get_dimj() * pars.stack.get_dimi() *
			                        sizeof(PRISM_FLOAT_PRECISION));
		}


		// Copy to pinned memory
		// The GPU computational model operates on Scompact in a different order than the CPU, and it is
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

		memcpy(PsiProbeInit_ph, &(*pars.psiProbeInit.begin()), pars.psiProbeInit.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>));
		memcpy(alphaInd_ph,  &(*pars.alphaInd.begin()),       pars.alphaInd.size()       * sizeof(PRISM_FLOAT_PRECISION));
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



		// pointers to read-only GPU memory (one copy per GPU)
		PRISM_CUDA_COMPLEX_FLOAT **PsiProbeInit_d = new PRISM_CUDA_COMPLEX_FLOAT*[pars.meta.NUM_GPUS];
		PRISM_FLOAT_PRECISION    **qxaReduce_d    = new PRISM_FLOAT_PRECISION*[pars.meta.NUM_GPUS];
		PRISM_FLOAT_PRECISION    **qyaReduce_d    = new PRISM_FLOAT_PRECISION*[pars.meta.NUM_GPUS];
		PRISM_FLOAT_PRECISION    **alphaInd_d 	  = new PRISM_FLOAT_PRECISION*[pars.meta.NUM_GPUS];
		size_t                   **yBeams_d 	  = new size_t*[pars.meta.NUM_GPUS];
		size_t                   **xBeams_d 	  = new size_t*[pars.meta.NUM_GPUS];
//		PRISM_CUDA_COMPLEX_FLOAT *PsiProbeInit_d[pars.meta.NUM_GPUS];
//		PRISM_FLOAT_PRECISION    *qxaReduce_d[pars.meta.NUM_GPUS];
//		PRISM_FLOAT_PRECISION    *qyaReduce_d[pars.meta.NUM_GPUS];
//		PRISM_FLOAT_PRECISION    *alphaInd_d[pars.meta.NUM_GPUS];
//		size_t                   *yBeams_d[pars.meta.NUM_GPUS];
//		size_t                   *xBeams_d[pars.meta.NUM_GPUS];

		// pointers to read/write GPU memory (one per stream)
		PRISM_CUDA_COMPLEX_FLOAT **permuted_Scompact_ds = new PRISM_CUDA_COMPLEX_FLOAT*[total_num_streams];
		PRISM_CUDA_COMPLEX_FLOAT **psi_ds 				= new PRISM_CUDA_COMPLEX_FLOAT*[total_num_streams];
		PRISM_CUDA_COMPLEX_FLOAT **phaseCoeffs_ds 		= new PRISM_CUDA_COMPLEX_FLOAT*[total_num_streams];
		PRISM_FLOAT_PRECISION    **psi_intensity_ds 	= new PRISM_FLOAT_PRECISION*[total_num_streams];
		PRISM_FLOAT_PRECISION    **integratedOutput_ds  = new PRISM_FLOAT_PRECISION*[total_num_streams];
		long                     **y_ds = new long*[total_num_streams];
		long                     **x_ds = new long*[total_num_streams];
//		PRISM_CUDA_COMPLEX_FLOAT *permuted_Scompact_ds[total_num_streams];
//		PRISM_CUDA_COMPLEX_FLOAT *psi_ds[total_num_streams];
//		PRISM_CUDA_COMPLEX_FLOAT *phaseCoeffs_ds[total_num_streams];
//		PRISM_FLOAT_PRECISION    *psi_intensity_ds[total_num_streams];
//		PRISM_FLOAT_PRECISION    *integratedOutput_ds[total_num_streams];
//		long                     *y_ds[total_num_streams];
//		long                     *x_ds[total_num_streams];

		// allocate memory on each GPU
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaMalloc((void **) &PsiProbeInit_d[g],      pars.psiProbeInit.size()  * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMalloc((void **) &qxaReduce_d[g], pars.qxaReduce.size() * sizeof(PRISM_FLOAT_PRECISION)));
			cudaErrchk(cudaMalloc((void **) &qyaReduce_d[g], pars.qyaReduce.size() * sizeof(PRISM_FLOAT_PRECISION)));
			cudaErrchk(cudaMalloc((void **) &alphaInd_d[g],  pars.alphaInd.size()       * sizeof(PRISM_FLOAT_PRECISION)));
			cudaErrchk(cudaMalloc((void **) &yBeams_d[g], pars.xyBeams.get_dimj()  * sizeof(size_t)));
			cudaErrchk(cudaMalloc((void **) &xBeams_d[g], pars.xyBeams.get_dimj()  * sizeof(size_t)));
		}

		// allocate memory per stream and zero it
		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaMalloc((void **) &permuted_Scompact_ds[s],
			                      pars.imageSizeReduce[0] * pars.imageSizeReduce[1] * pars.numberBeams * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
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
		// have more than one stream per GPU, then we try to interleave as much as possible
		int stream_id = 0;
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			stream_id = g;
			cudaErrchk(cudaSetDevice(g));
//			cudaErrchk(cudaMemcpyAsync(permuted_Scompact_d[g], &permuted_Scompact_ph[0],
//			                           pars.Scompact.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>),
//			                           cudaMemcpyHostToDevice, streams[stream_id]));

			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cout << "stream_id = " << stream_id << endl;
			cudaErrchk(cudaMemcpyAsync(PsiProbeInit_d[g], &PsiProbeInit_ph[0],
			                           pars.psiProbeInit.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>),
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
			                           pars.alphaInd.size() * sizeof(pars.alphaInd[0]), cudaMemcpyHostToDevice, streams[stream_id]));
			cout << "stream_id = " << stream_id << endl;
			cudaErrchk(cudaMemcpyAsync(yBeams_d[g], &yBeams_ph[0],
			                           pars.xyBeams.get_dimj() * sizeof(size_t), cudaMemcpyHostToDevice,
			                           streams[stream_id]));
			cudaErrchk(cudaMemcpyAsync(xBeams_d[g], &xBeams_ph[0],
			                           pars.xyBeams.get_dimj() * sizeof(size_t), cudaMemcpyHostToDevice,
			                           streams[stream_id]));
		}

		// wait for transfers to complete
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaDeviceSynchronize());
		}


		// launch threads that will consume work provided by getWorkID
		vector<thread> workers_GPU;
		workers_GPU.reserve(total_num_streams); // prevents multiple reallocations
		int stream_count = 0;
//		setWorkStartStop(0, pars.xp.size() * pars.yp.size());
		WorkDispatcher dispatcher(0, pars.xp.size() * pars.yp.size(), 1);

//		setWorkStartStop(0, 1);
		for (auto t = 0; t < total_num_streams; ++t) {

			int GPU_num = stream_count % pars.meta.NUM_GPUS; // determine which GPU handles this job

			cudaStream_t &current_stream = streams[stream_count];
			cout << "Launching GPU worker on stream #" << stream_count << " of GPU #" << GPU_num << '\n';

			// get pointers to the pre-copied arrays, making sure to get those on the current GPU
			PRISM_CUDA_COMPLEX_FLOAT *current_PsiProbeInit_d      = PsiProbeInit_d[GPU_num];
			PRISM_FLOAT_PRECISION *current_qxaReduce_d            = qxaReduce_d[GPU_num];
			PRISM_FLOAT_PRECISION *current_qyaReduce_d            = qyaReduce_d[GPU_num];
			size_t *current_yBeams_d                              = yBeams_d[GPU_num];
			size_t *current_xBeams_d                              = xBeams_d[GPU_num];
			PRISM_FLOAT_PRECISION *current_alphaInd_d             = alphaInd_d[GPU_num];

			// get pointers to per-stream arrays
			PRISM_CUDA_COMPLEX_FLOAT *current_permuted_Scompact_ds= permuted_Scompact_ds[stream_count];
			PRISM_CUDA_COMPLEX_FLOAT *current_psi_ds              = psi_ds[stream_count];
			PRISM_CUDA_COMPLEX_FLOAT *current_phaseCoeffs_ds      = phaseCoeffs_ds[stream_count];
			PRISM_FLOAT_PRECISION *current_psi_intensity_ds       = psi_intensity_ds[stream_count];
			PRISM_FLOAT_PRECISION *current_integratedOutput_ds    = integratedOutput_ds[stream_count];
			long *current_y_ds                                    = y_ds[stream_count];
			long *current_x_ds                                    = x_ds[stream_count];
			cufftHandle &current_cufft_plan                       = cufft_plan[stream_count];

			// get pointer to output pinned memory
			PRISM_FLOAT_PRECISION *current_output_ph              = output_ph[stream_count];

			// push_back is better whenever constructing a new object
			workers_GPU.push_back(thread([&pars, GPU_num, stream_count, current_permuted_Scompact_ds, permuted_Scompact_ph, &dispatcher,
					                                current_alphaInd_d, current_PsiProbeInit_d, current_qxaReduce_d, current_qyaReduce_d,
					                                current_yBeams_d, current_xBeams_d, current_psi_ds, current_phaseCoeffs_ds,
					                                current_psi_intensity_ds, current_y_ds, current_x_ds, current_integratedOutput_ds,
					                                current_output_ph, &current_cufft_plan, &current_stream]() {
				cudaErrchk(cudaSetDevice(GPU_num));
				size_t Nstart, Nstop, ay, ax;
                Nstart=Nstop=0;
//				while (getWorkID(pars, Nstart, Nstop)) { // synchronously get work assignment
				while (dispatcher.getWork(Nstart, Nstop)) { // synchronously get work assignment
					while (Nstart != Nstop) {
						ay = Nstart / pars.xp.size();
						ax = Nstart % pars.xp.size();
						buildSignal_GPU_streaming(pars, ay, ax, current_permuted_Scompact_ds, permuted_Scompact_ph,
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


		// Now launch CPU work
		if (pars.meta.also_do_CPU_work) {
			PRISM_FFTW_INIT_THREADS();
			PRISM_FFTW_PLAN_WITH_NTHREADS(pars.meta.NUM_THREADS);
			vector<thread> workers_CPU;
			workers_CPU.reserve(pars.meta.NUM_THREADS); // prevents multiple reallocations
			for (auto t = 0; t < pars.meta.NUM_THREADS; ++t) {
				cout << "Launching CPU worker thread #" << t << " to compute partial PRISM result\n";
				// push_back is better whenever constructing a new object
				workers_CPU.push_back(thread([&pars, &dispatcher, t]() {
					size_t Nstart, Nstop, ay, ax, early_CPU_stop;
                    Nstart=Nstop=0;
//					early_CPU_stop = pars.xp.size() * pars.yp.size() * (1-pars.meta.cpu_gpu_ratio);
					early_CPU_stop = pars.xp.size() * pars.yp.size() - (1./pars.meta.cpu_gpu_ratio);
//					early_CPU_stop = 0;
//					while (getWorkID(pars, Nstart, Nstop)) { // synchronously get work assignment
					if(dispatcher.getWork(Nstart, Nstop)) { // synchronously get work assignment
						Array2D<std::complex<PRISM_FLOAT_PRECISION> > psi = PRISM::zeros_ND<2, std::complex<PRISM_FLOAT_PRECISION> > (
								{{pars.imageSizeReduce[0], pars.imageSizeReduce[1]}});
						unique_lock<mutex> gatekeeper(fftw_plan_lock);


						PRISM_FFTW_PLAN plan = PRISM_FFTW_PLAN_DFT_2D(psi.get_dimj(), psi.get_dimi(),
																	  reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi[0]),
																	  reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi[0]),
																	  FFTW_FORWARD, FFTW_MEASURE);
						gatekeeper.unlock();
						do {
							while (Nstart != Nstop) {
								ay = Nstart / pars.xp.size();
								ax = Nstart % pars.xp.size();
								buildSignal_CPU(pars, ay, ax, plan, psi);
								++Nstart;
							}
if (Nstop >= early_CPU_stop) break;
						} while(dispatcher.getWork(Nstart, Nstop));
						gatekeeper.lock();
						PRISM_FFTW_DESTROY_PLAN(plan);
						gatekeeper.unlock();
					}
				}));
			}
			cout << "Waiting for CPU threads...\n";
			for (auto& t:workers_CPU)t.join();
			PRISM_FFTW_CLEANUP_THREADS();
		}


		// synchronize
		cout << "Waiting for GPU threads...\n";
		for (auto &t:workers_GPU)t.join();

		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g){
			cudaSetDevice(g);
			cudaDeviceSynchronize();
		}

		// free pinned memory
		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaFreeHost(output_ph[s]));
		}
		cudaErrchk(cudaFreeHost(permuted_Scompact_ph));
		cudaErrchk(cudaFreeHost(PsiProbeInit_ph));
		cudaErrchk(cudaFreeHost(qxaReduce_ph));
		cudaErrchk(cudaFreeHost(qyaReduce_ph));
		cudaErrchk(cudaFreeHost(xBeams_ph));
		cudaErrchk(cudaFreeHost(yBeams_ph));

		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaFree(permuted_Scompact_ds[g]));
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

		for (auto j = 0; j < pars.meta.NUM_GPUS; ++j) {
			cudaErrchk(cudaSetDevice(j));
			cudaErrchk(cudaDeviceReset());
		}
        delete[] streams;
        delete[] cufft_plan;
		delete[] output_ph;
		delete[] PsiProbeInit_d;
		delete[] qxaReduce_d;
		delete[] qyaReduce_d;
		delete[] alphaInd_d;
		delete[] yBeams_d;
		delete[] xBeams_d;
		delete[] permuted_Scompact_ds;
		delete[] psi_ds;
		delete[] phaseCoeffs_ds;
		delete[] psi_intensity_ds;
		delete[] integratedOutput_ds;
		delete[] y_ds;
		delete[] x_ds;
	}


	void buildSignal_GPU_singlexfer(Parameters<PRISM_FLOAT_PRECISION>&  pars,
	                                const size_t& ay,
	                                const size_t& ax,
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

		const PRISM_FLOAT_PRECISION yp = pars.yp[ay];
		const PRISM_FLOAT_PRECISION xp = pars.xp[ax];
		const size_t psi_size = pars.imageSizeReduce[0] * pars.imageSizeReduce[1];
		shiftIndices <<<(pars.imageSizeReduce[0] - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream>>> (
			y_ds, std::round(yp / pars.pixelSizeOutput[0]),pars.imageSize[0], pars.imageSizeReduce[0]);

		shiftIndices <<<(pars.imageSizeReduce[1] - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream>>> (
			x_ds, std::round(xp / pars.pixelSizeOutput[1]), pars.imageSize[1], pars.imageSizeReduce[1]);

		computePhaseCoeffs <<<(pars.numberBeams - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream>>>(
                   phaseCoeffs_ds, PsiProbeInit_d, qyaReduce_d, qxaReduce_d,
		           yBeams_d, xBeams_d, yp, xp, pars.yTiltShift, pars.xTiltShift, pars.imageSizeReduce[1], pars.numberBeams);


//if (ay==0&&ax==0)		cout << "pars.imageSizeReduce[0] = " << pars.imageSizeReduce[0] << endl;
//if (ay==0&&ax==0)		cout << "pars.imageSizeReduce[1] = " << pars.imageSizeReduce[1] << endl;


		// Choose a good launch configuration
		// Heuristically use 2^p / 2 as the block size where p is the first power of 2 greater than the number of elements to work on.
		// This balances having enough work per thread and enough blocks without having so many blocks that the shared memory doesn't last long
		size_t p = getNextPower2(pars.numberBeams);
		const size_t BlockSizeX = (size_t)std::max(1.0, pow(2,p) / 2);

		// Determine maximum threads per streaming multiprocessor based on the compute capability of the device
		size_t max_threads_per_sm;
		if (pars.deviceProperties.major > 3){
			max_threads_per_sm = 2048;
		} else if (pars.deviceProperties.major > 2) {
			max_threads_per_sm = 1536;
		} else {
			max_threads_per_sm = pars.deviceProperties.minor == 0 ? 768 : 1024;
		}

		// Estimate max number of blocks per streaming multiprocessor (threads are the limit)
		const size_t max_blocks_per_sm = max_threads_per_sm / BlockSizeX;

		// We find providing around 3 times as many blocks as the estimated maximum provides good performance
		const size_t target_blocks_per_sm = max_blocks_per_sm * 3;
		const size_t total_blocks         = target_blocks_per_sm * pars.deviceProperties.multiProcessorCount;

		// Determine the shape of the grid
		const PRISM_FLOAT_PRECISION aspect_ratio = pars.imageSizeReduce[1] / pars.imageSizeReduce[0];
		const size_t BlockSizeZ = std::floor(sqrt(total_blocks / aspect_ratio));
		const size_t BlockSizeY = aspect_ratio * BlockSizeZ;

//		if (ax == 0 & ay == 0) {
//			cout << "aspect_ratio = " << aspect_ratio << endl;
//			cout << "BlockSizeX = " << BlockSizeX << endl;
//			cout << "BlockSizeZ = " << BlockSizeZ << endl;
//			cout << "BlockSizeY = " << BlockSizeY << endl;
//			cout << "target_blocks_per_sm = " << target_blocks_per_sm << endl;
//			cout << "total_blocks = " << total_blocks << endl;
//			cout << " pars.deviceProperties.multiProcessorCount = " <<  pars.deviceProperties.multiProcessorCount << endl;
//		}
		dim3 grid(1, BlockSizeY, BlockSizeZ);
		dim3 block(BlockSizeX, 1, 1);
//		dim3 block(4, 1, 1);

		// Determine amount of shared memory needed
		const unsigned long smem = pars.numberBeams * sizeof(PRISM_CUDA_COMPLEX_FLOAT);

//		scaleReduceS<4> <<< grid, block, smem, stream >>> (
//				permuted_Scompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
//				pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]);

		// Launch kernel. Block size must be visible at compile time so we use a switch statement
		switch (BlockSizeX) {
			case 1024 :
			scaleReduceS<1024> <<< grid, block, smem, stream >>> (
					permuted_Scompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
					pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
			case 512 :
			scaleReduceS<512> <<< grid, block, smem, stream >>> (
					permuted_Scompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
					pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
			case 256 :
			scaleReduceS<256> <<< grid, block, smem, stream >>> (
					permuted_Scompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
					pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
			case 128 :
			scaleReduceS<128> <<< grid, block, smem, stream >>> (
					permuted_Scompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
					pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
			case 64 :
			scaleReduceS<64> <<< grid, block, smem, stream >>> (
					permuted_Scompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
					pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
			case 32 :
			scaleReduceS<32> <<< grid, block, smem, stream >>> (
					permuted_Scompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
					pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
			case 16 :
			scaleReduceS<16> <<< grid, block, smem, stream >>> (
					permuted_Scompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
					pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
			case 8 :
			scaleReduceS<8> <<< grid, block, smem, stream >>> (
					permuted_Scompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
					pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
			case 4 :
			scaleReduceS<4> <<< grid, block, smem, stream >>> (
					permuted_Scompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
					pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
			case 2 :
			scaleReduceS<2> <<< grid, block, smem, stream >>> (
					permuted_Scompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
					pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
			default :
			scaleReduceS<1> <<< grid, block, smem, stream >>> (
					permuted_Scompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
					pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
		}

		/*
		{
////			cudaDeviceSynchronize();
			complex<PRISM_FLOAT_PRECISION > ans;
			for (auto i = 0; i < pars.imageSizeReduce[0]; ++i){
				cudaErrchk(cudaMemcpy(&ans, (psi_ds + i), sizeof(complex<PRISM_FLOAT_PRECISION>),cudaMemcpyDeviceToHost));
				cout << "psi_ds[" << i << "] = " << ans << endl;
			}
		}
		*/

		// final fft
		cufftErrchk(PRISM_CUFFT_EXECUTE(cufft_plan, &psi_ds[0], &psi_ds[0], CUFFT_FORWARD));

		// convert to squared intensity
		abs_squared <<< (psi_size - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >>> (psi_intensity_ds, psi_ds, psi_size);

		// output calculation result
		formatOutput_GPU_integrate(pars, psi_intensity_ds, alphaInd_d, output_ph,
		                           integratedOutput_ds, ay, ax, pars.imageSizeReduce[0],
		                           pars.imageSizeReduce[1], stream, pars.scale);

	}


















	void buildSignal_GPU_streaming(Parameters<PRISM_FLOAT_PRECISION>&  pars,
	                               const size_t& ay,
	                               const size_t& ax,
	                               PRISM_CUDA_COMPLEX_FLOAT *permuted_Scompact_ds,
	                               const std::complex<PRISM_FLOAT_PRECISION> *permuted_Scompact_ph,
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
		// the coordinates y and x of the output image phi map to z and y of the permuted S compact matrix
		const PRISM_FLOAT_PRECISION yp = pars.yp[ay];
		const PRISM_FLOAT_PRECISION xp = pars.xp[ax];
		const size_t psi_size = pars.imageSizeReduce[0] * pars.imageSizeReduce[1];
		shiftIndices <<<(pars.imageSizeReduce[0] - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream>>> (
				y_ds, std::round(yp / pars.pixelSizeOutput[0]),pars.imageSize[0], pars.imageSizeReduce[0]);

		shiftIndices <<<(pars.imageSizeReduce[1] - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream>>> (
				x_ds, std::round(xp / pars.pixelSizeOutput[1]), pars.imageSize[1], pars.imageSizeReduce[1]);

		computePhaseCoeffs <<<(pars.numberBeams - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream>>>(
				phaseCoeffs_ds, PsiProbeInit_d, qyaReduce_d, qxaReduce_d,
				yBeams_d, xBeams_d, yp, xp, pars.yTiltShift, pars.xTiltShift, pars.imageSizeReduce[1], pars.numberBeams);

		long x1,y1;
		y1 = pars.yVec[0] + std::round(yp / pars.pixelSizeOutput[0]);
		x1 = pars.xVec[0] + std::round(xp / pars.pixelSizeOutput[1]);
//		if (ay == 0 && ax == 0) {
//			cout << "x1 = " << x1 << endl;
//			cout << "y1 = " << y1 << endl;
//			cout << " pars.xVec[0 = " <<  pars.xVec[0] << endl;
//			cout << " pars.yVec[0 = " <<  pars.yVec[0] << endl;
//
//		}
		if ( (y1 + pars.imageSizeReduce[0]) > pars.Scompact.get_dimj()){throw std::domain_error("Attempting to stream out of bounds value due to wrap around");}
		if ( (x1 + pars.imageSizeReduce[1]) > pars.Scompact.get_dimi()){throw std::domain_error("Attempting to stream out of bounds value due to wrap around");}

//		for (auto row = 0; row < pars.imageSizeReduce[0]; ++row){
//			cudaErrchk(cudaMemcpy(&permuted_Scompact_ds[row*pars.imageSizeReduce[1]*pars.numberBeams],
//			                      &permuted_Scompact_ph[(y1 + row)*pars.numberBeams*pars.Scompact.get_dimi() + x1 * pars.numberBeams],
//			                      pars.imageSizeReduce[1] * pars.numberBeams * sizeof(PRISM_CUDA_COMPLEX_FLOAT),
//			                      cudaMemcpyHostToDevice));
//		}
		cudaErrchk(cudaMemcpy2DAsync(permuted_Scompact_ds,
		                             pars.imageSizeReduce[1] * pars.numberBeams * sizeof(PRISM_CUDA_COMPLEX_FLOAT),
		                             &permuted_Scompact_ph[y1*pars.numberBeams*pars.Scompact.get_dimi() + x1 * pars.numberBeams],
		                             pars.Scompact.get_dimi() * pars.numberBeams * sizeof(PRISM_CUDA_COMPLEX_FLOAT), // corresponds to stride between permuted Scompact elements in k-direction
		                             pars.imageSizeReduce[1] * pars.numberBeams  * sizeof(PRISM_CUDA_COMPLEX_FLOAT),
		                             pars.imageSizeReduce[0],
		                             cudaMemcpyHostToDevice,
		                             stream));



//		cudaErrchk(cudaMemcpy2DAsync(permuted_Scompact_ds,
//		                             pars.imageSizeReduce[1] * sizeof(PRISM_CUDA_COMPLEX_FLOAT),
//		                             &permuted_Scompact_ph[0],
//		                             pars.Scompact.get_dimi() * pars.Scompact.get_dimk() * sizeof(PRISM_CUDA_COMPLEX_FLOAT), // corresponds to stride between permuted Scompact elements in k-direction
//		                             pars.imageSizeReduce[1] * sizeof(PRISM_CUDA_COMPLEX_FLOAT),
//		                             pars.imageSizeReduce[0],
//		                             cudaMemcpyHostToDevice,
//		                             stream));
//		cudaErrchk(cudaMemcpy(permuted_Scompact_ds,
//		                      &permuted_Scompact_ph[0],
//		                             pars.imageSizeReduce[1] *pars.imageSizeReduce[0]*pars.numberBeams* sizeof(PRISM_CUDA_COMPLEX_FLOAT),
//		                             cudaMemcpyHostToDevice));
//		Array3D<complex<float> > debug_sc = zeros_ND<3, complex<float> >({{pars.imageSizeReduce[0], pars.imageSizeReduce[1], pars.numberBeams}});
//	if (ax == 0 && ay == 0) {
//		cout << "pars.Scompact.at(0,0,0) = " << pars.Scompact.at(0,0,0) << endl;
//		cout << "pars.Scompact.at(1,0,0) = " << pars.Scompact.at(1,0,0) << endl;
//		cout << "pars.Scompact.at(2,0,0) = " << pars.Scompact.at(2,0,0) << endl;
//
//		cudaErrchk(cudaMemcpy(&debug_sc[0],
//		                      permuted_Scompact_ds,
//		                      pars.imageSizeReduce[0] * pars.imageSizeReduce[1] * pars.numberBeams *
//		                      sizeof(PRISM_CUDA_COMPLEX_FLOAT),
//		                      cudaMemcpyDeviceToHost));
//
//		for (auto i = 0; i < 10; ++i) {
//			cout << "debug_sc[" << i << "] = " << debug_sc[i] << endl;
//		}
//
//		cout << "pars.Scompact.at(0,y1,x1) = " << pars.Scompact.at(0,y1,x1) << endl;
//		cout << "pars.Scompact.at(1,y1,x1) = " << pars.Scompact.at(1,y1,x1) << endl;
//		cout << "pars.Scompact.at(2,y1,x1) = " << pars.Scompact.at(2,y1,x1) << endl;
//		cout << "pars.Scompact.at(0,y1,x1 + 1) = " << pars.Scompact.at(0,y1,x1 + 1) << endl;
//		cout << "pars.Scompact.at(0,y1+1,x1) = " << pars.Scompact.at(0,y1+1,x1) << endl;
//		cout << "pars.Scompact.at(1,y1+1,x1) = " << pars.Scompact.at(1,y1+1,x1) << endl;
//		cout << "pars.Scompact.at(2,y1+1,x1+1) = " << pars.Scompact.at(2,y1+1,x1+1) << endl;
//		cout << "pars.Scompact.at(0,0,1) = " << pars.Scompact.at(0,0,1) << endl;
//		cout << "pars.Scompact.at(0,1,0) = " << pars.Scompact.at(0,1,0) << endl;
//
//
//
//		cout << "debug_sc.at(0,0,0) = " << debug_sc.at(0,0,0) << endl;
//		cout << "debug_sc.at(1,0,0) = " << debug_sc.at(0,0,1) << endl;
//		cout << "debug_sc.at(2,0,0) = " << debug_sc.at(0,0,2) << endl;
//		cout << "debug_sc.at(0,0+1,0) = " << debug_sc.at(1,0,0) << endl;
//		cout << "debug_sc.at(1,0+1,0) = " << debug_sc.at(1,0,1) << endl;
//		cout << "debug_sc.at(2,0+1,0+1) = " << debug_sc.at(1,1,2) << endl;
//		cout << "debug_sc.at(0,y1,x1 + 1) = " << debug_sc.at(0,1,0) << endl;
//	}
//		if (ay==0&&ax==0)		cout << "pars.imageSizeReduce[0] = " << pars.imageSizeReduce[0] << endl;
//		if (ay==0&&ax==0)		cout << "pars.imageSizeReduce[1] = " << pars.imageSizeReduce[1] << endl;

		// re-center the indices
		zeroIndices <<<(pars.imageSizeReduce[0] - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream>>> (
				y_ds, pars.imageSizeReduce[0]);

		zeroIndices <<<(pars.imageSizeReduce[1] - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream>>> (
				x_ds, pars.imageSizeReduce[1]);

//		if (ax ==0 && ay == 0) {
//			long y_ds_t, x_ds_t;
//			cudaErrchk(cudaMemcpy(&y_ds_t, y_ds, sizeof(long), cudaMemcpyDeviceToHost));
//			cudaErrchk(cudaMemcpy(&x_ds_t, x_ds, sizeof(long), cudaMemcpyDeviceToHost));
//			cout << "y_ds_t = " << y_ds_t << endl;
//			cout << "x_ds_t = " << x_ds_t << endl;
//		}
		// Choose a good launch configuration
		// Heuristically use 2^p / 2 as the block size where p is the first power of 2 greater than the number of elements to work on.
		// This balances having enough work per thread and enough blocks without having so many blocks that the shared memory doesn't last long
		size_t p = getNextPower2(pars.numberBeams);
		const size_t BlockSizeX = (size_t)std::max(1.0, pow(2,p) / 2);

		// Determine maximum threads per streaming multiprocessor based on the compute capability of the device
		size_t max_threads_per_sm;
		if (pars.deviceProperties.major > 3){
			max_threads_per_sm = 2048;
		} else if (pars.deviceProperties.major > 2) {
			max_threads_per_sm = 1536;
		} else {
			max_threads_per_sm = pars.deviceProperties.minor == 0 ? 768 : 1024;
		}

		// Estimate max number of blocks per streaming multiprocessor (threads are the limit)
		const size_t max_blocks_per_sm = max_threads_per_sm / BlockSizeX;

		// We find providing around 3 times as many blocks as the estimated maximum provides good performance
		const size_t target_blocks_per_sm = max_blocks_per_sm * 3;
		const size_t total_blocks         = target_blocks_per_sm * pars.deviceProperties.multiProcessorCount;

		// Determine the shape of the grid
		const PRISM_FLOAT_PRECISION aspect_ratio = pars.imageSizeReduce[1] / pars.imageSizeReduce[0];
		const size_t BlockSizeZ = std::floor(sqrt(total_blocks / aspect_ratio));
		const size_t BlockSizeY = aspect_ratio * BlockSizeZ;

//		if (ax == 0 & ay == 0) {
//			cout << "aspect_ratio = " << aspect_ratio << endl;
//			cout << "BlockSizeX = " << BlockSizeX << endl;
//			cout << "BlockSizeZ = " << BlockSizeZ << endl;
//			cout << "BlockSizeY = " << BlockSizeY << endl;
//			cout << "target_blocks_per_sm = " << target_blocks_per_sm << endl;
//			cout << "total_blocks = " << total_blocks << endl;
//			cout << " pars.deviceProperties.multiProcessorCount = " <<  pars.deviceProperties.multiProcessorCount << endl;
//		}
		dim3 grid(1, BlockSizeY, BlockSizeZ);
		dim3 block(BlockSizeX, 1, 1);
//		dim3 block(4, 1, 1);

		// Determine amount of shared memory needed
		const unsigned long smem = pars.numberBeams * sizeof(PRISM_CUDA_COMPLEX_FLOAT);

//		scaleReduceS<4> <<< grid, block, smem, stream >>> (
//				permuted_Scompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
//				pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]);

		// Launch kernel. Block size must be visible at compile time so we use a switch statement
		switch (BlockSizeX) {
			case 1024 :
				scaleReduceS<1024> <<< grid, block, smem, stream >>> (
						permuted_Scompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
						pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
			case 512 :
				scaleReduceS<512> <<< grid, block, smem, stream >>> (
						permuted_Scompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
						pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
			case 256 :
				scaleReduceS<256> <<< grid, block, smem, stream >>> (
						permuted_Scompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
						pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
			case 128 :
				scaleReduceS<128> <<< grid, block, smem, stream >>> (
						permuted_Scompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
						pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
			case 64 :
				scaleReduceS<64> <<< grid, block, smem, stream >>> (
						permuted_Scompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
						pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
			case 32 :
				scaleReduceS<32> <<< grid, block, smem, stream >>> (
						permuted_Scompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
						pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
			case 16 :
				scaleReduceS<16> <<< grid, block, smem, stream >>> (
						permuted_Scompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
						pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
			case 8 :
				scaleReduceS<8> <<< grid, block, smem, stream >>> (
						permuted_Scompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
						pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
			case 4 :
				scaleReduceS<4> <<< grid, block, smem, stream >>> (
						permuted_Scompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
						pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
			case 2 :
				scaleReduceS<2> <<< grid, block, smem, stream >>> (
						permuted_Scompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
						pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
			default :
				scaleReduceS<1> <<< grid, block, smem, stream >>> (
						permuted_Scompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
						pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
		}

		/*
		{
////			cudaDeviceSynchronize();
			complex<PRISM_FLOAT_PRECISION > ans;
			for (auto i = 0; i < pars.imageSizeReduce[0]; ++i){
				cudaErrchk(cudaMemcpy(&ans, (psi_ds + i), sizeof(complex<PRISM_FLOAT_PRECISION>),cudaMemcpyDeviceToHost));
				cout << "psi_ds[" << i << "] = " << ans << endl;
			}
		}
		*/

		// final fft
		cufftErrchk(PRISM_CUFFT_EXECUTE(cufft_plan, &psi_ds[0], &psi_ds[0], CUFFT_FORWARD));

		// convert to squared intensity
		abs_squared <<< (psi_size - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >>> (psi_intensity_ds, psi_ds, psi_size);

		// output calculation result
		formatOutput_GPU_integrate(pars, psi_intensity_ds, alphaInd_d, output_ph,
		                           integratedOutput_ds, ay, ax, pars.imageSizeReduce[0],
		                           pars.imageSizeReduce[1], stream, pars.scale);

	}










}
