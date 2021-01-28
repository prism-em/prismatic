// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// Prismatic is distributed under the GNU General Public License (GPL)
// If you use Prismatic, we kindly ask that you cite the following papers:

// 1. Ophus, C.: A fast image simulation algorithm for scanning
//    transmission electron microscopy. Advanced Structural and
//    Chemical Imaging 3(1), 13 (2017)

// 2. Pryor, Jr., A., Ophus, C., and Miao, J.: A Streaming Multi-GPU
//    Implementation of Image Simulation Algorithms for Scanning
//	  Transmission Electron Microscopy. arXiv:1706.08563 (2017)

// For variable naming, the suffixes are "_d" for "device" (1 copy per GPU), "_ds" for "device stream (1 copy per stream), "_ph" for "pinned host"

#include <iostream>
#include <complex>
#include <thread>
#include <vector>
#include "WorkDispatcher.h"
#include "PRISM03_calcOutput.cuh"
#include "PRISM03_calcOutput.h"
#include "utility.cuh"
#include "fileIO.cuh"

namespace Prismatic {
	extern std::mutex fftw_plan_lock;

	// define some constants
	__device__ __constant__ float pi_f = PI;
	__device__ __constant__ cuFloatComplex i_f = {0, 1};
	__device__ __constant__ cuFloatComplex pi_cx_f = {PI, 0};
	__device__ __constant__ cuFloatComplex minus_2pii_f = {0, -2 * PI};
	__device__ __constant__ double pi = PI;
	__device__ __constant__ cuDoubleComplex i = {0, 1};
	__device__ __constant__ cuDoubleComplex pi_cx = {PI, 0};
	__device__ __constant__ cuDoubleComplex minus_2pii = {0, -2 * PI};

	inline void createStreamsAndPlans3(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
                                      CudaParameters<PRISMATIC_FLOAT_PRECISION> &cuda_pars) {
		const int total_num_streams = pars.meta.numGPUs * pars.meta.numStreamsPerGPU;

		// create CUDA streams and cuFFT plans
		cuda_pars.streams = new cudaStream_t[total_num_streams];;
		cuda_pars.cufftPlans = new cufftHandle[total_num_streams];
		
		for (auto j = 0; j < total_num_streams; ++j) {
			cudaSetDevice(j % pars.meta.numGPUs);
			cudaErrchk(cudaStreamCreate(&cuda_pars.streams[j]));
			cufftErrchk(cufftPlan2d(&cuda_pars.cufftPlans[j], pars.imageSizeReduce[0], pars.imageSizeReduce[1], PRISMATIC_CUFFT_PLAN_TYPE));
			cufftErrchk(cufftSetStream(cuda_pars.cufftPlans[j], cuda_pars.streams[j]));
		}
	}

	inline void allocatePinnedHostMemory_singlexfer3(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
	                                          CudaParameters<PRISMATIC_FLOAT_PRECISION> &cuda_pars) {
		// Allocate pinned memory buffers 

		const int total_num_streams = pars.meta.numGPUs * pars.meta.numStreamsPerGPU;
		cuda_pars.output_ph = new PRISMATIC_FLOAT_PRECISION *[total_num_streams]; // one output array per stream
		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaMallocHost((void **) &cuda_pars.output_ph[s],     pars.output.get_dimi()   * sizeof(PRISMATIC_FLOAT_PRECISION)));
		}
		cudaErrchk(cudaMallocHost((void **) &cuda_pars.permutedScompact_ph, pars.Scompact.size()     * sizeof(std::complex<PRISMATIC_FLOAT_PRECISION>)));
		cudaErrchk(cudaMallocHost((void **) &cuda_pars.PsiProbeInit_ph,      pars.psiProbeInit.size() * sizeof(std::complex<PRISMATIC_FLOAT_PRECISION>)));
		cudaErrchk(cudaMallocHost((void **) &cuda_pars.qxaReduce_ph,         pars.qxaReduce.size()    * sizeof(PRISMATIC_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **) &cuda_pars.qyaReduce_ph,         pars.qyaReduce.size()    * sizeof(PRISMATIC_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **) &cuda_pars.alphaInd_ph,          pars.alphaInd.size()     * sizeof(PRISMATIC_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **) &cuda_pars.xBeams_ph,            pars.xyBeams.get_dimj()  * sizeof(size_t)));
		cudaErrchk(cudaMallocHost((void **) &cuda_pars.yBeams_ph,            pars.xyBeams.get_dimj()  * sizeof(size_t)));
	}

	inline void allocatePinnedHostMemory_streaming3(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
	                                         CudaParameters<PRISMATIC_FLOAT_PRECISION> &cuda_pars) {

		const int total_num_streams = pars.meta.numGPUs * pars.meta.numStreamsPerGPU;
		cuda_pars.output_ph = new PRISMATIC_FLOAT_PRECISION *[total_num_streams]; // one output array per stream
		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaMallocHost((void **) &cuda_pars.output_ph[s],
			pars.output.get_dimi() * sizeof(PRISMATIC_FLOAT_PRECISION)));
		}

		// allocate pinned memory
		cudaErrchk(cudaMallocHost((void **) &cuda_pars.permutedScompact_ph, pars.Scompact.size()     * sizeof(std::complex<PRISMATIC_FLOAT_PRECISION>)));
		cudaErrchk(cudaMallocHost((void **) &cuda_pars.PsiProbeInit_ph,      pars.psiProbeInit.size() * sizeof(std::complex<PRISMATIC_FLOAT_PRECISION>)));
		cudaErrchk(cudaMallocHost((void **) &cuda_pars.qxaReduce_ph,         pars.qxaReduce.size()    * sizeof(PRISMATIC_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **) &cuda_pars.qyaReduce_ph,         pars.qyaReduce.size()    * sizeof(PRISMATIC_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **) &cuda_pars.alphaInd_ph,          pars.alphaInd.size()     * sizeof(PRISMATIC_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **) &cuda_pars.xBeams_ph,            pars.xyBeams.get_dimj()  * sizeof(size_t)));
		cudaErrchk(cudaMallocHost((void **) &cuda_pars.yBeams_ph,            pars.xyBeams.get_dimj()  * sizeof(size_t)));
	}

	inline void copyToPinnedMemory_singlexfer3(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
	                                          CudaParameters<PRISMATIC_FLOAT_PRECISION> &cuda_pars) {
		// Copy data to pinned memory buffers

		const int total_num_streams = pars.meta.numGPUs * pars.meta.numStreamsPerGPU;
		// copy host memory to pinned
		for (auto s = 0; s < total_num_streams; ++s) {
			memset(cuda_pars.output_ph[s], 0, pars.output.get_dimi() * sizeof(PRISMATIC_FLOAT_PRECISION));
		}

		// Copy to pinned memory
		// The GPU computational model operates on Scompact in a different order than the CPU, and it is
		// more optimal to permute the dimensions so that the consecutive elements represent different
		// beams on the GPU as opposed to consecutive x-probe positions on the CPU
		{
			auto S_ptr = cuda_pars.permutedScompact_ph;
			for (auto jj = 0; jj < pars.Scompact.get_dimj(); ++jj) {
				for (auto ii = 0; ii < pars.Scompact.get_dimi(); ++ii) {
					for (auto kk = 0; kk < pars.Scompact.get_dimk(); ++kk) {
						*S_ptr++ = pars.Scompact.at(kk, jj, ii);
					}
				}
			}
		}

		memcpy(cuda_pars.PsiProbeInit_ph, &(*pars.psiProbeInit.begin()), pars.psiProbeInit.size() * sizeof(std::complex<PRISMATIC_FLOAT_PRECISION>));
		memcpy(cuda_pars.alphaInd_ph,     &(*pars.alphaInd.begin()),     pars.alphaInd.size()     * sizeof(PRISMATIC_FLOAT_PRECISION));
		memcpy(cuda_pars.qxaReduce_ph,    &(*pars.qxaReduce.begin()),    pars.qxaReduce.size()    * sizeof(PRISMATIC_FLOAT_PRECISION));
		memcpy(cuda_pars.qyaReduce_ph,    &(*pars.qyaReduce.begin()),    pars.qyaReduce.size()    * sizeof(PRISMATIC_FLOAT_PRECISION));

		{
			auto x_ptr = cuda_pars.xBeams_ph;
			auto y_ptr = cuda_pars.yBeams_ph;
			for (auto jj = 0; jj < pars.xyBeams.get_dimj(); ++jj) {
				*y_ptr++ = pars.xyBeams.at(jj, 0);
				*x_ptr++ = pars.xyBeams.at(jj, 1);
			}
		}
	}

	inline void copyToPinnedMemory_streaming3(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
	                                         CudaParameters<PRISMATIC_FLOAT_PRECISION> &cuda_pars) {
		// Copy data to pinned memory buffers

		const int total_num_streams = pars.meta.numGPUs * pars.meta.numStreamsPerGPU;
		// copy host memory to pinned
		for (auto s = 0; s < total_num_streams; ++s) {
			memset(cuda_pars.output_ph[s], 0, pars.output.get_dimi() * sizeof(PRISMATIC_FLOAT_PRECISION));
		}

		// Copy to pinned memory
		// The GPU computational model operates on Scompact in a different order than the CPU, and it is
		// more optimal to permute the dimensions so that the consecutive elements represent different
		// beams on the GPU as opposed to consecutive x-probe positions on the CPU
		{
			auto S_ptr = cuda_pars.permutedScompact_ph;
			for (auto jj = 0; jj < pars.Scompact.get_dimj(); ++jj){
				for (auto ii = 0; ii < pars.Scompact.get_dimi(); ++ii){
					for (auto kk = 0; kk < pars.Scompact.get_dimk(); ++kk){
						*S_ptr++ = pars.Scompact.at(kk, jj, ii);
					}
				}
			}
		}

		memcpy(cuda_pars.PsiProbeInit_ph, &(*pars.psiProbeInit.begin()), pars.psiProbeInit.size() * sizeof(std::complex<PRISMATIC_FLOAT_PRECISION>));
		memcpy(cuda_pars.alphaInd_ph,     &(*pars.alphaInd.begin()),     pars.alphaInd.size()     * sizeof(PRISMATIC_FLOAT_PRECISION));
		memcpy(cuda_pars.qxaReduce_ph,    &(*pars.qxaReduce.begin()),    pars.qxaReduce.size()    * sizeof(PRISMATIC_FLOAT_PRECISION));
		memcpy(cuda_pars.qyaReduce_ph,    &(*pars.qyaReduce.begin()),    pars.qyaReduce.size()    * sizeof(PRISMATIC_FLOAT_PRECISION));

		{
			auto x_ptr = cuda_pars.xBeams_ph;
			auto y_ptr = cuda_pars.yBeams_ph;
			for (auto jj = 0; jj < pars.xyBeams.get_dimj(); ++jj){
				*y_ptr++ = pars.xyBeams.at(jj,0);
				*x_ptr++ = pars.xyBeams.at(jj,1);
			}
		}

	}

	inline void allocateDeviceMemory_singlexfer3(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
	                                         CudaParameters<PRISMATIC_FLOAT_PRECISION> &cuda_pars) {
		const int total_num_streams = pars.meta.numGPUs * pars.meta.numStreamsPerGPU;

		// pointers to read-only GPU memory (one copy per GPU)
		cuda_pars.permutedScompact_d = new PRISMATIC_CUDA_COMPLEX_FLOAT *[pars.meta.numGPUs];
		cuda_pars.PsiProbeInit_d      = new PRISMATIC_CUDA_COMPLEX_FLOAT *[pars.meta.numGPUs];
		cuda_pars.qxaReduce_d         = new PRISMATIC_FLOAT_PRECISION *[pars.meta.numGPUs];
		cuda_pars.qyaReduce_d         = new PRISMATIC_FLOAT_PRECISION *[pars.meta.numGPUs];
		cuda_pars.alphaInd_d          = new PRISMATIC_FLOAT_PRECISION *[pars.meta.numGPUs];
		cuda_pars.yBeams_d            = new size_t *[pars.meta.numGPUs];
		cuda_pars.xBeams_d            = new size_t *[pars.meta.numGPUs];

		// pointers to read/write GPU memory (one copy per stream)
		cuda_pars.psi_ds              = new PRISMATIC_CUDA_COMPLEX_FLOAT *[total_num_streams];
		cuda_pars.phaseCoeffs_ds      = new PRISMATIC_CUDA_COMPLEX_FLOAT *[total_num_streams];
		cuda_pars.y_ds                = new long *[total_num_streams];
		cuda_pars.x_ds                = new long *[total_num_streams];
		cuda_pars.psiIntensity_ds    = new PRISMATIC_FLOAT_PRECISION *[total_num_streams];
		cuda_pars.integratedOutput_ds = new PRISMATIC_FLOAT_PRECISION *[total_num_streams];

		// allocate memory on each GPU
		for (auto g = 0; g < pars.meta.numGPUs; ++g) {
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.permutedScompact_d[g], pars.Scompact.size()     * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.PsiProbeInit_d[g],      pars.psiProbeInit.size() * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.qxaReduce_d[g],         pars.qxaReduce.size()    * sizeof(PRISMATIC_FLOAT_PRECISION)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.qyaReduce_d[g],         pars.qyaReduce.size()    * sizeof(PRISMATIC_FLOAT_PRECISION)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.alphaInd_d[g],          pars.alphaInd.size()     * sizeof(PRISMATIC_FLOAT_PRECISION)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.yBeams_d[g],            pars.xyBeams.get_dimj()  * sizeof(size_t)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.xBeams_d[g],            pars.xyBeams.get_dimj()  * sizeof(size_t)));
		}

		// allocate memory per stream and zero it
		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaSetDevice(s % pars.meta.numGPUs));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.psi_ds[s],              pars.imageSizeReduce[0]    * pars.imageSizeReduce[1] * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.phaseCoeffs_ds[s],      pars.numberBeams           * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.y_ds[s],                pars.imageSizeReduce[0]    * sizeof(long)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.x_ds[s],                pars.imageSizeReduce[1]    * sizeof(long)));

			cudaErrchk(cudaMemset(cuda_pars.psi_ds[s],              0, pars.imageSizeReduce[0]    * pars.imageSizeReduce[1] * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMemset(cuda_pars.phaseCoeffs_ds[s],      0, pars.numberBeams           * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMemset(cuda_pars.y_ds[s],                0, pars.imageSizeReduce[0]    * sizeof(long)));
			cudaErrchk(cudaMemset(cuda_pars.x_ds[s],                0, pars.imageSizeReduce[1]    * sizeof(long)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.psiIntensity_ds[s],    pars.imageSizeReduce[0]    * pars.imageSizeReduce[1] * sizeof(PRISMATIC_FLOAT_PRECISION)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.integratedOutput_ds[s], pars.detectorAngles.size() * sizeof(PRISMATIC_FLOAT_PRECISION)));
			cudaErrchk(cudaMemset(cuda_pars.psiIntensity_ds[s],    0, pars.imageSizeReduce[0]    * pars.imageSizeReduce[1] * sizeof(PRISMATIC_FLOAT_PRECISION)));
			cudaErrchk(cudaMemset(cuda_pars.integratedOutput_ds[s], 0, pars.detectorAngles.size() * sizeof(PRISMATIC_FLOAT_PRECISION)));
		}
	}

	inline void allocateDeviceMemory_streaming3(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
	                                        CudaParameters<PRISMATIC_FLOAT_PRECISION> &cuda_pars) {
		const int total_num_streams = pars.meta.numGPUs * pars.meta.numStreamsPerGPU;

		// pointers to read-only GPU memory (one copy per GPU)
		cuda_pars.PsiProbeInit_d = new PRISMATIC_CUDA_COMPLEX_FLOAT*[pars.meta.numGPUs];
		cuda_pars.qxaReduce_d    = new PRISMATIC_FLOAT_PRECISION*[pars.meta.numGPUs];
		cuda_pars.qyaReduce_d    = new PRISMATIC_FLOAT_PRECISION*[pars.meta.numGPUs];
		cuda_pars.alphaInd_d 	 = new PRISMATIC_FLOAT_PRECISION*[pars.meta.numGPUs];
		cuda_pars.yBeams_d 	     = new size_t*[pars.meta.numGPUs];
		cuda_pars.xBeams_d 	     = new size_t*[pars.meta.numGPUs];

		// pointers to read/write GPU memory (one per stream)
		cuda_pars.permutedScompact_d   = new PRISMATIC_CUDA_COMPLEX_FLOAT*[total_num_streams];
		cuda_pars.psi_ds 				= new PRISMATIC_CUDA_COMPLEX_FLOAT*[total_num_streams];
		cuda_pars.phaseCoeffs_ds 		= new PRISMATIC_CUDA_COMPLEX_FLOAT*[total_num_streams];
		cuda_pars.y_ds                  = new long*[total_num_streams];
		cuda_pars.x_ds                  = new long*[total_num_streams];
		cuda_pars.psiIntensity_ds    = new PRISMATIC_FLOAT_PRECISION *[total_num_streams];
		cuda_pars.integratedOutput_ds = new PRISMATIC_FLOAT_PRECISION *[total_num_streams];

		// allocate memory on each GPU
		for (auto g = 0; g < pars.meta.numGPUs; ++g) {
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.PsiProbeInit_d[g], pars.psiProbeInit.size()  * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.qxaReduce_d[g],    pars.qxaReduce.size()     * sizeof(PRISMATIC_FLOAT_PRECISION)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.qyaReduce_d[g],    pars.qyaReduce.size()     * sizeof(PRISMATIC_FLOAT_PRECISION)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.alphaInd_d[g],     pars.alphaInd.size()      * sizeof(PRISMATIC_FLOAT_PRECISION)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.yBeams_d[g],       pars.xyBeams.get_dimj()   * sizeof(size_t)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.xBeams_d[g],       pars.xyBeams.get_dimj()   * sizeof(size_t)));
		}

		// allocate memory per stream and zero it
		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaSetDevice(s % pars.meta.numGPUs));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.permutedScompact_d[s], pars.imageSizeReduce[0] * pars.imageSizeReduce[1] * pars.numberBeams * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.psi_ds[s],              pars.imageSizeReduce[0] * pars.imageSizeReduce[1] * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.phaseCoeffs_ds[s],      pars.numberBeams * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.y_ds[s],                pars.imageSizeReduce[0] * sizeof(long)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.x_ds[s],                pars.imageSizeReduce[1] * sizeof(long)));
			
			cudaErrchk(cudaMemset(cuda_pars.psi_ds[s],              0, pars.imageSizeReduce[0]    * pars.imageSizeReduce[1] * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMemset(cuda_pars.phaseCoeffs_ds[s],      0, pars.numberBeams           * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMemset(cuda_pars.y_ds[s],                0, pars.imageSizeReduce[0]    * sizeof(long)));
			cudaErrchk(cudaMemset(cuda_pars.x_ds[s],                0, pars.imageSizeReduce[1]    * sizeof(long)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.psiIntensity_ds[s],    pars.imageSizeReduce[0]    * pars.imageSizeReduce[1] * sizeof(PRISMATIC_FLOAT_PRECISION)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.integratedOutput_ds[s], pars.detectorAngles.size() * sizeof(PRISMATIC_FLOAT_PRECISION)));
			cudaErrchk(cudaMemset(cuda_pars.psiIntensity_ds[s],    0, pars.imageSizeReduce[0]    * pars.imageSizeReduce[1] * sizeof(PRISMATIC_FLOAT_PRECISION)));
			cudaErrchk(cudaMemset(cuda_pars.integratedOutput_ds[s], 0, pars.detectorAngles.size() * sizeof(PRISMATIC_FLOAT_PRECISION)));
		}
	}

	inline void copyToGPUMemory_singlexfer3(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
	                                       CudaParameters<PRISMATIC_FLOAT_PRECISION> &cuda_pars) {

		// Copy memory to each GPU asynchronously from the pinned host memory spaces.
		// The streams are laid out so that consecutive streams represent different GPUs. If we
		// have more than one stream per GPU, then we try to interleave as much as possible
		const int total_num_streams = pars.meta.numGPUs * pars.meta.numStreamsPerGPU;
		int stream_id = 0;
		for (auto g = 0; g < pars.meta.numGPUs; ++g) {
			stream_id = g;
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaMemcpyAsync(cuda_pars.permutedScompact_d[g], &cuda_pars.permutedScompact_ph[0], pars.Scompact.size() * sizeof(std::complex<PRISMATIC_FLOAT_PRECISION>), cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));

			stream_id = (stream_id + pars.meta.numGPUs) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(cuda_pars.PsiProbeInit_d[g], &cuda_pars.PsiProbeInit_ph[0], pars.psiProbeInit.size() * sizeof(std::complex<PRISMATIC_FLOAT_PRECISION>), cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));

			stream_id = (stream_id + pars.meta.numGPUs) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(cuda_pars.qxaReduce_d[g], &cuda_pars.qxaReduce_ph[0], pars.qxaReduce.size() * sizeof(PRISMATIC_FLOAT_PRECISION), cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));

			stream_id = (stream_id + pars.meta.numGPUs) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(cuda_pars.qyaReduce_d[g], &cuda_pars.qyaReduce_ph[0], pars.qyaReduce.size() * sizeof(PRISMATIC_FLOAT_PRECISION), cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));

			stream_id = (stream_id + pars.meta.numGPUs) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(cuda_pars.alphaInd_d[g], &cuda_pars.alphaInd_ph[0], pars.alphaInd.size() * sizeof(PRISMATIC_FLOAT_PRECISION), cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));

			stream_id = (stream_id + pars.meta.numGPUs) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(cuda_pars.yBeams_d[g], &cuda_pars.yBeams_ph[0], pars.xyBeams.get_dimj() * sizeof(size_t), cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));

			stream_id = (stream_id + pars.meta.numGPUs) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(cuda_pars.xBeams_d[g], &cuda_pars.xBeams_ph[0], pars.xyBeams.get_dimj() * sizeof(size_t), cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));

		}
		// wait for transfers to complete
		for (auto g = 0; g < pars.meta.numGPUs; ++g) {
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaDeviceSynchronize());
		}
	}

	inline void copyToGPUMemory_streaming3(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
	                                      CudaParameters<PRISMATIC_FLOAT_PRECISION> &cuda_pars) {

		// Copy memory to each GPU asynchronously from the pinned host memory spaces.
		// The streams are laid out so that consecutive streams represent different GPUs. If we
		// have more than one stream per GPU, then we try to interleave as much as possible
		const int total_num_streams = pars.meta.numGPUs * pars.meta.numStreamsPerGPU;
		int stream_id = 0;
		for (auto g = 0; g < pars.meta.numGPUs; ++g) {
			stream_id = g;
			cudaErrchk(cudaSetDevice(g));

			stream_id = (stream_id + pars.meta.numGPUs) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(cuda_pars.PsiProbeInit_d[g], &cuda_pars.PsiProbeInit_ph[0], pars.psiProbeInit.size() * sizeof(std::complex<PRISMATIC_FLOAT_PRECISION>), cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));

			stream_id = (stream_id + pars.meta.numGPUs) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(cuda_pars.qxaReduce_d[g], &cuda_pars.qxaReduce_ph[0], pars.qxaReduce.size() * sizeof(PRISMATIC_FLOAT_PRECISION), cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));

			stream_id = (stream_id + pars.meta.numGPUs) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(cuda_pars.qyaReduce_d[g], &cuda_pars.qyaReduce_ph[0], pars.qyaReduce.size() * sizeof(PRISMATIC_FLOAT_PRECISION), cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));

			stream_id = (stream_id + pars.meta.numGPUs) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(cuda_pars.alphaInd_d[g], &cuda_pars.alphaInd_ph[0], pars.alphaInd.size() * sizeof(pars.alphaInd[0]), cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));

			stream_id = (stream_id + pars.meta.numGPUs) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(cuda_pars.yBeams_d[g], &cuda_pars.yBeams_ph[0], pars.xyBeams.get_dimj() * sizeof(size_t), cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));

			stream_id = (stream_id + pars.meta.numGPUs) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(cuda_pars.xBeams_d[g], &cuda_pars.xBeams_ph[0], pars.xyBeams.get_dimj() * sizeof(size_t), cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));
		}

		// wait for transfers to complete
		for (auto g = 0; g < pars.meta.numGPUs; ++g) {
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaDeviceSynchronize());
		}
	}

	inline void launchWorkers_singlexfer3(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
	                                     CudaParameters<PRISMATIC_FLOAT_PRECISION> &cuda_pars) {
		using namespace std;
		const int total_num_streams = pars.meta.numGPUs * pars.meta.numStreamsPerGPU;
		// launch threads that will consume work provided by getWorkID
		vector<thread> workers_GPU;
		workers_GPU.reserve(total_num_streams); // prevents multiple reallocations
		int stream_count = 0;
		const size_t PRISMATIC_PRINT_FREQUENCY_PROBES = max((size_t) 1, pars.numProbes / 10); // for printing status
		WorkDispatcher dispatcher(0, pars.numProbes); // create work dispatcher

		for (auto t = 0; t < total_num_streams; ++t) {

			int GPU_num = stream_count % pars.meta.numGPUs; // determine which GPU handles this job
			cudaStream_t &current_stream = cuda_pars.streams[stream_count];
			cout << "Launching GPU worker on stream #" << stream_count << " of GPU #" << GPU_num << endl;

			// get pointers to the pre-copied arrays, making sure to get those on the current GPU
			PRISMATIC_CUDA_COMPLEX_FLOAT *current_permutedScompact_d = cuda_pars.permutedScompact_d[GPU_num];
			PRISMATIC_CUDA_COMPLEX_FLOAT *current_PsiProbeInit_d      = cuda_pars.PsiProbeInit_d[GPU_num];
			PRISMATIC_FLOAT_PRECISION *current_qxaReduce_d            = cuda_pars.qxaReduce_d[GPU_num];
			PRISMATIC_FLOAT_PRECISION *current_qyaReduce_d            = cuda_pars.qyaReduce_d[GPU_num];
			size_t *current_yBeams_d                              = cuda_pars.yBeams_d[GPU_num];
			size_t *current_xBeams_d                              = cuda_pars.xBeams_d[GPU_num];
			PRISMATIC_FLOAT_PRECISION *current_alphaInd_d             = cuda_pars.alphaInd_d[GPU_num];

			// get pointers to per-stream arrays
			PRISMATIC_CUDA_COMPLEX_FLOAT *current_psi_ds           = cuda_pars.psi_ds[stream_count];
			PRISMATIC_CUDA_COMPLEX_FLOAT *current_phaseCoeffs_ds   = cuda_pars.phaseCoeffs_ds[stream_count];
			long *current_y_ds                                 = cuda_pars.y_ds[stream_count];
			long *current_x_ds                                 = cuda_pars.x_ds[stream_count];
			cufftHandle &current_cufft_plan                    = cuda_pars.cufftPlans[stream_count];
			// get pointer to output pinned memory
			PRISMATIC_FLOAT_PRECISION *current_psiIntensity_ds  = cuda_pars.psiIntensity_ds[stream_count];
			PRISMATIC_FLOAT_PRECISION *current_integratedOutput_ds = cuda_pars.integratedOutput_ds[stream_count];
			PRISMATIC_FLOAT_PRECISION *current_output_ph = cuda_pars.output_ph[stream_count];

			workers_GPU.push_back(
					thread([&pars, GPU_num, stream_count, current_permutedScompact_d, &dispatcher, current_alphaInd_d, current_PsiProbeInit_d, current_qxaReduce_d, current_qyaReduce_d,
							       current_yBeams_d, current_xBeams_d, current_psi_ds, current_phaseCoeffs_ds,
								   current_psiIntensity_ds, current_y_ds, current_x_ds, current_integratedOutput_ds, current_output_ph,
							       &current_cufft_plan, &current_stream, &PRISMATIC_PRINT_FREQUENCY_PROBES, &cuda_pars]() {
						cudaErrchk(cudaSetDevice(GPU_num));


#ifndef NDEBUG
						{
//					 check memory usage on the GPU
							std::lock_guard<mutex> lock(Prismatic::memLock);
							size_t free_mem, total_mem;
							free_mem = total_mem = 0;
							cudaErrchk(cudaMemGetInfo(&free_mem, &total_mem));
							pars.maxGPUMem = std::max(total_mem - free_mem, pars.maxGPUMem);
//					cout << "maxGPUMem = " << pars.maxGPUMem << endl;
						}
#endif // NDEBUG

						// main work loop
						size_t Nstart, Nstop, ay, ax;
						Nstart = Nstop = 0;
						while (dispatcher.getWork(Nstart, Nstop)) { // synchronously get work assignment
							while (Nstart < Nstop) {
								if (Nstart % PRISMATIC_PRINT_FREQUENCY_PROBES == 0 | Nstart == 100) {
									cout << "Computing Probe Position #" << Nstart << "/"
									     << pars.numProbes << endl;
								}
								ay = (pars.meta.arbitraryProbes) ? Nstart : Nstart / pars.numXprobes;
								ax = (pars.meta.arbitraryProbes) ? Nstart : Nstart % pars.numXprobes;
								buildSignal_GPU_singlexfer(pars, ay, ax, current_permutedScompact_d,
									current_PsiProbeInit_d, current_qxaReduce_d,
									current_qyaReduce_d,
									current_yBeams_d, current_xBeams_d, current_alphaInd_d,
									current_psi_ds,
									current_phaseCoeffs_ds, current_psiIntensity_ds,
									current_y_ds,
									current_x_ds, current_output_ph, current_integratedOutput_ds,
									current_cufft_plan, current_stream, cuda_pars);
#ifdef PRISMATIC_BUILDING_GUI
								pars.progressbar->signalOutputUpdate(Nstart, pars.numProbes);
#endif
//						buildSignal_CPU(pars, ay, ax, yTiltShift, xTiltShift, alphaInd, PsiProbeInit);
								++Nstart;
							}
						}
					}));
			++stream_count;
		}


		// Now launch CPU work
		if (pars.meta.alsoDoCPUWork) {
			PRISMATIC_FFTW_INIT_THREADS();
			PRISMATIC_FFTW_PLAN_WITH_NTHREADS(pars.meta.numThreads);
			vector <thread> workers_CPU;
			workers_CPU.reserve(pars.meta.numThreads); // prevents multiple reallocations
			for (auto t = 0; t < pars.meta.numThreads; ++t) {
				cout << "Launching CPU worker thread #" << t << " to compute partial PRISM result\n";
				workers_CPU.push_back(thread([&pars, &dispatcher, t, &PRISMATIC_PRINT_FREQUENCY_PROBES]() {
					size_t Nstart, Nstop, ay, ax, early_CPU_stop;
					Nstart = Nstop = 0;
					if (pars.meta.numGPUs > 0) {
						// if there are no GPUs, make sure to do all work on CPU
						early_CPU_stop = (size_t) std::max((PRISMATIC_FLOAT_PRECISION) 0.0,
						                                   pars.numProbes - pars.meta.earlyCPUStopCount);
					} else {
						early_CPU_stop = pars.numProbes;
					}
					cout << "early_CPU_stop= " << early_CPU_stop << endl;
					if (dispatcher.getWork(Nstart, Nstop, 1, early_CPU_stop)) { // synchronously get work assignment
						Array2D <std::complex<PRISMATIC_FLOAT_PRECISION>> psi = Prismatic::zeros_ND<2, std::complex<PRISMATIC_FLOAT_PRECISION> >(
								{{pars.imageSizeReduce[0], pars.imageSizeReduce[1]}});
						unique_lock <mutex> gatekeeper(fftw_plan_lock);
						PRISMATIC_FFTW_PLAN plan = PRISMATIC_FFTW_PLAN_DFT_2D(psi.get_dimj(), psi.get_dimi(),
						                                              reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&psi[0]),
						                                              reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&psi[0]),
						                                              FFTW_FORWARD, FFTW_MEASURE);
						gatekeeper.unlock();

						// main work loop
						do {
							while (Nstart < Nstop) {
								if (Nstart % PRISMATIC_PRINT_FREQUENCY_PROBES == 0 | Nstart == 100) {
									cout << "Computing Probe Position #" << Nstart << "/"
									     << pars.numProbes << endl;
								}
								ay = (pars.meta.arbitraryProbes) ? Nstart : Nstart / pars.numXprobes;
								ax = (pars.meta.arbitraryProbes) ? Nstart : Nstart % pars.numXprobes;
								buildSignal_CPU(pars, ay, ax, plan, psi);
#ifdef PRISMATIC_BUILDING_GUI
								pars.progressbar->signalOutputUpdate(Nstart, pars.numProbes);
#endif
								++Nstart;
							}
							if (Nstop >= early_CPU_stop) break;
						} while (dispatcher.getWork(Nstart, Nstop, 1, early_CPU_stop));
						gatekeeper.lock();
						PRISMATIC_FFTW_DESTROY_PLAN(plan);
						gatekeeper.unlock();
					}
				}));
			}
			cout << "Waiting for CPU threads...\n";
			for (auto &t:workers_CPU)t.join();
			PRISMATIC_FFTW_CLEANUP_THREADS();
		}

		// synchronize
		cout << "Waiting for GPU threads...\n";
		for (auto &t:workers_GPU)t.join();

		for (auto g = 0; g < pars.meta.numGPUs; ++g){
			cudaSetDevice(g);
			cudaDeviceSynchronize();
		}
	}

	void launchWorkers_streaming3(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
	                             CudaParameters<PRISMATIC_FLOAT_PRECISION> &cuda_pars) {
		using namespace std;
		const int total_num_streams = pars.meta.numGPUs * pars.meta.numStreamsPerGPU;

		// launch threads that will consume work provided by getWorkID
		vector<thread> workers_GPU;
		workers_GPU.reserve(total_num_streams); // prevents multiple reallocations
		int stream_count = 0;
		const size_t PRISMATIC_PRINT_FREQUENCY_PROBES = max((size_t)1,pars.numProbes / 10); // for printing status
		WorkDispatcher dispatcher(0, pars.numProbes);

		for (auto t = 0; t < total_num_streams; ++t) {

			int GPU_num = stream_count % pars.meta.numGPUs; // determine which GPU handles this job

			cudaStream_t &current_stream = cuda_pars.streams[stream_count];
			cout << "Launching GPU worker on stream #" << stream_count << " of GPU #" << GPU_num << endl;

			// get pointers to the pre-copied arrays, making sure to get those on the current GPU
			PRISMATIC_CUDA_COMPLEX_FLOAT *current_PsiProbeInit_d       = cuda_pars.PsiProbeInit_d[GPU_num];
			PRISMATIC_FLOAT_PRECISION *current_qxaReduce_d             = cuda_pars.qxaReduce_d[GPU_num];
			PRISMATIC_FLOAT_PRECISION *current_qyaReduce_d             = cuda_pars.qyaReduce_d[GPU_num];
			size_t *current_yBeams_d                               = cuda_pars.yBeams_d[GPU_num];
			size_t *current_xBeams_d                               = cuda_pars.xBeams_d[GPU_num];
			PRISMATIC_FLOAT_PRECISION *current_alphaInd_d              = cuda_pars.alphaInd_d[GPU_num];

			// get pointers to per-stream arrays
			PRISMATIC_CUDA_COMPLEX_FLOAT *current_permutedScompact_ds = cuda_pars.permutedScompact_d[stream_count];
			PRISMATIC_CUDA_COMPLEX_FLOAT *current_psi_ds               = cuda_pars.psi_ds[stream_count];
			PRISMATIC_CUDA_COMPLEX_FLOAT *current_phaseCoeffs_ds       = cuda_pars.phaseCoeffs_ds[stream_count];
			long *current_y_ds                                     = cuda_pars.y_ds[stream_count];
			long *current_x_ds                                     = cuda_pars.x_ds[stream_count];
			cufftHandle &current_cufft_plan                        = cuda_pars.cufftPlans[stream_count];

			// get pointer to output pinned memory
			PRISMATIC_FLOAT_PRECISION *current_psiIntensity_ds = cuda_pars.psiIntensity_ds[stream_count];
			PRISMATIC_FLOAT_PRECISION *current_integratedOutput_ds = cuda_pars.integratedOutput_ds[stream_count];
			PRISMATIC_FLOAT_PRECISION *current_output_ph = cuda_pars.output_ph[stream_count];
			
			workers_GPU.push_back(thread([&pars, GPU_num, stream_count, current_permutedScompact_ds,  &dispatcher,
					                             current_alphaInd_d, current_PsiProbeInit_d, current_qxaReduce_d, current_qyaReduce_d,
					                             current_yBeams_d, current_xBeams_d, current_psi_ds, current_phaseCoeffs_ds,
					                             current_psiIntensity_ds, current_y_ds, current_x_ds, current_integratedOutput_ds, current_output_ph,
					                             &current_cufft_plan, &current_stream, &PRISMATIC_PRINT_FREQUENCY_PROBES, &cuda_pars]() {
				cudaErrchk(cudaSetDevice(GPU_num));

#ifndef NDEBUG
				{
//					 check memory usage on the GPU
					std::lock_guard<mutex> lock(Prismatic::memLock);
					size_t free_mem, total_mem;
					free_mem=total_mem=0;
					cudaErrchk(cudaMemGetInfo(&free_mem, &total_mem));
					pars.maxGPUMem = std::max(total_mem - free_mem, pars.maxGPUMem);
//					cout << "maxGPUMem = " << pars.maxGPUMem << endl;
				}
#endif // NDEBUG

				size_t Nstart, Nstop, ay, ax;
				Nstart=Nstop=0;
				while (dispatcher.getWork(Nstart, Nstop)) { // synchronously get work assignment
					while (Nstart < Nstop) {
						if (Nstart % PRISMATIC_PRINT_FREQUENCY_PROBES == 0 | Nstart == 100){
							cout << "Computing Probe Position #" << Nstart << "/" << pars.numProbes << endl;
						}
						ay = (pars.meta.arbitraryProbes) ? Nstart : Nstart / pars.numXprobes;
						ax = (pars.meta.arbitraryProbes) ? Nstart : Nstart % pars.numXprobes;
						buildSignal_GPU_streaming(pars, ay, ax, current_permutedScompact_ds, cuda_pars.permutedScompact_ph,
													current_PsiProbeInit_d, current_qxaReduce_d, current_qyaReduce_d,
													current_yBeams_d, current_xBeams_d, current_alphaInd_d, current_psi_ds,
													current_phaseCoeffs_ds, current_psiIntensity_ds, current_y_ds,
													current_x_ds, current_output_ph, current_integratedOutput_ds, current_cufft_plan, current_stream,  cuda_pars );
#ifdef PRISMATIC_BUILDING_GUI
						pars.progressbar->signalOutputUpdate(Nstart, pars.numProbes);
#endif
						++Nstart;
					}
				}
			}));
			++stream_count;
		}


		// Now launch CPU work
		if (pars.meta.alsoDoCPUWork) {
			PRISMATIC_FFTW_INIT_THREADS();
			PRISMATIC_FFTW_PLAN_WITH_NTHREADS(pars.meta.numThreads);
			vector<thread> workers_CPU;
			workers_CPU.reserve(pars.meta.numThreads); // prevents multiple reallocations
			for (auto t = 0; t < pars.meta.numThreads; ++t) {
				cout << "Launching CPU worker thread #" << t << " to compute partial PRISM result\n";
				workers_CPU.push_back(thread([&pars, &dispatcher, t, &PRISMATIC_PRINT_FREQUENCY_PROBES]() {
					size_t Nstart, Nstop, ay, ax, early_CPU_stop;
					Nstart=Nstop=0;
					if (pars.meta.numGPUs > 0){
						// if there are no GPUs, make sure to do all work on CPU
						early_CPU_stop = (size_t)std::max((PRISMATIC_FLOAT_PRECISION)0.0, pars.numProbes - pars.meta.earlyCPUStopCount);
					} else {
						early_CPU_stop = pars.numProbes;
					}
//					while (getWorkID(pars, Nstart, Nstop)) { // synchronously get work assignment
					if(dispatcher.getWork(Nstart, Nstop, 1, early_CPU_stop)) { // synchronously get work assignment
						Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> > psi = Prismatic::zeros_ND<2, std::complex<PRISMATIC_FLOAT_PRECISION> > (
								{{pars.imageSizeReduce[0], pars.imageSizeReduce[1]}});
						unique_lock<mutex> gatekeeper(fftw_plan_lock);


						PRISMATIC_FFTW_PLAN plan = PRISMATIC_FFTW_PLAN_DFT_2D(psi.get_dimj(), psi.get_dimi(),
						                                              reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&psi[0]),
						                                              reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&psi[0]),
						                                              FFTW_FORWARD, FFTW_MEASURE);
						gatekeeper.unlock();

						// main work loop
						do {
							while (Nstart < Nstop) {
								if (Nstart % PRISMATIC_PRINT_FREQUENCY_PROBES == 0 | Nstart == 100){
									cout << "Computing Probe Position #" << Nstart << "/" << pars.numProbes << endl;
								}
								ay = (pars.meta.arbitraryProbes) ? Nstart : Nstart / pars.numXprobes;
								ax = (pars.meta.arbitraryProbes) ? Nstart : Nstart % pars.numXprobes;
								buildSignal_CPU(pars, ay, ax, plan, psi);
#ifdef PRISMATIC_BUILDING_GUI
								pars.progressbar->signalOutputUpdate(Nstart, pars.numProbes);
#endif
								++Nstart;
							}
							if (Nstop >= early_CPU_stop) break;
						} while(dispatcher.getWork(Nstart, Nstop, 1, early_CPU_stop));
						gatekeeper.lock();
						PRISMATIC_FFTW_DESTROY_PLAN(plan);
						gatekeeper.unlock();
					}
				}));
			}
			cout << "Waiting for CPU threads...\n";
			for (auto& t:workers_CPU)t.join();
			PRISMATIC_FFTW_CLEANUP_THREADS();
		}

		// synchronize
		cout << "Waiting for GPU threads...\n";
		for (auto &t:workers_GPU)t.join();

		cout << "Synchronizing" << endl;
		for (auto g = 0; g < pars.meta.numGPUs; ++g){
			cudaSetDevice(g);
			cudaDeviceSynchronize();
		}
	}

	inline void cleanupMemory3(Parameters<PRISMATIC_FLOAT_PRECISION>& pars,
	                          CudaParameters<PRISMATIC_FLOAT_PRECISION>& cuda_pars){
		const int total_num_streams = pars.meta.numGPUs * pars.meta.numStreamsPerGPU;
		// free pinned memory
		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaFreeHost(cuda_pars.output_ph[s]));
		}
		cudaErrchk(cudaFreeHost(cuda_pars.permutedScompact_ph));
		cudaErrchk(cudaFreeHost(cuda_pars.PsiProbeInit_ph));
		cudaErrchk(cudaFreeHost(cuda_pars.qxaReduce_ph));
		cudaErrchk(cudaFreeHost(cuda_pars.qyaReduce_ph));
		cudaErrchk(cudaFreeHost(cuda_pars.xBeams_ph));
		cudaErrchk(cudaFreeHost(cuda_pars.yBeams_ph));

		for (auto g = 0; g < pars.meta.numGPUs; ++g) {
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaFree(cuda_pars.permutedScompact_d[g]));
			cudaErrchk(cudaFree(cuda_pars.PsiProbeInit_d[g]));
			cudaErrchk(cudaFree(cuda_pars.qxaReduce_d[g]));
			cudaErrchk(cudaFree(cuda_pars.qyaReduce_d[g]));
			cudaErrchk(cudaFree(cuda_pars.yBeams_d[g]));
			cudaErrchk(cudaFree(cuda_pars.xBeams_d[g]));
			cudaErrchk(cudaFree(cuda_pars.alphaInd_d[g]));
		}

		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaSetDevice(s % pars.meta.numGPUs));
			cudaErrchk(cudaFree(cuda_pars.psi_ds[s]));
			cudaErrchk(cudaFree(cuda_pars.phaseCoeffs_ds[s]));
			cudaErrchk(cudaFree(cuda_pars.y_ds[s]));
			cudaErrchk(cudaFree(cuda_pars.x_ds[s]));
			cufftErrchk(cufftDestroy(cuda_pars.cufftPlans[s]));
			cudaErrchk(cudaFree(cuda_pars.psiIntensity_ds[s]));
			cudaErrchk(cudaFree(cuda_pars.integratedOutput_ds[s]));
		}

		for (auto j = 0; j < pars.meta.numGPUs; ++j) {
			cudaErrchk(cudaSetDevice(j));
			cudaErrchk(cudaDeviceReset());
		}

		delete[] cuda_pars.streams;
		delete[] cuda_pars.cufftPlans;
		delete[] cuda_pars.PsiProbeInit_d;
		delete[] cuda_pars.qxaReduce_d;
		delete[] cuda_pars.qyaReduce_d;
		delete[] cuda_pars.alphaInd_d;
		delete[] cuda_pars.yBeams_d;
		delete[] cuda_pars.xBeams_d;
		delete[] cuda_pars.permutedScompact_d;
		delete[] cuda_pars.psi_ds;
		delete[] cuda_pars.phaseCoeffs_ds;
		delete[] cuda_pars.y_ds;
		delete[] cuda_pars.x_ds;
		delete[] cuda_pars.output_ph;
		delete[] cuda_pars.psiIntensity_ds;
		delete[] cuda_pars.integratedOutput_ds;
	}


	//There are a number of template specializations that follow. They exploit the fact that scaleReduceS provides one
	// element of shared memory per thread. Therefore, for every reduction operation other than the one involving BlockSize_numBeams/2 threads
	// the bounds checking operation can be skipped -- the latter threads just operate on meaningless values
	template <size_t BlockSize_numBeams>
	__device__  void warpReduce_cx(volatile cuFloatComplex* sdata, int idx){
		// When 32 or fewer threads remain, there is only a single warp remaining and no need to synchronize; however,
		// the volatile keyword is necessary otherwise the compiler will optimize these operations into registers
		// and the result will be incorrect
		if (BlockSize_numBeams >= 64){
			sdata[idx].x += sdata[idx + 32].x;
			sdata[idx].y += sdata[idx + 32].y;
		}
		if (BlockSize_numBeams >= 32){
			sdata[idx].x += sdata[idx + 16].x;
			sdata[idx].y += sdata[idx + 16].y;
		}
		if (BlockSize_numBeams >= 16){
			sdata[idx].x += sdata[idx + 8].x;
			sdata[idx].y += sdata[idx + 8].y;
		}
		if (BlockSize_numBeams >= 8){
			sdata[idx].x += sdata[idx + 4].x;
			sdata[idx].y += sdata[idx + 4].y;
		}
		if (BlockSize_numBeams >= 4){
			sdata[idx].x += sdata[idx + 2].x;
			sdata[idx].y += sdata[idx + 2].y;
		}
		if (BlockSize_numBeams >= 2){
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

	template <>
	__device__  void warpReduce_cx<1>(volatile cuFloatComplex* sdata, int idx){
	}



	template <size_t BlockSize_numBeams>
	__device__  void warpReduce_cx(volatile cuDoubleComplex* sdata, int idx){
		// When 32 or fewer threads remain, there is only a single warp remaining and no need to synchronize; however,
		// the volatile keyword is necessary otherwise the compiler will optimize these operations into registers
		// and the result will be incorrect
		if (BlockSize_numBeams >= 64){
			sdata[idx].x += sdata[idx + 32].x;
			sdata[idx].y += sdata[idx + 32].y;
		}
		if (BlockSize_numBeams >= 32){
			sdata[idx].x += sdata[idx + 16].x;
			sdata[idx].y += sdata[idx + 16].y;
		}
		if (BlockSize_numBeams >= 16){
			sdata[idx].x += sdata[idx + 8].x;
			sdata[idx].y += sdata[idx + 8].y;
		}
		if (BlockSize_numBeams >= 8){
			sdata[idx].x += sdata[idx + 4].x;
			sdata[idx].y += sdata[idx + 4].y;
		}
		if (BlockSize_numBeams >= 4){
			sdata[idx].x += sdata[idx + 2].x;
			sdata[idx].y += sdata[idx + 2].y;
		}
		if (BlockSize_numBeams >= 2){
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
	template <>
	__device__  void warpReduce_cx<1>(volatile cuDoubleComplex* sdata, int idx){
		// When 32 or fewer threads remain, there is only a single warp remaining and no need to synchronize; however,
		// the volatile keyword is necessary otherwise the compiler will optimize these operations into registers
		// and the result will be incorrect
	}

	template <size_t BlockSizeX>
	__global__ void scaleReduceS(const cuFloatComplex *permutedScompact_d,
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
		int y   = threadIdx.y + blockDim.y * blockIdx.y;
		int z   = threadIdx.z + blockDim.z * blockIdx.z;

		// determine grid size for stepping through the array
		int gridSizeY = gridDim.y * blockDim.y;
		int gridSizeZ = gridDim.z * blockDim.z;

		// guarantee the shared memory is initialized to 0 so we can accumulate without bounds checking
		scaled_values[threadIdx.x] = make_cuFloatComplex(0,0);
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

//		// each block processes several reductions strided by the grid size
		int y_saved = y;
		while (z < dimj_psi){
			y=y_saved; // reset y
			while(y < dimi_psi){
				//	 read in first values
				if (idx < numberBeams) {
					scaled_values[idx] = cuCmulf(permutedScompact_d[z_ds[z]*numberBeams*dimj_S + y_ds[y]*numberBeams + idx],
					                             coeff_cache[idx]);
					__syncthreads();
				}

//		 step through global memory accumulating until values have been reduced to BlockSizeX elements in shared memory
				size_t offset = BlockSizeX;
				while (offset < numberBeams){
					if (idx + offset < numberBeams){
						scaled_values[idx] = cuCaddf(scaled_values[idx],
						                             cuCmulf( permutedScompact_d[z_ds[z]*numberBeams*dimj_S + y_ds[y]*numberBeams + idx + offset],
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
				if (idx < 32 & BlockSizeX <= numberBeams){
					warpReduce_cx<BlockSizeX>(scaled_values, idx);

				} else {
					warpReduce_cx<1>(scaled_values, idx);
				}

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

	template <size_t BlockSize_numBeams, size_t BlockSize_alongArray>
__global__ void scaleReduceS(const cuFloatComplex *permutedScompact_d,
                             const cuFloatComplex *phaseCoeffs_ds,
                             cuFloatComplex *psi_ds,
                             const long *z_ds,
                             const long* y_ds,
                             const size_t numberBeams,
                             const size_t dimk_S,
                             const size_t dimj_S,
                             const size_t dimj_psi,
                             const size_t dimi_psi){
		// This code is heavily modeled after Mark Harris's presentation on optimized parallel reduction
		// http://developer.download.nvidia.com/compute/cuda/1.1-Beta/x86_website/projects/reduction/doc/reduction.pdf

		// shared memory
		__shared__ cuFloatComplex scaled_values[BlockSize_numBeams * BlockSize_alongArray]; // for holding the values to reduce
		extern __shared__ cuFloatComplex coeff_cache []; // cache the coefficients to prevent repeated global reads

		// for the permuted Scompact matrix, the x direction runs along the number of beams, leaving y and z to represent the
		//	2D array of reduced values in psi

		int beam_idx        = threadIdx.x; // index along the beam direction, there is only ever one block along this dimension
		int array_idx       = threadIdx.y + blockDim.y * blockIdx.y; // index along the (flattened) array direction
		int t_id            = threadIdx.x + BlockSize_numBeams*threadIdx.y; // linear index of the thread within the block
		size_t array_offset = threadIdx.y * BlockSize_numBeams; // offset of where each block of shared memory begins for each reduction job
//		size_t t_id = array_offset + beam_idx;
		// determine grid size for stepping through the array

		int gridSize_alongArray = gridDim.y * blockDim.y;

		// guarantee the shared memory is initialized to 0 so we can accumulate without bounds checking
		scaled_values[t_id] = make_cuFloatComplex(0, 0);
		__syncthreads();

		// read the coefficients into shared memory once
		size_t offset_phase_idx = 0;
		const size_t inc = BlockSize_numBeams*BlockSize_alongArray;
		while (offset_phase_idx < numberBeams){
			if (t_id  < numberBeams){
				coeff_cache[t_id] = phaseCoeffs_ds[t_id + offset_phase_idx];
			}
			offset_phase_idx += inc;
		}
		__syncthreads();

//		// each block processes several reductions strided by the grid size
		while (array_idx < dimj_psi * dimi_psi){
			int y   = array_idx % dimi_psi;
			int z   = array_idx / dimi_psi;

	//	 read in first values
		if (beam_idx < numberBeams) {
			scaled_values[t_id] = cuCmulf(permutedScompact_d[z_ds[z]*numberBeams*dimj_S + y_ds[y]*numberBeams + beam_idx],
			                                                 coeff_cache[beam_idx]);
			__syncthreads();
		}

//		 step through global memory accumulating until values have been reduced to BlockSize_numBeams elements in shared memory
		size_t offset = BlockSize_numBeams;
		while (offset < numberBeams){
			if (beam_idx + offset < numberBeams){
				scaled_values[t_id] = cuCaddf(scaled_values[t_id],
				                                          cuCmulf(permutedScompact_d[z_ds[z]*numberBeams*dimj_S + y_ds[y]*numberBeams + beam_idx + offset],
				                                                  coeff_cache[beam_idx + offset]));
			}
			offset += BlockSize_numBeams;
			__syncthreads();
		}

		// At this point we have exactly BlockSize_numBeams elements to reduce from shared memory which we will add by recursively
		// dividing the array in half

		// Take advantage of templates. Because BlockSize_numBeams is passed at compile time, all of these comparisons are also
		// evaluated at compile time
		if (BlockSize_numBeams >= 1024){
			if (beam_idx < 512){
				scaled_values[t_id] = cuCaddf(scaled_values[t_id], scaled_values[t_id + 512]);
			}
			__syncthreads();
		}

		if (BlockSize_numBeams >= 512){
			if (beam_idx < 256){
				scaled_values[t_id] = cuCaddf(scaled_values[t_id], scaled_values[t_id + 256]);
			}
			__syncthreads();
		}

		if (BlockSize_numBeams >= 256){
			if (beam_idx < 128){
				scaled_values[t_id] = cuCaddf(scaled_values[t_id], scaled_values[t_id + 128]);
			}
			__syncthreads();
		}

		if (BlockSize_numBeams >= 128){
			if (beam_idx < 64){
				scaled_values[t_id] = cuCaddf(scaled_values[t_id], scaled_values[t_id + 64]);
			}
			__syncthreads();
		}
			if (BlockSize_numBeams >= 64){
				if (beam_idx < 32){
					scaled_values[t_id] = cuCaddf(scaled_values[t_id], scaled_values[t_id + 32]);
				}
				__syncthreads();
			}

			if (BlockSize_numBeams >= 32){
				if (beam_idx < 16){
					scaled_values[t_id] = cuCaddf(scaled_values[t_id], scaled_values[t_id + 16]);
				}
				__syncthreads();
			}

			if (BlockSize_numBeams >= 16){
				if (beam_idx < 8){
					scaled_values[t_id] = cuCaddf(scaled_values[t_id], scaled_values[t_id + 8]);
				}
				__syncthreads();
			}

			if (BlockSize_numBeams >= 8){
				if (beam_idx < 4){
					scaled_values[t_id] = cuCaddf(scaled_values[t_id], scaled_values[t_id + 4]);
				}
				__syncthreads();
			}

			if (BlockSize_numBeams >= 4){
				if (beam_idx < 2){
					scaled_values[t_id] = cuCaddf(scaled_values[t_id], scaled_values[t_id + 2]);
				}
				__syncthreads();
			}
			if (BlockSize_numBeams >= 2){
				if (beam_idx < 1){
					scaled_values[t_id] = cuCaddf(scaled_values[t_id], scaled_values[t_id + 1]);
				}
				__syncthreads();
			}
		// write out the result
		if (beam_idx == 0)psi_ds[z*dimi_psi + y] = scaled_values[array_offset];
		// increment
		array_idx+=gridSize_alongArray;
			__syncthreads();
	}
}

	template <size_t BlockSize_numBeams, size_t BlockSize_alongArray>
	__global__ void scaleReduceS(const cuDoubleComplex *permutedScompact_d,
	                             const cuDoubleComplex *phaseCoeffs_ds,
	                             cuDoubleComplex *psi_ds,
	                             const long *z_ds,
	                             const long* y_ds,
	                             const size_t numberBeams,
	                             const size_t dimk_S,
	                             const size_t dimj_S,
	                             const size_t dimj_psi,
	                             const size_t dimi_psi){
		// This code is heavily modeled after Mark Harris's presentation on optimized parallel reduction
		// http://developer.download.nvidia.com/compute/cuda/1.1-Beta/x86_website/projects/reduction/doc/reduction.pdf

		// shared memory
		__shared__ cuDoubleComplex scaled_values[BlockSize_numBeams * BlockSize_alongArray]; // for holding the values to reduce
		extern __shared__ cuDoubleComplex coeff_cache_double[]; // cache the coefficients to prevent repeated global reads

		// for the permuted Scompact matrix, the x direction runs along the number of beams, leaving y and z to represent the
		//	2D array of reduced values in psi

		int beam_idx    = threadIdx.x; // index along the beam direction, there is only ever one block along this dimension
		int array_idx   = threadIdx.y + blockDim.y * blockIdx.y; // index along the (flattened) array direction
		int t_id        = threadIdx.x + BlockSize_numBeams*threadIdx.y; // linear index of the thread within the block
		size_t array_offset = threadIdx.y * BlockSize_numBeams; // offset of where each block of shared memory begins for each reduction job
		// determine grid size for stepping through the array

		int gridSize_alongArray = gridDim.y * blockDim.y;

		// guarantee the shared memory is initialized to 0 so we can accumulate without bounds checking
		scaled_values[t_id] = make_cuDoubleComplex(0,0);
		__syncthreads();

		// read the coefficients into shared memory once
		size_t offset_phase_idx = 0;
		const size_t inc = BlockSize_numBeams * BlockSize_alongArray;
		while (offset_phase_idx < numberBeams){
			if (t_id < numberBeams){
				coeff_cache_double[t_id] = phaseCoeffs_ds[t_id + offset_phase_idx];
			}
			offset_phase_idx += inc;
		}
		__syncthreads();

//		// each block processes several reductions strided by the grid size
		while (array_idx < dimj_psi * dimi_psi){
			int y   = array_idx % dimi_psi;
			int z   = array_idx / dimi_psi;

			//	 read in first values
			if (beam_idx < numberBeams) {
				scaled_values[t_id] = cuCmul(permutedScompact_d[z_ds[z]*numberBeams*dimj_S + y_ds[y]*numberBeams + beam_idx],
				                              coeff_cache_double[beam_idx]);
				__syncthreads();
			}

//		 step through global memory accumulating until values have been reduced to BlockSize_numBeams elements in shared memory
			size_t offset = BlockSize_numBeams;
			while (offset < numberBeams){
				if (beam_idx + offset < numberBeams){
					scaled_values[t_id] = cuCadd(scaled_values[t_id],
					                              cuCmul(permutedScompact_d[z_ds[z]*numberBeams*dimj_S + y_ds[y]*numberBeams + beam_idx + offset],
					                                      coeff_cache_double[beam_idx + offset]));
				}
				offset += BlockSize_numBeams;
				__syncthreads();
			}

			// At this point we have exactly BlockSize_numBeams elements to reduce from shared memory which we will add by recursively
			// dividing the array in half

			// Take advantage of templates. Because BlockSize_numBeams is passed at compile time, all of these comparisons are also
			// evaluated at compile time
			if (BlockSize_numBeams >= 1024){
				if (beam_idx < 512){
					scaled_values[t_id] = cuCadd(scaled_values[t_id], scaled_values[t_id + 512]);
				}
				__syncthreads();
			}

			if (BlockSize_numBeams >= 512){
				if (beam_idx < 256){
					scaled_values[t_id] = cuCadd(scaled_values[t_id], scaled_values[t_id + 256]);
				}
				__syncthreads();
			}

			if (BlockSize_numBeams >= 256){
				if (beam_idx < 128){
					scaled_values[t_id] = cuCadd(scaled_values[t_id], scaled_values[t_id + 128]);
				}
				__syncthreads();
			}

			if (BlockSize_numBeams >= 128){
				if (beam_idx < 64){
					scaled_values[t_id] = cuCadd(scaled_values[t_id], scaled_values[t_id + 64]);
				}
				__syncthreads();
			}
			if (BlockSize_numBeams >= 64){
				if (beam_idx < 32){
					scaled_values[t_id] = cuCadd(scaled_values[t_id], scaled_values[t_id + 32]);
				}
				__syncthreads();
			}

			if (BlockSize_numBeams >= 32){
				if (beam_idx < 16){
					scaled_values[t_id] = cuCadd(scaled_values[t_id], scaled_values[t_id + 16]);
				}
				__syncthreads();
			}

			if (BlockSize_numBeams >= 16){
				if (beam_idx < 8){
					scaled_values[t_id] = cuCadd(scaled_values[t_id], scaled_values[t_id + 8]);
				}
				__syncthreads();
			}

			if (BlockSize_numBeams >= 8){
				if (beam_idx < 4){
					scaled_values[t_id] = cuCadd(scaled_values[t_id], scaled_values[t_id + 4]);
				}
				__syncthreads();
			}

			if (BlockSize_numBeams >= 4){
				if (beam_idx < 2){
					scaled_values[t_id] = cuCadd(scaled_values[t_id], scaled_values[t_id + 2]);
				}
				__syncthreads();
			}
			if (BlockSize_numBeams >= 2){
				if (beam_idx < 1){
					scaled_values[t_id] = cuCadd(scaled_values[t_id], scaled_values[t_id + 1]);
				}
				__syncthreads();
			}
			// write out the result
			if (beam_idx == 0)psi_ds[z*dimi_psi + y] = scaled_values[array_offset];
			// increment
			array_idx+=gridSize_alongArray;
			__syncthreads();
		}
	}

	template <size_t BlockSize_numBeams>
	// double precision version, see float version above for comments
	__global__ void scaleReduceS(const cuDoubleComplex *permutedScompact_d,
	                             const cuDoubleComplex *phaseCoeffs_ds,
	                             cuDoubleComplex *psi_ds,
	                             const long *z_ds,
	                             const long* y_ds,
	                             const size_t numberBeams,
	                             const size_t dimk_S,
	                             const size_t dimj_S,
	                             const size_t dimj_psi,
	                             const size_t dimi_psi){
		__shared__ cuDoubleComplex scaled_values[BlockSize_numBeams];
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
			offset_phase_idx += BlockSize_numBeams;
		}
		__syncthreads();
		int y_saved = y;
		while (z < dimj_psi){
			y=y_saved;
			while(y < dimi_psi){
				if (idx < numberBeams) {
					scaled_values[idx] = cuCmul(permutedScompact_d[z_ds[z]*numberBeams*dimj_S + y_ds[y]*numberBeams + idx],
					                             coeff_cache_double[idx]);
					__syncthreads();
				}
				size_t offset = BlockSize_numBeams;
				while (offset < numberBeams){
					if (idx + offset < numberBeams){
						scaled_values[idx] = cuCadd(scaled_values[idx],
						                             cuCmul( permutedScompact_d[z_ds[z]*numberBeams*dimj_S + y_ds[y]*numberBeams + idx + offset],
						                                      coeff_cache_double[idx + offset]));
					}
					offset += BlockSize_numBeams;
					__syncthreads();
				}
				if (BlockSize_numBeams >= 1024){
					if (idx < 512){
						scaled_values[idx] = cuCadd(scaled_values[idx], scaled_values[idx + 512]);
					}
					__syncthreads();
				}
				if (BlockSize_numBeams >= 512){
					if (idx < 256){
						scaled_values[idx] = cuCadd(scaled_values[idx], scaled_values[idx + 256]);
					}
					__syncthreads();
				}
				if (BlockSize_numBeams >= 256){
					if (idx < 128){
						scaled_values[idx] = cuCadd(scaled_values[idx], scaled_values[idx + 128]);
					}
					__syncthreads();
				}
				if (BlockSize_numBeams >= 128){
					if (idx < 64){
						scaled_values[idx] = cuCadd(scaled_values[idx], scaled_values[idx + 64]);
					}
					__syncthreads();
				}
				if (idx < 32)warpReduce_cx<BlockSize_numBeams>(scaled_values, idx);
				if (idx == 0)psi_ds[z*dimi_psi + y] = scaled_values[0];
				y+=gridSizeY;
				__syncthreads();
			}
			z+=gridSizeZ;
			__syncthreads();
		}
	}

	using namespace std;

	void buildPRISMOutput_GPU_singlexfer(Parameters<PRISMATIC_FLOAT_PRECISION> &pars){
#ifdef PRISMATIC_BUILDING_GUI
		pars.progressbar->signalDescriptionMessage("Computing final output (PRISM)");
#endif
		// construct the PRISM output array using GPUs

		CudaParameters<PRISMATIC_FLOAT_PRECISION> cuda_pars;

		// create CUDA streams and cuFFT plans
		createStreamsAndPlans3(pars, cuda_pars);
		// create page-locked (pinned) host memory buffers
		allocatePinnedHostMemory_singlexfer3(pars, cuda_pars);
		
		// copy data to pinned buffers
		copyToPinnedMemory_singlexfer3(pars, cuda_pars);
		
		// allocate memory on the GPUs
		allocateDeviceMemory_singlexfer3(pars, cuda_pars);
		
		// copy memory to GPUs
		copyToGPUMemory_singlexfer3(pars, cuda_pars);
		
		// launch GPU and CPU workers
		launchWorkers_singlexfer3(pars, cuda_pars);
		
		// free memory on the host/device
		cleanupMemory3(pars, cuda_pars);
	}

	void buildPRISMOutput_GPU_streaming(Parameters<PRISMATIC_FLOAT_PRECISION> &pars){
#ifdef PRISMATIC_BUILDING_GUI
		pars.progressbar->signalDescriptionMessage("Computing final output (PRISM)");
#endif
		CudaParameters<PRISMATIC_FLOAT_PRECISION> cuda_pars;
		// construct the PRISM output array using GPUs

		// create CUDA streams and cuFFT plans
		createStreamsAndPlans3(pars, cuda_pars);

		// allocate pinned memory
		allocatePinnedHostMemory_streaming3(pars, cuda_pars);

		// copy data to pinned buffers
		copyToPinnedMemory_streaming3(pars, cuda_pars);

		// allocate memory on the GPUs
		allocateDeviceMemory_streaming3(pars, cuda_pars);

		// copy memory to GPUs
		copyToGPUMemory_streaming3(pars, cuda_pars);

		// launch GPU and CPU workers
		launchWorkers_streaming3(pars, cuda_pars);

		// free memory on the host/device
		cleanupMemory3(pars, cuda_pars);
	}

	void buildSignal_GPU_singlexfer(Parameters<PRISMATIC_FLOAT_PRECISION>&  pars,
	                                const size_t& ay,
	                                const size_t& ax,
	                                const PRISMATIC_CUDA_COMPLEX_FLOAT *permutedScompact_d,
	                                const PRISMATIC_CUDA_COMPLEX_FLOAT *PsiProbeInit_d,
	                                const PRISMATIC_FLOAT_PRECISION *qxaReduce_d,
	                                const PRISMATIC_FLOAT_PRECISION *qyaReduce_d,
	                                const size_t *yBeams_d,
	                                const size_t *xBeams_d,
	                                const PRISMATIC_FLOAT_PRECISION *alphaInd_d,
	                                PRISMATIC_CUDA_COMPLEX_FLOAT *psi_ds,
	                                PRISMATIC_CUDA_COMPLEX_FLOAT *phaseCoeffs_ds,
	                                PRISMATIC_FLOAT_PRECISION *psiIntensity_ds,
	                                long  *y_ds,
	                                long  *x_ds,
	                                PRISMATIC_FLOAT_PRECISION *output_ph,
	                                PRISMATIC_FLOAT_PRECISION *integratedOutput_ds,
	                                const cufftHandle &cufft_plan,
	                                const cudaStream_t& stream,
	                                CudaParameters<PRISMATIC_FLOAT_PRECISION>& cuda_pars){

		const PRISMATIC_FLOAT_PRECISION yp = pars.yp[ay];
		const PRISMATIC_FLOAT_PRECISION xp = pars.xp[ax];

		const size_t psi_size = pars.imageSizeReduce[0] * pars.imageSizeReduce[1];
		shiftIndices <<<((long)pars.imageSizeReduce[0] - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream>>> (
				y_ds, (long)std::round(yp / (PRISMATIC_FLOAT_PRECISION)pars.pixelSizeOutput[0]), (long)pars.imageSizeOutput[0], (long)pars.imageSizeReduce[0]);

		shiftIndices <<<((long)pars.imageSizeReduce[1] - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream>>> (
				x_ds, (long)std::round(xp / (PRISMATIC_FLOAT_PRECISION)pars.pixelSizeOutput[1]), (long)pars.imageSizeOutput[1], (long)pars.imageSizeReduce[1]);

		computePhaseCoeffs <<<(pars.numberBeams - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream>>>(
                   phaseCoeffs_ds, PsiProbeInit_d, qyaReduce_d, qxaReduce_d,
		           yBeams_d, xBeams_d, yp, xp, pars.yTiltShift, pars.xTiltShift, pars.imageSizeReduce[1], pars.numberBeams);

		// Choose a good launch configuration
		// Heuristically use 2^p / 2 as the block size where p is the first power of 2 greater than the number of elements to work on.
		// This balances having enough work per thread and enough blocks without having so many blocks that the shared memory doesn't last long
		size_t p = getNextPower2(pars.numberBeams);
		const size_t BlockSize_numBeams = min((size_t)pars.deviceProperties.maxThreadsPerBlock,(size_t)std::max(1.0, pow(2,p) / 2));

		// Determine maximum threads per streaming multiprocessor based on the compute capability of the device
		size_t max_threads_per_sm;
		if (pars.deviceProperties.major > 3){
			max_threads_per_sm = 2048;
		} else if (pars.deviceProperties.major > 2) {
			max_threads_per_sm = 1536;
		} else {
			max_threads_per_sm = pars.deviceProperties.minor == 0 ? 768 : 1024;
		}

		// Estimate max number of blocks per streaming multiprocessor
		const size_t max_blocks_per_sm = std::min((size_t)32, max_threads_per_sm / BlockSize_numBeams);

		// We find providing around 3 times as many blocks as the estimated maximum provides good performance
		const size_t target_blocks_per_sm = max_blocks_per_sm * 3;
		const size_t total_blocks         = target_blocks_per_sm * pars.deviceProperties.multiProcessorCount;

		// Determine the shape of the grid
		
		if (BlockSize_numBeams >= 64){
			const PRISMATIC_FLOAT_PRECISION aspect_ratio = (PRISMATIC_FLOAT_PRECISION)pars.imageSizeReduce[1] / (PRISMATIC_FLOAT_PRECISION)pars.imageSizeReduce[0];
			const size_t GridSizeZ = std::floor(sqrt(total_blocks / aspect_ratio));
			const size_t GridSizeY = aspect_ratio * GridSizeZ;
			dim3 grid(1, GridSizeY, GridSizeZ);
			dim3 block(BlockSize_numBeams, 1, 1);
			// Determine amount of shared memory needed
			const unsigned long smem = pars.numberBeams * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT);

			// Launch kernel. Block size must be visible at compile time so we use a switch statement
			switch (BlockSize_numBeams) {
				case 1024 :
					scaleReduceS<1024> <<< grid, block, smem, stream >>> (
							permutedScompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
							pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
				case 512 :
					scaleReduceS<512> <<< grid, block, smem, stream >>> (
							permutedScompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
							pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
				case 256 :
					scaleReduceS<256> <<< grid, block, smem, stream >>> (
							permutedScompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
							pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
				case 128 :
					scaleReduceS<128> <<< grid, block, smem, stream >>> (
							permutedScompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
							pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
				case 64 :
					scaleReduceS<64> <<< grid, block, smem, stream >>> (
							permutedScompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
							pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
				case 32 :
					scaleReduceS<32> <<< grid, block, smem, stream >>> (
							permutedScompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
							pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
				case 16 :
					scaleReduceS<16> <<< grid, block, smem, stream >>> (
							permutedScompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
							pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
				case 8 :
					scaleReduceS<8> <<< grid, block, smem, stream >>> (
							permutedScompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
							pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
				case 4 :
					scaleReduceS<4> <<< grid, block, smem, stream >>> (
							permutedScompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
							pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
				default :
					scaleReduceS<2> <<< grid, block, smem, stream >>> (
							permutedScompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
							pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]); break;
			}
		} else {
			const size_t BlockSize_alongArray = 512 / BlockSize_numBeams;
			const size_t GridSize_alongArray = pars.deviceProperties.multiProcessorCount * max_blocks_per_sm * 3;
			dim3 grid(1, GridSize_alongArray, 1);
			dim3 block(BlockSize_numBeams, BlockSize_alongArray, 1);
			// Determine amount of shared memory needed
			const unsigned long smem = pars.numberBeams * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT);
			// Launch kernel. Block size must be visible at compile time so we use a switch statement
			switch (BlockSize_numBeams) {
				case 1024 :
					scaleReduceS<1024, 1> << < grid, block, smem, stream >> > (
							permutedScompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
							pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
				case 512 :
					scaleReduceS<512, 1> << < grid, block, smem, stream >> > (
							permutedScompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
							pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
				case 256 :
					scaleReduceS<256, 2> << < grid, block, smem, stream >> > (
							permutedScompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
							pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
				case 128 :
					scaleReduceS<128, 4> << < grid, block, smem, stream >> > (
							permutedScompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
							pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
				case 64 :
					scaleReduceS<64, 8> << < grid, block, smem, stream >> > (
							permutedScompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
							pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
				case 32 :
					scaleReduceS<32, 16> << < grid, block, smem, stream >> > (
							permutedScompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
							pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
				case 16 :
					scaleReduceS<16, 32> << < grid, block, smem, stream >> > (
							permutedScompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
							pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
				case 8 :
					scaleReduceS<8, 64> << < grid, block, smem, stream >> > (
							permutedScompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
							pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
				case 4 :
					scaleReduceS<4, 128> << < grid, block, smem, stream >> > (
							permutedScompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
							pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
				default :
					scaleReduceS<1, 512> << < grid, block, smem, stream >> > (
							permutedScompact_d, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.Scompact.get_dimj(),
							pars.Scompact.get_dimi(), pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
			}
		}

		// final fft
		cufftErrchk(PRISMATIC_CUFFT_EXECUTE(cufft_plan, &psi_ds[0], &psi_ds[0], CUFFT_FORWARD));

		// convert to squared intensity
		abs_squared <<< (psi_size - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >>> (psiIntensity_ds, psi_ds, psi_size);

		// output calculation result
		size_t write_ay = (pars.meta.arbitraryProbes) ? 0 : ay;
		if(pars.meta.saveComplexOutputWave)
		{
			formatOutput_GPU_c_integrate(pars, psi_ds, psiIntensity_ds, alphaInd_d, output_ph,
				integratedOutput_ds, qyaReduce_d, qxaReduce_d, 0, write_ay, ax, pars.imageSizeReduce[0],
				pars.imageSizeReduce[1], stream, pars.scale);
		}
		else
		{
			formatOutput_GPU_integrate(pars, psiIntensity_ds, alphaInd_d, output_ph,
				integratedOutput_ds, qyaReduce_d, qxaReduce_d, 0, write_ay, ax, pars.imageSizeReduce[0],
				pars.imageSizeReduce[1], stream, pars.scale);
		}
	}

	void buildSignal_GPU_streaming(Parameters<PRISMATIC_FLOAT_PRECISION>&  pars,
	                               const size_t& ay,
	                               const size_t& ax,
	                               PRISMATIC_CUDA_COMPLEX_FLOAT *permutedScompact_ds,
	                               const std::complex<PRISMATIC_FLOAT_PRECISION> *permutedScompact_ph,
	                               const PRISMATIC_CUDA_COMPLEX_FLOAT *PsiProbeInit_d,
	                               const PRISMATIC_FLOAT_PRECISION *qxaReduce_d,
	                               const PRISMATIC_FLOAT_PRECISION *qyaReduce_d,
	                               const size_t *yBeams_d,
	                               const size_t *xBeams_d,
	                               const PRISMATIC_FLOAT_PRECISION *alphaInd_d,
	                               PRISMATIC_CUDA_COMPLEX_FLOAT *psi_ds,
	                               PRISMATIC_CUDA_COMPLEX_FLOAT *phaseCoeffs_ds,
	                               PRISMATIC_FLOAT_PRECISION *psiIntensity_ds,
	                               long  *y_ds,
	                               long  *x_ds,
	                               PRISMATIC_FLOAT_PRECISION *output_ph,
	                               PRISMATIC_FLOAT_PRECISION *integratedOutput_ds,
	                               const cufftHandle &cufft_plan,
	                               const cudaStream_t& stream,
	                               CudaParameters<PRISMATIC_FLOAT_PRECISION>& cuda_pars){

		// the coordinates y and x of the output image phi map to z and y of the permuted S compact matrix
		const PRISMATIC_FLOAT_PRECISION yp = pars.yp[ay];
		const PRISMATIC_FLOAT_PRECISION xp = pars.xp[ax];
		const size_t psi_size = pars.imageSizeReduce[0] * pars.imageSizeReduce[1];

		shiftIndices <<<(pars.imageSizeReduce[0] - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream>>> (
				y_ds, std::round(yp / pars.pixelSizeOutput[0]),pars.imageSize[0], pars.imageSizeReduce[0]);

		shiftIndices <<<(pars.imageSizeReduce[1] - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream>>> (
				x_ds, std::round(xp / pars.pixelSizeOutput[1]), pars.imageSize[1], pars.imageSizeReduce[1]);

		computePhaseCoeffs <<<(pars.numberBeams - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream>>>(
				phaseCoeffs_ds, PsiProbeInit_d, qyaReduce_d, qxaReduce_d,
				yBeams_d, xBeams_d, yp, xp, pars.yTiltShift, pars.xTiltShift, pars.imageSizeReduce[1], pars.numberBeams);

		// Copy the relevant portion of the Scompact matrix. This can be accomplished with ideally one but at most 4 strided 3-D memory copies
		// depending on whether or not the coordinates wrap around.
		long x1,y1;
		y1 = pars.yVec[0] +  std::round(yp / (PRISMATIC_FLOAT_PRECISION)pars.pixelSizeOutput[0]);
		x1 = pars.xVec[0] +  std::round(xp / (PRISMATIC_FLOAT_PRECISION)pars.pixelSizeOutput[1]);


		// determine where in the coordinate list wrap-around occurs (if at all)
		long xsplit, ysplit, nx2, ny2, xstart1, xstart2, ystart1, ystart2;
		xsplit = (x1 < 0) ? -x1 : (x1 + pars.xVec.size() > pars.Scompact.get_dimi()) ? pars.Scompact.get_dimi() - x1 : pars.xVec.size();
		ysplit = (y1 < 0) ? -y1 : (y1 + pars.yVec.size() > pars.Scompact.get_dimj()) ? pars.Scompact.get_dimj() - y1 : pars.yVec.size();

		nx2 = pars.xVec.size() - xsplit;
		ny2 = pars.yVec.size() - ysplit;

		xstart1 = ((long) pars.imageSizeOutput[1] + (x1 % (long) pars.imageSizeOutput[1])) %
		           (long) pars.imageSizeOutput[1];
		xstart2 = ((long) pars.imageSizeOutput[1] + (x1 + xsplit % (long) pars.imageSizeOutput[1])) %
		           (long) pars.imageSizeOutput[1];
		ystart1 = ((long) pars.imageSizeOutput[0] + (y1 % (long) pars.imageSizeOutput[0])) %
		           (long) pars.imageSizeOutput[0];
		ystart2 = ((long) pars.imageSizeOutput[0] + (y1 + ysplit % (long) pars.imageSizeOutput[0])) %
		           (long) pars.imageSizeOutput[0];

		cudaErrchk(cudaMemcpy2DAsync(permutedScompact_ds,
		                             pars.imageSizeReduce[1] * pars.numberBeams * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT),
		                             &permutedScompact_ph[ystart1 * pars.numberBeams * pars.Scompact.get_dimi() +
		                                                   xstart1 * pars.numberBeams],
		                             pars.Scompact.get_dimi() * pars.numberBeams * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT), // corresponds to stride between permuted Scompact elements in k-direction
		                             xsplit * pars.numberBeams * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT),
		                             ysplit,
		                             cudaMemcpyHostToDevice,
		                             stream));
		if (nx2 > 0 ) {
			cudaErrchk(cudaMemcpy2DAsync(&permutedScompact_ds[xsplit * pars.numberBeams],
			                             pars.imageSizeReduce[1] * pars.numberBeams * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT),
			                             &permutedScompact_ph[ystart1 * pars.numberBeams * pars.Scompact.get_dimi() +
			                                                   xstart2 * pars.numberBeams],
			                             pars.Scompact.get_dimi() * pars.numberBeams * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT), // corresponds to stride between permuted Scompact elements in k-direction
			                             nx2 * pars.numberBeams * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT),
			                             ysplit,
			                             cudaMemcpyHostToDevice,
			                             stream));
		}
		if (ny2 > 0 ) {
			cudaErrchk(cudaMemcpy2DAsync(&permutedScompact_ds[ysplit * pars.imageSizeReduce[1] * pars.numberBeams],
			                             pars.imageSizeReduce[1] * pars.numberBeams * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT),
			                             &permutedScompact_ph[ystart2 * pars.numberBeams * pars.Scompact.get_dimi() +
			                                                   xstart1 * pars.numberBeams],
			                             pars.Scompact.get_dimi() * pars.numberBeams * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT), // corresponds to stride between permuted Scompact elements in k-direction
			                             xsplit * pars.numberBeams * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT),
			                             ny2,
			                             cudaMemcpyHostToDevice,
			                             stream));
		}
		if (ny2 > 0 & nx2 > 0) {
			cudaErrchk(cudaMemcpy2DAsync(&permutedScompact_ds[ysplit * pars.imageSizeReduce[1] * pars.numberBeams +
			                                                   xsplit * pars.numberBeams],
			                             pars.imageSizeReduce[1] * pars.numberBeams * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT),
			                             &permutedScompact_ph[ystart2 * pars.numberBeams * pars.Scompact.get_dimi() +
			                                                   xstart2 * pars.numberBeams],
			                             pars.Scompact.get_dimi() * pars.numberBeams * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT), // corresponds to stride between permuted Scompact elements in k-direction
			                             nx2 * pars.numberBeams * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT),
			                             ny2,
			                             cudaMemcpyHostToDevice,
			                             stream));
		}

		// The data is now copied and we can proceed with the actual calculation

		// re-center the indices
		resetIndices <<<(pars.imageSizeReduce[0] - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream>>> (
				y_ds, pars.imageSizeReduce[0]);

		resetIndices <<<(pars.imageSizeReduce[1] - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream>>> (
				x_ds, pars.imageSizeReduce[1]);

		// Choose a good launch configuration
		// Heuristically use 2^p / 2 as the block size where p is the first power of 2 greater than the number of elements to work on.
		// This balances having enough work per thread and enough blocks without having so many blocks that the shared memory doesn't last long

		size_t p = getNextPower2(pars.numberBeams);
		const size_t BlockSize_numBeams = min((size_t)pars.deviceProperties.maxThreadsPerBlock,(size_t)std::max(1.0, pow(2,p) / 2));

		// Determine maximum threads per streaming multiprocessor based on the compute capability of the device
		size_t max_threads_per_sm;
		if (pars.deviceProperties.major > 3){
			max_threads_per_sm = 2048;
		} else if (pars.deviceProperties.major > 2) {
			max_threads_per_sm = 1536;
		} else {
			max_threads_per_sm = pars.deviceProperties.minor == 0 ? 768 : 1024;
		}

		// Estimate max number of simultaneous blocks per streaming multiprocessor
		const size_t max_blocks_per_sm = std::min((size_t)32, max_threads_per_sm / BlockSize_numBeams);

		// We find providing around 3 times as many blocks as the estimated maximum provides good performance
		const size_t target_blocks_per_sm = max_blocks_per_sm * 3;
		const size_t total_blocks         = target_blocks_per_sm * pars.deviceProperties.multiProcessorCount;

		// Determine amount of shared memory needed
		const unsigned long smem = pars.numberBeams * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT);

		if (BlockSize_numBeams >= 64) {

			const PRISMATIC_FLOAT_PRECISION aspect_ratio = (PRISMATIC_FLOAT_PRECISION)pars.imageSizeReduce[1] / (PRISMATIC_FLOAT_PRECISION)pars.imageSizeReduce[0];
			const size_t GridSizeZ = std::floor(sqrt(total_blocks / aspect_ratio));
			const size_t GridSizeY = aspect_ratio * GridSizeZ;
			dim3 grid(1, GridSizeY, GridSizeZ);
			dim3 block(BlockSize_numBeams, 1, 1);
//		// Launch kernel. Block size must be visible at compile time so we use a switch statement
			switch (BlockSize_numBeams) {
				case 1024 :
					scaleReduceS<1024> << < grid, block, smem, stream >> > (
							permutedScompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
							pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
				case 512 :
					scaleReduceS<512> << < grid, block, smem, stream >> > (
							permutedScompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
							pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
				case 256 :
					scaleReduceS<256> << < grid, block, smem, stream >> > (
							permutedScompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
							pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
				case 128 :
					scaleReduceS<128> << < grid, block, smem, stream >> > (
							permutedScompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
							pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
				case 64 :
					scaleReduceS<64> << < grid, block, smem, stream >> > (
							permutedScompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
							pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
				case 32 :
					scaleReduceS<32> << < grid, block, smem, stream >> > (
							permutedScompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
							pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
				case 16 :
					scaleReduceS<16> << < grid, block, smem, stream >> > (
							permutedScompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
							pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
				case 8 :
					scaleReduceS<8> << < grid, block, smem, stream >> > (
							permutedScompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
							pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
				case 4 :
					scaleReduceS<4> << < grid, block, smem, stream >> > (
							permutedScompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
							pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
				default :
					scaleReduceS<2> << < grid, block, smem, stream >> > (
							permutedScompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
							pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
			}
		} else {
			const size_t BlockSize_alongArray = 512 / BlockSize_numBeams;
			const size_t GridSize_alongArray = total_blocks;
			dim3 grid(1, GridSize_alongArray, 1);
			dim3 block(BlockSize_numBeams, BlockSize_alongArray, 1);
			// Determine amount of shared memory needed
			const unsigned long smem = pars.numberBeams * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT);

			// Launch kernel. Block size must be visible at compile time so we use a switch statement
			switch (BlockSize_numBeams) {
				case 1024 :
					scaleReduceS<1024, 1> << < grid, block, smem, stream >> > (
                            permutedScompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
                            pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
				case 512 :
					scaleReduceS<512, 1> << < grid, block, smem, stream >> > (
                            permutedScompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
                            pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
				case 256 :
					scaleReduceS<256, 2> << < grid, block, smem, stream >> > (
                            permutedScompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
                            pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
				case 128 :
					scaleReduceS<128, 4> << < grid, block, smem, stream >> > (
                            permutedScompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
                            pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
				case 64 :
					scaleReduceS<64, 8> << < grid, block, smem, stream >> > (
                            permutedScompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
                            pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
				case 32 :
					scaleReduceS<32, 16> << < grid, block, smem, stream >> > (
                            permutedScompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
                            pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
				case 16 :
					scaleReduceS<16, 32> << < grid, block, smem, stream >> > (
                            permutedScompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
                            pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
				case 8 :
					scaleReduceS<8, 64> << < grid, block, smem, stream >> > (
                            permutedScompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
                            pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
				case 4 :
					scaleReduceS<4, 128> << < grid, block, smem, stream >> > (
                            permutedScompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
                            pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
				default :
					scaleReduceS<1, 512> << < grid, block, smem, stream >> > (
                            permutedScompact_ds, phaseCoeffs_ds, psi_ds, y_ds, x_ds, pars.numberBeams, pars.imageSizeReduce[0],
                            pars.imageSizeReduce[1], pars.imageSizeReduce[0], pars.imageSizeReduce[1]);break;
			}
		}
		// final fft
		cufftErrchk(PRISMATIC_CUFFT_EXECUTE(cufft_plan, &psi_ds[0], &psi_ds[0], CUFFT_FORWARD));

		// convert to squared intensity
		abs_squared <<< (psi_size - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >>> (psiIntensity_ds, psi_ds, psi_size);

		// output calculation result
		size_t write_ay = (pars.meta.arbitraryProbes) ? 0 : ay;
		if(pars.meta.saveComplexOutputWave)
		{
			formatOutput_GPU_c_integrate(pars, psi_ds, psiIntensity_ds, alphaInd_d, output_ph,
				integratedOutput_ds, qyaReduce_d, qxaReduce_d, 0, write_ay, ax, pars.imageSizeReduce[0],
				pars.imageSizeReduce[1], stream, pars.scale);	
		}
		else
		{
			formatOutput_GPU_integrate(pars, psiIntensity_ds, alphaInd_d, output_ph,
				integratedOutput_ds, qyaReduce_d, qxaReduce_d, 0, write_ay, ax, pars.imageSizeReduce[0],
				pars.imageSizeReduce[1], stream, pars.scale);
		}
	}

}