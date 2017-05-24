// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

// For variable naming, the suffixes are "_d" for "device" (1 copy per GPU), "_ds" for "device stream (1 copy per stream), "_ph" for "pinned host"

#include "PRISM02.cuh"
#include "PRISM02.h"
#include <thread>
#include "WorkDispatcher.h"
#include "cufft.h"
#include "utility.cuh"
#include "params.cuh"


namespace PRISM {
	using namespace std;
	inline void createStreamsAndPlans2(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                                  CudaParameters<PRISM_FLOAT_PRECISION> &cuda_pars){
		// create CUDA streams
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;
		cuda_pars.streams 		    = new cudaStream_t[total_num_streams];
		cuda_pars.cufft_plans		= new cufftHandle[total_num_streams];
		cuda_pars.cufft_plans_small = new cufftHandle[total_num_streams];

		// batch parameters for cuFFT
		const int rank      = 2;
		int n[]             = {(int)pars.imageSize[0], (int)pars.imageSize[1]};
		const int howmany   = pars.meta.batch_size_GPU;
		int idist           = n[0]*n[1];
		int odist           = n[0]*n[1];
		int istride         = 1;
		int ostride         = 1;
		int *inembed        = n;
		int *onembed        = n;

		int n_small[]       = {(int)pars.qyInd.size(), (int)pars.qxInd.size()};
		int idist_small     = n_small[0]*n_small[1];
		int odist_small     = n_small[0]*n_small[1];
		int *inembed_small  = n_small;
		int *onembed_small  = n_small;

		// create cuFFT plans and CUDA streams
		for (auto j = 0; j < total_num_streams; ++j) {
			cudaSetDevice(j % pars.meta.NUM_GPUS);
			cudaErrchk(cudaStreamCreate(&cuda_pars.streams[j]));
			cufftErrchk(cufftPlanMany(&cuda_pars.cufft_plans[j], rank, n, inembed, istride, idist, onembed, ostride, odist, PRISM_CUFFT_PLAN_TYPE, howmany));
			cufftErrchk(cufftPlanMany(&cuda_pars.cufft_plans_small[j], rank, n_small, inembed_small, istride, idist_small, onembed_small, ostride, odist_small, PRISM_CUFFT_PLAN_TYPE, howmany));
			cufftErrchk(cufftSetStream(cuda_pars.cufft_plans[j], cuda_pars.streams[j]));
			cufftErrchk(cufftSetStream(cuda_pars.cufft_plans_small[j], cuda_pars.streams[j]));
		}
	}

	inline void allocatePinnedHostMemory_singlexfer2(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                                          CudaParameters<PRISM_FLOAT_PRECISION> &cuda_pars){
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;

		// allocate pinned memory
		cuda_pars.Scompact_slice_ph = new std::complex<PRISM_FLOAT_PRECISION>*[total_num_streams];
		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaMallocHost((void **) &cuda_pars.Scompact_slice_ph[s],
			                          pars.Scompact.get_dimj() * pars.Scompact.get_dimi() *
			                          sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
		}
		cudaErrchk(cudaMallocHost((void **) &cuda_pars.trans_ph,      pars.transmission.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
		cudaErrchk(cudaMallocHost((void **) &cuda_pars.prop_ph,       pars.prop.size()         * sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
		cudaErrchk(cudaMallocHost((void **) &cuda_pars.qxInd_ph,      pars.qxInd.size()        * sizeof(size_t)));
		cudaErrchk(cudaMallocHost((void **) &cuda_pars.qyInd_ph,      pars.qyInd.size()        * sizeof(size_t)));
		cudaErrchk(cudaMallocHost((void **) &cuda_pars.beamsIndex_ph, pars.beamsIndex.size()   * sizeof(size_t)));

	}

	inline void allocatePinnedHostMemory_streaming2(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                                         CudaParameters<PRISM_FLOAT_PRECISION> &cuda_pars){
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;

		// allocate pinned memory
		cuda_pars.Scompact_slice_ph = new std::complex<PRISM_FLOAT_PRECISION>*[total_num_streams];
		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaMallocHost((void **) &cuda_pars.Scompact_slice_ph[s],
			                          pars.Scompact.get_dimj() * pars.Scompact.get_dimi() *
			                          sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
		}
		cudaErrchk(cudaMallocHost((void **) &cuda_pars.trans_ph,      pars.transmission.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
		cudaErrchk(cudaMallocHost((void **) &cuda_pars.prop_ph,       pars.prop.size()         * sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
		cudaErrchk(cudaMallocHost((void **) &cuda_pars.qxInd_ph,      pars.qxInd.size()        * sizeof(size_t)));
		cudaErrchk(cudaMallocHost((void **) &cuda_pars.qyInd_ph,      pars.qyInd.size()        * sizeof(size_t)));
		cudaErrchk(cudaMallocHost((void **) &cuda_pars.beamsIndex_ph, pars.beamsIndex.size()   * sizeof(size_t)));
	}

	inline void copyToPinnedMemory_singlexfer2(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                                          CudaParameters<PRISM_FLOAT_PRECISION> &cuda_pars){
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;

		// copy host memory to pinned
		for (auto s = 0; s < total_num_streams; ++s) {
			memset(cuda_pars.Scompact_slice_ph[s], 0, pars.Scompact.get_dimj() * pars.Scompact.get_dimi() *
			                                          sizeof(std::complex<PRISM_FLOAT_PRECISION>));
		}
		memcpy(cuda_pars.trans_ph,      &pars.transmission[0], pars.transmission.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>));
		memcpy(cuda_pars.prop_ph,       &pars.prop[0],         pars.prop.size()         * sizeof(std::complex<PRISM_FLOAT_PRECISION>));
		memcpy(cuda_pars.qxInd_ph,      &pars.qxInd[0],        pars.qxInd.size()        * sizeof(size_t));
		memcpy(cuda_pars.qyInd_ph,      &pars.qyInd[0],        pars.qyInd.size()        * sizeof(size_t));
		memcpy(cuda_pars.beamsIndex_ph, &pars.beamsIndex[0],   pars.beamsIndex.size()   * sizeof(size_t));
	}

	inline void copyToPinnedMemory_streaming2(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                                         CudaParameters<PRISM_FLOAT_PRECISION> &cuda_pars){
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;

		// copy host memory to pinned
		for (auto s = 0; s < total_num_streams; ++s) {
			memset(cuda_pars.Scompact_slice_ph[s], 0, pars.Scompact.get_dimj() * pars.Scompact.get_dimi() *
			                                          sizeof(std::complex<PRISM_FLOAT_PRECISION>));
		}
		memcpy(cuda_pars.trans_ph,      &pars.transmission[0], pars.transmission.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>));
		memcpy(cuda_pars.prop_ph,       &pars.prop[0],         pars.prop.size()         * sizeof(std::complex<PRISM_FLOAT_PRECISION>));
		memcpy(cuda_pars.qxInd_ph,      &pars.qxInd[0],        pars.qxInd.size()        * sizeof(size_t));
		memcpy(cuda_pars.qyInd_ph,      &pars.qyInd[0],        pars.qyInd.size()        * sizeof(size_t));
		memcpy(cuda_pars.beamsIndex_ph, &pars.beamsIndex[0],   pars.beamsIndex.size()   * sizeof(size_t));
	}

	inline void allocateDeviceMemory_singlexfer2(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                                            CudaParameters<PRISM_FLOAT_PRECISION> &cuda_pars){
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;

		// allocate memory on the device

		// pointers to read-only GPU memory (one copy per GPU)
		cuda_pars.trans_d       = new PRISM_CUDA_COMPLEX_FLOAT*[pars.meta.NUM_GPUS];
		cuda_pars.prop_d        = new PRISM_CUDA_COMPLEX_FLOAT*[pars.meta.NUM_GPUS];
		cuda_pars.qxInd_d       = new size_t*[pars.meta.NUM_GPUS];
		cuda_pars.qyInd_d       = new size_t*[pars.meta.NUM_GPUS];
		cuda_pars.beamsIndex_d  = new size_t*[pars.meta.NUM_GPUS];

//		// pointers to read/write GPU memory (one per stream)
		cuda_pars.psi_ds       = new PRISM_CUDA_COMPLEX_FLOAT*[total_num_streams];
		cuda_pars.psi_small_ds = new PRISM_CUDA_COMPLEX_FLOAT*[total_num_streams];

		// allocate memory on each GPU
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.trans_d[g],      pars.transmission.size() * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.prop_d[g],       pars.prop.size()         * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.qxInd_d[g],      pars.qxInd.size()        * sizeof(size_t)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.qyInd_d[g],      pars.qyInd.size()        * sizeof(size_t)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.beamsIndex_d[g], pars.beamsIndex.size()   * sizeof(size_t)));
		}

		// allocate memory per stream and 0 it
		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaSetDevice(s % pars.meta.NUM_GPUS));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.psi_ds[s],
			                      pars.meta.batch_size_GPU*pars.imageSize[0] * pars.imageSize[1] * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.psi_small_ds[s],
			                      pars.meta.batch_size_GPU*pars.qxInd.size() * pars.qyInd.size() * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMemset(cuda_pars.psi_ds[s], 0,
			                      pars.meta.batch_size_GPU*pars.imageSize[0] * pars.imageSize[1] * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMemset(cuda_pars.psi_small_ds[s], 0,
			                      pars.meta.batch_size_GPU*pars.qxInd.size() * pars.qyInd.size() * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
		}
	}

	inline void allocateDeviceMemory_streaming2(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                                           CudaParameters<PRISM_FLOAT_PRECISION> &cuda_pars){
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;
		// pointers to read-only GPU memory (one copy per GPU)
		cuda_pars.prop_d  = new PRISM_CUDA_COMPLEX_FLOAT*[pars.meta.NUM_GPUS];
		cuda_pars.qxInd_d                   = new size_t*[pars.meta.NUM_GPUS];
		cuda_pars.qyInd_d                   = new size_t*[pars.meta.NUM_GPUS];
		cuda_pars.beamsIndex_d              = new size_t*[pars.meta.NUM_GPUS];

		// pointers to read/write GPU memory (one per stream)
		cuda_pars.psi_ds       = new PRISM_CUDA_COMPLEX_FLOAT*[total_num_streams];
		cuda_pars.psi_small_ds = new PRISM_CUDA_COMPLEX_FLOAT*[total_num_streams];
		cuda_pars.trans_d    = new PRISM_CUDA_COMPLEX_FLOAT*[total_num_streams];

		// allocate memory on each GPU
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			cudaErrchk(cudaSetDevice(g));

			cudaErrchk(cudaMalloc((void **) &cuda_pars.prop_d[g], pars.prop.size() * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.qxInd_d[g], pars.qxInd.size() * sizeof(size_t)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.qyInd_d[g], pars.qyInd.size() * sizeof(size_t)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.beamsIndex_d[g], pars.beamsIndex.size() * sizeof(size_t)));
		}

		// allocate memory per stream and 0 it
		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaSetDevice(s % pars.meta.NUM_GPUS));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.trans_d[s],
			                      pars.imageSize[0] * pars.imageSize[1] * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.psi_ds[s],
			                      pars.meta.batch_size_GPU*pars.imageSize[0] * pars.imageSize[1] * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.psi_small_ds[s],
			                      pars.meta.batch_size_GPU*pars.qxInd.size() * pars.qyInd.size() * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMemset(cuda_pars.psi_ds[s], 0,
			                      pars.meta.batch_size_GPU*pars.imageSize[0] * pars.imageSize[1] * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMemset(cuda_pars.psi_small_ds[s], 0,
			                      pars.meta.batch_size_GPU*pars.qxInd.size() * pars.qyInd.size() * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
		}
	}

	inline void copyToDeviceMemory_singlexfer2(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                                          CudaParameters<PRISM_FLOAT_PRECISION> &cuda_pars){
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;
		// Copy memory to each GPU asynchronously from the pinned host memory spaces.
		// The streams are laid out so that consecutive streams represent different GPUs. If we
		// have more than one stream per GPU, then we want to interleave as much as possible
		int stream_id = 0;
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			stream_id = g;
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaMemcpyAsync(cuda_pars.trans_d[g], &cuda_pars.trans_ph[0],
			                           pars.transmission.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>),
			                           cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));

			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(cuda_pars.prop_d[g], &cuda_pars.prop_ph[0],
			                           pars.prop.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>),
			                           cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(cuda_pars.qxInd_d[g], &cuda_pars.qxInd_ph[0],
			                           pars.qxInd.size() * sizeof(size_t), cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(cuda_pars.qyInd_d[g], &cuda_pars.qyInd_ph[0],
			                           pars.qyInd.size() * sizeof(size_t), cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(cuda_pars.beamsIndex_d[g], &cuda_pars.beamsIndex_ph[0],
			                           pars.beamsIndex.size() * sizeof(size_t), cudaMemcpyHostToDevice,
			                           cuda_pars.streams[stream_id]));
		}

		// make sure transfers are complete
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			cudaSetDevice(g);
			cudaDeviceSynchronize();
		}
	}

	inline void copyToDeviceMemory_streaming2(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                                         CudaParameters<PRISM_FLOAT_PRECISION> &cuda_pars){
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;
		// Copy memory to each GPU asynchronously from the pinned host memory spaces.
		// The streams are laid out so that consecutive streams represent different GPUs. If we
		// have more than one stream per GPU, then we want to interleave as much as possible
		int stream_id = 0;
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			stream_id = g;
			cudaErrchk(cudaSetDevice(g));
//			cudaErrchk(cudaMemcpyAsync(trans_ds[g], &trans_ph[0],
//			                           pars.imageSize[0] * pars.imageSize[1] * sizeof(std::complex<PRISM_FLOAT_PRECISION>),
//			                           cudaMemcpyHostToDevice, streams[stream_id]));

			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(cuda_pars.prop_d[g], &cuda_pars.prop_ph[0],
			                           pars.prop.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>),
			                           cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(cuda_pars.qxInd_d[g], &cuda_pars.qxInd_ph[0],
			                           pars.qxInd.size() * sizeof(size_t), cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(cuda_pars.qyInd_d[g], &cuda_pars.qyInd_ph[0],
			                           pars.qyInd.size() * sizeof(size_t), cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(cuda_pars.beamsIndex_d[g], &cuda_pars.beamsIndex_ph[0],
			                           pars.beamsIndex.size() * sizeof(size_t), cudaMemcpyHostToDevice,
			                           cuda_pars.streams[stream_id]));
		}

		// make sure transfers are complete
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			cudaSetDevice(g);
			cudaDeviceSynchronize();
		}
	}
	inline void launchWorkers_singlexfer2(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                                              CudaParameters<PRISM_FLOAT_PRECISION> &cuda_pars){
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;
		// launch GPU work
		vector<thread> workers_GPU;
		workers_GPU.reserve(total_num_streams); // prevents multiple reallocations
		int stream_count = 0;
		const size_t PRISM_PRINT_FREQUENCY_BEAMS = max((size_t)1,pars.numberBeams / 10); // for printing status
		WorkDispatcher dispatcher(0, pars.numberBeams);
		for (auto t = 0; t < total_num_streams; ++t) {

			int GPU_num = stream_count % pars.meta.NUM_GPUS; // determine which GPU handles this job
			cudaSetDevice(GPU_num);
			cudaStream_t &current_stream = cuda_pars.streams[stream_count];
			cout << "Launching GPU worker on stream #" << stream_count << " of GPU #" << GPU_num << endl;

			// get pointers to the pre-copied arrays, making sure to get those on the current GPU
			PRISM_CUDA_COMPLEX_FLOAT *current_trans_d = cuda_pars.trans_d[GPU_num];
			PRISM_CUDA_COMPLEX_FLOAT *current_prop_d = cuda_pars.prop_d[GPU_num];
			size_t *current_qxInd_d = cuda_pars.qxInd_d[GPU_num];
			size_t *current_qyInd_d = cuda_pars.qyInd_d[GPU_num];
			size_t *current_beamsIndex = cuda_pars.beamsIndex_d[GPU_num];
			// get pointers to per-stream arrays
			PRISM_CUDA_COMPLEX_FLOAT *current_psi_ds = cuda_pars.psi_ds[stream_count];
			PRISM_CUDA_COMPLEX_FLOAT *current_psi_small_ds = cuda_pars.psi_small_ds[stream_count];
			cufftHandle &current_cufft_plan = cuda_pars.cufft_plans[stream_count];
			cufftHandle &current_cufft_plan_small = cuda_pars.cufft_plans_small[stream_count];
			complex<PRISM_FLOAT_PRECISION> *current_S_slice_ph = cuda_pars.Scompact_slice_ph[stream_count];

			workers_GPU.push_back(thread([&pars, current_trans_d, current_prop_d, current_qxInd_d, current_qyInd_d, &dispatcher,
					                             current_psi_ds, current_psi_small_ds, &current_cufft_plan, &current_cufft_plan_small,
					                             current_S_slice_ph, current_beamsIndex, GPU_num, stream_count, &current_stream, &PRISM_PRINT_FREQUENCY_BEAMS]() {
				cudaErrchk(cudaSetDevice(GPU_num));

#ifndef NDEBUG
				{
//					 check memory usage on the GPU
					std::lock_guard<mutex> lock(PRISM::mem_lock);
					size_t free_mem, total_mem;
					free_mem=total_mem=0;
					cudaErrchk(cudaMemGetInfo(&free_mem, &total_mem));
					pars.max_mem = std::max(total_mem - free_mem, pars.max_mem);
//					cout << "max_mem = " << pars.max_mem << endl;
				}
#endif // NDEBUG

				size_t currentBeam, stopBeam;
				currentBeam=stopBeam=0;
//				while (getWorkID(pars, currentBeam, stopBeam)) {
				while (dispatcher.getWork(currentBeam, stopBeam, pars.meta.batch_size_GPU)) {
					while (currentBeam < stopBeam) {
						if (currentBeam % PRISM_PRINT_FREQUENCY_BEAMS < pars.meta.batch_size_GPU | currentBeam == 100){
							cout << "Computing Plane Wave #" << currentBeam << "/" << pars.numberBeams << endl;
						}
//						propagatePlaneWave_GPU_singlexfer(pars,
//						                                  current_trans_d,
//						                                  current_psi_ds,
//						                                  current_psi_small_ds,
//						                                  current_S_slice_ph,
//						                                  current_qyInd_d,
//						                                  current_qxInd_d,
//						                                  current_prop_d,
//						                                  current_beamsIndex,
//						                                  currentBeam,
//						                                  current_cufft_plan,
//						                                  current_cufft_plan_small,
//						                                  current_stream);
						propagatePlaneWave_GPU_singlexfer_batch(pars,
						                                        current_trans_d,
						                                        current_psi_ds,
						                                        current_psi_small_ds,
						                                        current_S_slice_ph,
						                                        current_qyInd_d,
						                                        current_qxInd_d,
						                                        current_prop_d,
						                                        current_beamsIndex,
						                                        currentBeam,
						                                        stopBeam,
						                                        current_cufft_plan,
						                                        current_cufft_plan_small,
						                                        current_stream);
#ifdef PRISM_BUILDING_GUI
						pars.progressbar->signalScompactUpdate(currentBeam, pars.numberBeams);
#endif
//						++currentBeam;
						currentBeam=stopBeam;
					}
				}
				cout << "GPU worker on stream #" << stream_count << " of GPU #" << GPU_num << "finished\n";
			}));
			++stream_count;
		}

		if (pars.meta.also_do_CPU_work){

			// launch CPU work
			vector<thread> workers_CPU;
			workers_CPU.reserve(pars.meta.NUM_THREADS); // prevents multiple reallocations
			mutex fftw_plan_lock;
			pars.meta.batch_size_CPU = min(pars.meta.batch_size_target_CPU, max((size_t)1, pars.numberBeams / pars.meta.NUM_THREADS));
			cout << "PRISM02 pars.meta.batch_size_CPU = " << pars.meta.batch_size_CPU << endl;
			PRISM_FFTW_INIT_THREADS();
			PRISM_FFTW_PLAN_WITH_NTHREADS(pars.meta.NUM_THREADS);
			for (auto t = 0; t < pars.meta.NUM_THREADS; ++t) {
				cout << "Launching thread #" << t << " to compute beams\n";
				workers_CPU.push_back(thread([&pars, &fftw_plan_lock, &dispatcher, &PRISM_PRINT_FREQUENCY_BEAMS]() {

					size_t currentBeam, stopBeam, early_CPU_stop;
					currentBeam=stopBeam=0;
					if (pars.meta.NUM_GPUS > 0){
						// if there are no GPUs, make sure to do all work on CPU
//						early_CPU_stop = (size_t)std::max((PRISM_FLOAT_PRECISION)0.0,pars.numberBeams - pars.meta.gpu_cpu_ratio);
						early_CPU_stop = (size_t)std::max((PRISM_FLOAT_PRECISION)0.0,pars.numberBeams - pars.meta.gpu_cpu_ratio * pars.meta.batch_size_CPU);
					} else {
						early_CPU_stop = pars.numberBeams;
					}
					if (dispatcher.getWork(currentBeam, stopBeam, pars.meta.batch_size_CPU, early_CPU_stop)) {
						// allocate array for psi just once per thread
						Array1D<complex<PRISM_FLOAT_PRECISION> > psi_stack = zeros_ND<1, complex<PRISM_FLOAT_PRECISION> >(
								{{pars.imageSize[0]*pars.imageSize[1]*pars.meta.batch_size_CPU}});

//						 setup batch FFTW parameters
						const int rank = 2;
						int n[] = {(int)pars.imageSize[0], (int)pars.imageSize[1]};
						const int howmany = pars.meta.batch_size_CPU;
						int idist = n[0]*n[1];
						int odist = n[0]*n[1];
						int istride = 1;
						int ostride = 1;
						int *inembed = n;
						int *onembed = n;



						unique_lock<mutex> gatekeeper(fftw_plan_lock);
						PRISM_FFTW_PLAN plan_forward = PRISM_FFTW_PLAN_DFT_BATCH(rank, n, howmany,
						                                                         reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi_stack[0]), inembed,
						                                                         istride, idist,
						                                                         reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi_stack[0]), onembed,
						                                                         ostride, odist,
						                                                         FFTW_FORWARD, FFTW_MEASURE);
						PRISM_FFTW_PLAN plan_inverse = PRISM_FFTW_PLAN_DFT_BATCH(rank, n, howmany,
						                                                         reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi_stack[0]), inembed,
						                                                         istride, idist,
						                                                         reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi_stack[0]), onembed,
						                                                         ostride, odist,
						                                                         FFTW_BACKWARD, FFTW_MEASURE);

						gatekeeper.unlock(); // unlock it so we only block as long as necessary to deal with plans

						// main work loop
						do { // synchronously get work assignment
							while (currentBeam < stopBeam) {
								if (currentBeam % PRISM_PRINT_FREQUENCY_BEAMS < pars.meta.batch_size_CPU | currentBeam == 100){
									cout << "Computing Plane Wave #" << currentBeam << "/" << pars.numberBeams << endl;
								}
								// re-zero psi each iteration
								memset((void *) &psi_stack[0], 0, psi_stack.size() * sizeof(complex<PRISM_FLOAT_PRECISION>));
//								propagatePlaneWave_CPU(pars, currentBeam, psi, plan_forward, plan_inverse, fftw_plan_lock);
								propagatePlaneWave_CPU_batch(pars, currentBeam, stopBeam, psi_stack, plan_forward, plan_inverse, fftw_plan_lock);
#ifdef PRISM_BUILDING_GUI
								pars.progressbar->signalScompactUpdate(currentBeam, pars.numberBeams);
#endif
//                                currentBeam = stopBeam;
								++currentBeam;
							}
							if (currentBeam >= early_CPU_stop) break;
						} while (dispatcher.getWork(currentBeam, stopBeam, pars.meta.batch_size_CPU, early_CPU_stop));
						// clean up
						gatekeeper.lock();
						PRISM_FFTW_DESTROY_PLAN(plan_forward);
						PRISM_FFTW_DESTROY_PLAN(plan_inverse);
						gatekeeper.unlock();
					}
				}));
			}
			for (auto &t:workers_CPU)t.join();
			PRISM_FFTW_CLEANUP_THREADS();
		}

		for (auto &t:workers_GPU)t.join();
	}
	 void launchWorkers_streaming(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                                    CudaParameters<PRISM_FLOAT_PRECISION> &cuda_pars){
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;
		// launch GPU work
		vector<thread> workers_GPU;
		workers_GPU.reserve(total_num_streams); // prevents multiple reallocations
		int stream_count = 0;
		const size_t PRISM_PRINT_FREQUENCY_BEAMS = max((size_t)1,pars.numberBeams / 10); // for printing status
		WorkDispatcher dispatcher(0, pars.numberBeams);
		for (auto t = 0; t < total_num_streams; ++t) {
			int GPU_num = stream_count % pars.meta.NUM_GPUS; // determine which GPU handles this job
			cudaSetDevice(GPU_num);
			cudaStream_t &current_stream = cuda_pars.streams[stream_count];
			// get pointers to the pre-copied arrays, making sure to get those on the current GPU
			PRISM_CUDA_COMPLEX_FLOAT *current_prop_d = cuda_pars.prop_d[GPU_num];
			size_t *current_qxInd_d = cuda_pars.qxInd_d[GPU_num];
			size_t *current_qyInd_d = cuda_pars.qyInd_d[GPU_num];
			size_t *current_beamsIndex = cuda_pars.beamsIndex_d[GPU_num];
			// get pointers to per-stream arrays
			PRISM_CUDA_COMPLEX_FLOAT *current_trans_ds = cuda_pars.trans_d[stream_count];
			PRISM_CUDA_COMPLEX_FLOAT *current_psi_ds = cuda_pars.psi_ds[stream_count];
			PRISM_CUDA_COMPLEX_FLOAT *current_psi_small_ds = cuda_pars.psi_small_ds[stream_count];
			cufftHandle &current_cufft_plan = cuda_pars.cufft_plans[stream_count];
			cufftHandle &current_cufft_plan_small = cuda_pars.cufft_plans_small[stream_count];
			complex<PRISM_FLOAT_PRECISION> *current_S_slice_ph = cuda_pars.Scompact_slice_ph[stream_count];

			workers_GPU.push_back(thread([&pars, current_trans_ds, current_prop_d, current_qxInd_d, current_qyInd_d, &dispatcher,
					                             current_psi_ds, current_psi_small_ds, &current_cufft_plan, &current_cufft_plan_small,
					                             current_S_slice_ph, current_beamsIndex, GPU_num, stream_count, &current_stream, &PRISM_PRINT_FREQUENCY_BEAMS, &cuda_pars]() {
				cudaErrchk(cudaSetDevice(GPU_num));

#ifndef NDEBUG
				{
//					 check memory usage on the GPU
					std::lock_guard<mutex> lock(PRISM::mem_lock);
					size_t free_mem, total_mem;
					free_mem=total_mem=0;
					cudaErrchk(cudaMemGetInfo(&free_mem, &total_mem));
					pars.max_mem = std::max(total_mem - free_mem, pars.max_mem);
					cout << "max_mem = " << pars.max_mem << endl;
				}
#endif // NDEBUG

				size_t currentBeam, stopBeam;
				currentBeam=stopBeam=0;
//				while (getWorkID(pars, currentBeam, stopBeam)) {
				while (dispatcher.getWork(currentBeam, stopBeam, pars.meta.batch_size_GPU)) {
					while (currentBeam < stopBeam) {
						if (currentBeam % PRISM_PRINT_FREQUENCY_BEAMS < pars.meta.batch_size_GPU | currentBeam == 100){
							cout << "Computing Plane Wave #" << currentBeam << "/" << pars.numberBeams << endl;
						}
//						propagatePlaneWave_GPU_streaming(pars,
//						                                 current_trans_ds,
//						                                 trans_ph,
//						                                 current_psi_ds,
//						                                 current_psi_small_ds,
//						                                 current_S_slice_ph,
//						                                 current_qyInd_d,
//						                                 current_qxInd_d,
//						                                 current_prop_d,
//						                                 current_beamsIndex,
//						                                 currentBeam,
//						                                 current_cufft_plan,
//						                                 current_cufft_plan_small,
//						                                 current_stream);
						propagatePlaneWave_GPU_streaming_batch(pars,
						                                       current_trans_ds,
						                                       cuda_pars.trans_ph,
						                                       current_psi_ds,
						                                       current_psi_small_ds,
						                                       current_S_slice_ph,
						                                       current_qyInd_d,
						                                       current_qxInd_d,
						                                       current_prop_d,
						                                       current_beamsIndex,
						                                       currentBeam,
						                                       stopBeam,
						                                       current_cufft_plan,
						                                       current_cufft_plan_small,
						                                       current_stream);
//						++currentBeam;
						currentBeam=stopBeam;
#ifdef PRISM_BUILDING_GUI
						pars.progressbar->signalScompactUpdate(currentBeam, pars.numberBeams);
#endif
					}
				}
				cout << "GPU worker on stream #" << stream_count << " of GPU #" << GPU_num << " finished\n";
			}));
			++stream_count;
		}

		if (pars.meta.also_do_CPU_work){

			// launch CPU work
			vector<thread> workers_CPU;
			workers_CPU.reserve(pars.meta.NUM_THREADS); // prevents multiple reallocations
			mutex fftw_plan_lock;
			pars.meta.batch_size_CPU = min(pars.meta.batch_size_target_CPU, max((size_t)1, pars.numberBeams / pars.meta.NUM_THREADS));
			cout << "PRISM02 pars.meta.batch_size_CPU = " << pars.meta.batch_size_CPU << endl;
			PRISM_FFTW_INIT_THREADS();
			PRISM_FFTW_PLAN_WITH_NTHREADS(pars.meta.NUM_THREADS);
			for (auto t = 0; t < pars.meta.NUM_THREADS; ++t) {
				cout << "Launching thread #" << t << " to compute beams\n";
				workers_CPU.push_back(thread([&pars, &fftw_plan_lock, &dispatcher, &PRISM_PRINT_FREQUENCY_BEAMS]() {

					size_t currentBeam, stopBeam, early_CPU_stop;
					currentBeam=stopBeam=0;
					if (pars.meta.NUM_GPUS > 0){
						// if there are no GPUs, make sure to do all work on CPU
						//early_CPU_stop = (size_t)std::max((PRISM_FLOAT_PRECISION)0.0,pars.numberBeams - pars.meta.gpu_cpu_ratio);
						early_CPU_stop = (size_t)std::max((PRISM_FLOAT_PRECISION)0.0,pars.numberBeams - pars.meta.gpu_cpu_ratio*pars.meta.batch_size_CPU);
					} else {
						early_CPU_stop = pars.numberBeams;
					}
					if (dispatcher.getWork(currentBeam, stopBeam, pars.meta.batch_size_CPU, early_CPU_stop)) {
						// allocate array for psi just once per thread
//						Array2D<complex<PRISM_FLOAT_PRECISION> > psi = zeros_ND<2, complex<PRISM_FLOAT_PRECISION> >(
//								{{pars.imageSize[0], pars.imageSize[1]}});
						Array1D<complex<PRISM_FLOAT_PRECISION> > psi_stack = zeros_ND<1, complex<PRISM_FLOAT_PRECISION> >(
								{{pars.imageSize[0]*pars.imageSize[1]*pars.meta.batch_size_CPU}});

						// setup batch FFTW parameters
						const int rank = 2;
						int n[] = {(int)pars.imageSize[0], (int)pars.imageSize[1]};
						const int howmany = pars.meta.batch_size_CPU;
						int idist = n[0]*n[1];
						int odist = n[0]*n[1];
						int istride = 1;
						int ostride = 1;
						int *inembed = n;
						int *onembed = n;

						unique_lock<mutex> gatekeeper(fftw_plan_lock);
						PRISM_FFTW_PLAN plan_forward = PRISM_FFTW_PLAN_DFT_BATCH(rank, n, howmany,
						                                                         reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi_stack[0]), inembed,
						                                                         istride, idist,
						                                                         reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi_stack[0]), onembed,
						                                                         ostride, odist,
						                                                         FFTW_FORWARD, FFTW_MEASURE);
						PRISM_FFTW_PLAN plan_inverse = PRISM_FFTW_PLAN_DFT_BATCH(rank, n, howmany,
						                                                         reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi_stack[0]), inembed,
						                                                         istride, idist,
						                                                         reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi_stack[0]), onembed,
						                                                         ostride, odist,
						                                                         FFTW_BACKWARD, FFTW_MEASURE);

						gatekeeper.unlock(); // unlock it so we only block as long as necessary to deal with plans

						// main work loop
						do { // synchronously get work assignment
							while (currentBeam < stopBeam) {
								if (currentBeam % PRISM_PRINT_FREQUENCY_BEAMS < pars.meta.batch_size_CPU | currentBeam == 100){
									cout << "Computing Plane Wave #" << currentBeam << "/" << pars.numberBeams << endl;
								}
								// re-zero psi each iteration
								memset((void *) &psi_stack[0], 0, psi_stack.size() * sizeof(complex<PRISM_FLOAT_PRECISION>));
//								propagatePlaneWave_CPU(pars, currentBeam, psi, plan_forward, plan_inverse, fftw_plan_lock);
								propagatePlaneWave_CPU_batch(pars, currentBeam, stopBeam, psi_stack, plan_forward, plan_inverse, fftw_plan_lock);
#ifdef PRISM_BUILDING_GUI
								pars.progressbar->signalScompactUpdate(currentBeam, pars.numberBeams);
#endif
								currentBeam = stopBeam;
//								++currentBeam;
							}
							if (currentBeam >= early_CPU_stop) break;
						} while (dispatcher.getWork(currentBeam, stopBeam, pars.meta.batch_size_CPU, early_CPU_stop));
						// clean up
						gatekeeper.lock();
						PRISM_FFTW_DESTROY_PLAN(plan_forward);
						PRISM_FFTW_DESTROY_PLAN(plan_inverse);
						gatekeeper.unlock();
					}
				}));
			}
			for (auto &t:workers_CPU)t.join();
			PRISM_FFTW_CLEANUP_THREADS();
		}
		for (auto &t:workers_GPU)t.join();
	}

	inline void cleanupMemory2(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                          CudaParameters<PRISM_FLOAT_PRECISION> &cuda_pars){
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaFree(cuda_pars.trans_d[g]));
			cudaErrchk(cudaFree(cuda_pars.prop_d[g]));
			cudaErrchk(cudaFree(cuda_pars.qxInd_d[g]));
			cudaErrchk(cudaFree(cuda_pars.qyInd_d[g]));
			cudaErrchk(cudaFree(cuda_pars.beamsIndex_d[g]));
		}

		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaSetDevice(s % pars.meta.NUM_GPUS));
			cudaErrchk(cudaFree(cuda_pars.psi_ds[s]));
			cudaErrchk(cudaFree(cuda_pars.psi_small_ds[s]));
			cufftErrchk(cufftDestroy(cuda_pars.cufft_plans[s]));
			cufftErrchk(cufftDestroy(cuda_pars.cufft_plans_small[s]));
		}

		// free pinned memory
		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaFreeHost(cuda_pars.Scompact_slice_ph[s]));
		}
		cudaErrchk(cudaFreeHost(cuda_pars.trans_ph));
		cudaErrchk(cudaFreeHost(cuda_pars.prop_ph));
		cudaErrchk(cudaFreeHost(cuda_pars.qxInd_ph));
		cudaErrchk(cudaFreeHost(cuda_pars.qyInd_ph));
		cudaErrchk(cudaFreeHost(cuda_pars.beamsIndex_ph));


		// destroy CUDA streams
		for (auto j = 0; j < total_num_streams; ++j){
			cudaSetDevice(j % pars.meta.NUM_GPUS);
			cudaErrchk(cudaStreamDestroy(cuda_pars.streams[j]));
		}

		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g){
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaDeviceReset());
		}

		delete[] cuda_pars.streams;
		delete[] cuda_pars.cufft_plans;
		delete[] cuda_pars.cufft_plans_small;
		delete[] cuda_pars.trans_d;
		delete[] cuda_pars.prop_d;
		delete[] cuda_pars.qxInd_d;
		delete[] cuda_pars.qyInd_d;
		delete[] cuda_pars.beamsIndex_d;
		delete[] cuda_pars.psi_ds;
		delete[] cuda_pars.psi_small_ds;
		delete[] cuda_pars.Scompact_slice_ph;
	}

	void propagatePlaneWave_GPU_singlexfer(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                                       PRISM_CUDA_COMPLEX_FLOAT* trans_d,
		  	                               PRISM_CUDA_COMPLEX_FLOAT* psi_d,
			                               PRISM_CUDA_COMPLEX_FLOAT* psi_small_d,
			                               complex<PRISM_FLOAT_PRECISION>* Scompact_slice_ph,
			                               const size_t* qyInd_d,
			                               const size_t* qxInd_d,
			                               const PRISM_CUDA_COMPLEX_FLOAT* prop_d,
			                               const size_t* beamsIndex,
			                               const size_t beamNumber,
			                               const cufftHandle& plan,
			                               const cufftHandle& plan_small,
			                               cudaStream_t& stream){

		const size_t psi_size = pars.imageSize[0] * pars.imageSize[1];
		const size_t psi_small_size = pars.qxInd.size() * pars.qyInd.size();
		initializePsi_oneNonzero<<< (psi_size - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream>>>(psi_d, psi_size, pars.beamsIndex[beamNumber]);


		for (auto planeNum = 0; planeNum < pars.numPlanes; ++planeNum) {
			cufftErrchk(PRISM_CUFFT_EXECUTE(plan, &psi_d[0], &psi_d[0], CUFFT_INVERSE));
			multiply_cx<<<(psi_size-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_d, &trans_d[planeNum*psi_size], psi_size);
			divide_inplace<<<(psi_size-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_d, PRISM_MAKE_CU_COMPLEX(psi_size, 0), psi_size);
			cufftErrchk(PRISM_CUFFT_EXECUTE(plan, &psi_d[0], &psi_d[0], CUFFT_FORWARD));
			multiply_cx<<<(psi_size-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_d, prop_d, psi_size);
		}

		array_subset<<<(pars.qyInd.size()*pars.qxInd.size()-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>> (
		psi_d, psi_small_d, qyInd_d, qxInd_d, pars.imageSize[1], pars.qyInd.size(), pars.qxInd.size());

		PRISM_CUFFT_EXECUTE(plan_small,&psi_small_d[0], &psi_small_d[0], CUFFT_INVERSE);
        divide_inplace<<<(psi_small_size-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_small_d, PRISM_MAKE_CU_COMPLEX(psi_small_size, 0),psi_small_size);

		cudaErrchk(cudaMemcpyAsync(Scompact_slice_ph,&psi_small_d[0],psi_small_size * sizeof(PRISM_CUDA_COMPLEX_FLOAT),cudaMemcpyDeviceToHost,stream));
		cudaStreamSynchronize(stream);
		memcpy(&pars.Scompact[beamNumber * pars.Scompact.get_dimj() * pars.Scompact.get_dimi()], &Scompact_slice_ph[0], psi_small_size * sizeof(PRISM_CUDA_COMPLEX_FLOAT));

	}

	void propagatePlaneWave_GPU_singlexfer_batch(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                                             PRISM_CUDA_COMPLEX_FLOAT* trans_d,
	                                             PRISM_CUDA_COMPLEX_FLOAT* psi_d,
	                                             PRISM_CUDA_COMPLEX_FLOAT* psi_small_d,
	                                             complex<PRISM_FLOAT_PRECISION>* Scompact_slice_ph,
	                                             const size_t* qyInd_d,
	                                             const size_t* qxInd_d,
	                                             const PRISM_CUDA_COMPLEX_FLOAT* prop_d,
	                                             const size_t* beamsIndex,
	                                             const size_t beamNumber,
	                                             const size_t stopBeam,
	                                             const cufftHandle& plan,
	                                             const cufftHandle& plan_small,
	                                             cudaStream_t& stream){

		const size_t psi_size        = pars.imageSize[0] * pars.imageSize[1];
		const size_t psi_small_size = pars.qxInd.size() * pars.qyInd.size();
		for (auto batch_idx = 0; batch_idx < (stopBeam-beamNumber); ++batch_idx) {
			// initialize psi
			initializePsi_oneNonzero<<< (psi_size - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream>>>(psi_d + batch_idx*psi_size, psi_size, pars.beamsIndex[beamNumber + batch_idx]);
		}
		for (auto planeNum = 0; planeNum < pars.numPlanes; ++planeNum) {
			cufftErrchk(PRISM_CUFFT_EXECUTE(plan, &psi_d[0], &psi_d[0], CUFFT_INVERSE));
			for (auto batch_idx = 0; batch_idx < (stopBeam-beamNumber); ++batch_idx) {
				multiply_cx << < (psi_size - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >> >
						(psi_d + batch_idx*psi_size, &trans_d[planeNum * psi_size], psi_size);
				divide_inplace<<<(psi_size - 1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>
						(psi_d + batch_idx*psi_size, PRISM_MAKE_CU_COMPLEX(psi_size, 0), psi_size);
			}
			cufftErrchk(PRISM_CUFFT_EXECUTE(plan, &psi_d[0], &psi_d[0], CUFFT_FORWARD));
			for (auto batch_idx = 0; batch_idx < (stopBeam-beamNumber); ++batch_idx) {
				multiply_cx << < (psi_size - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >> > (psi_d + batch_idx*psi_size, prop_d, psi_size);
			}
		}

		for (auto batch_idx = 0; batch_idx < (stopBeam-beamNumber); ++batch_idx) {
			array_subset << < (pars.qyInd.size() * pars.qxInd.size() - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0,
					stream >> > (psi_d + batch_idx*psi_size, psi_small_d + batch_idx*psi_small_size, qyInd_d, qxInd_d, pars.imageSize[1], pars.qyInd.size(), pars.qxInd.size());
		}

		PRISM_CUFFT_EXECUTE(plan_small,&psi_small_d[0], &psi_small_d[0], CUFFT_INVERSE);
		for (auto batch_idx = 0; batch_idx < (stopBeam-beamNumber); ++batch_idx) {
			divide_inplace<<<(psi_small_size-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_small_d + batch_idx*psi_small_size, PRISM_MAKE_CU_COMPLEX(psi_small_size, 0), psi_small_size);
		}
		for (auto batch_idx = 0; batch_idx < (stopBeam-beamNumber); ++batch_idx) {
			cudaErrchk(cudaMemcpyAsync(Scompact_slice_ph, &psi_small_d[0 + batch_idx*psi_small_size],
			                           psi_small_size * sizeof(PRISM_CUDA_COMPLEX_FLOAT),
			                           cudaMemcpyDeviceToHost, stream));
			cudaStreamSynchronize(stream);
			memcpy(&pars.Scompact[(beamNumber + batch_idx) * pars.Scompact.get_dimj() * pars.Scompact.get_dimi()],
			       &Scompact_slice_ph[0], psi_small_size * sizeof(PRISM_CUDA_COMPLEX_FLOAT));
		}
	}

	void propagatePlaneWave_GPU_streaming(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                                      PRISM_CUDA_COMPLEX_FLOAT* trans_d,
	                                      const std::complex<PRISM_FLOAT_PRECISION> *trans_ph,
	                                      PRISM_CUDA_COMPLEX_FLOAT* psi_d,
	                                      PRISM_CUDA_COMPLEX_FLOAT* psi_small_d,
	                                      complex<PRISM_FLOAT_PRECISION>* Scompact_slice_ph,
	                                      const size_t* qyInd_d,
	                                      const size_t* qxInd_d,
	                                      const PRISM_CUDA_COMPLEX_FLOAT* prop_d,
	                                      const size_t* beamsIndex,
	                                      const size_t beamNumber,
	                                      const cufftHandle& plan,
	                                      const cufftHandle& plan_small,
	                                      cudaStream_t& stream){
		// In this version, each slice of the transmission matrix is streamed to the device

		const size_t psi_size = pars.imageSize[0] * pars.imageSize[1];
		const size_t psi_small_size = pars.qxInd.size() * pars.qyInd.size();
		initializePsi_oneNonzero<<< (psi_size - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream>>>(psi_d, psi_size, pars.beamsIndex[beamNumber]);

		for (auto planeNum = 0; planeNum < pars.numPlanes ; ++planeNum) {
			cudaErrchk(cudaMemcpyAsync(trans_d, &trans_ph[planeNum*psi_size], psi_size * sizeof(PRISM_CUDA_COMPLEX_FLOAT), cudaMemcpyHostToDevice, stream));
			cufftErrchk(PRISM_CUFFT_EXECUTE(plan, &psi_d[0], &psi_d[0], CUFFT_INVERSE));
			multiply_cx<<<(psi_size - 1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_d, trans_d, psi_size);
			divide_inplace<<<(psi_size - 1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_d, PRISM_MAKE_CU_COMPLEX(psi_size, 0), psi_size);
			cufftErrchk(PRISM_CUFFT_EXECUTE(plan, &psi_d[0], &psi_d[0], CUFFT_FORWARD));
			multiply_cx<<<(psi_size - 1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_d, prop_d, psi_size);
		}
		array_subset<<<(pars.qyInd.size()*pars.qxInd.size()-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>> (
				psi_d, psi_small_d, qyInd_d, qxInd_d, pars.imageSize[1], pars.qyInd.size(), pars.qxInd.size());

		PRISM_CUFFT_EXECUTE(plan_small,&psi_small_d[0], &psi_small_d[0], CUFFT_INVERSE);
		divide_inplace<<<(psi_small_size-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_small_d, PRISM_MAKE_CU_COMPLEX(psi_small_size, 0),psi_small_size);

		cudaErrchk(cudaMemcpyAsync(Scompact_slice_ph,&psi_small_d[0],psi_small_size * sizeof(PRISM_CUDA_COMPLEX_FLOAT),cudaMemcpyDeviceToHost,stream));
		cudaStreamSynchronize(stream);
		memcpy(&pars.Scompact[beamNumber * pars.Scompact.get_dimj() * pars.Scompact.get_dimi()], &Scompact_slice_ph[0], psi_small_size * sizeof(PRISM_CUDA_COMPLEX_FLOAT));
	}


	void propagatePlaneWave_GPU_streaming_batch(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                                            PRISM_CUDA_COMPLEX_FLOAT* trans_d,
	                                            const std::complex<PRISM_FLOAT_PRECISION> *trans_ph,
	                                            PRISM_CUDA_COMPLEX_FLOAT* psi_d,
	                                            PRISM_CUDA_COMPLEX_FLOAT* psi_small_d,
	                                            complex<PRISM_FLOAT_PRECISION>* Scompact_slice_ph,
	                                            const size_t* qyInd_d,
	                                            const size_t* qxInd_d,
	                                            const PRISM_CUDA_COMPLEX_FLOAT* prop_d,
	                                            const size_t* beamsIndex,
	                                            const size_t beamNumber,
	                                            const size_t stopBeam,
	                                            const cufftHandle& plan,
	                                            const cufftHandle& plan_small,
	                                            cudaStream_t& stream){
		// In this version, each slice of the transmission matrix is streamed to the device

		const size_t psi_size        = pars.imageSize[0] * pars.imageSize[1];
		const size_t psi_small_size = pars.qxInd.size() * pars.qyInd.size();
		for (auto batch_idx = 0; batch_idx < (stopBeam-beamNumber); ++batch_idx) {
			// initialize psi -- for PRISM this is just a delta function in Fourier space located depending on which plane wave it is
			initializePsi_oneNonzero<<< (psi_size - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream>>>(psi_d + batch_idx*psi_size, psi_size, pars.beamsIndex[beamNumber + batch_idx]);
		}

		for (auto planeNum = 0; planeNum < pars.numPlanes ; ++planeNum) {
			cudaErrchk(cudaMemcpyAsync(trans_d, &trans_ph[planeNum*psi_size], psi_size * sizeof(PRISM_CUDA_COMPLEX_FLOAT), cudaMemcpyHostToDevice, stream));
			cufftErrchk(PRISM_CUFFT_EXECUTE(plan, &psi_d[0], &psi_d[0], CUFFT_INVERSE));
			for (auto batch_idx = 0; batch_idx < (stopBeam-beamNumber); ++batch_idx) {
				multiply_cx << < (psi_size - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >> >
						(psi_d + batch_idx*psi_size, trans_d, psi_size); // transmit
				divide_inplace << < (psi_size - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >> >
						(psi_d + batch_idx*psi_size, PRISM_MAKE_CU_COMPLEX(psi_size, 0), psi_size); // normalize the FFT
			}
			cufftErrchk(PRISM_CUFFT_EXECUTE(plan, &psi_d[0], &psi_d[0], CUFFT_FORWARD));
			for (auto batch_idx = 0; batch_idx < (stopBeam-beamNumber); ++batch_idx) {
				multiply_cx << < (psi_size - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >> >
						(psi_d + batch_idx*psi_size, prop_d, psi_size); // propagate
			}
		}

		for (auto batch_idx = 0; batch_idx < (stopBeam-beamNumber); ++batch_idx) {
			// take relevant subset of the full array
			array_subset << < (pars.qyInd.size() * pars.qxInd.size() - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0,
					stream >> > (psi_d + batch_idx*psi_size, psi_small_d + batch_idx*psi_small_size, qyInd_d, qxInd_d, pars.imageSize[1], pars.qyInd.size(), pars.qxInd.size());
		}

		// final FFT
		PRISM_CUFFT_EXECUTE(plan_small,&psi_small_d[0], &psi_small_d[0], CUFFT_INVERSE);
		for (auto batch_idx = 0; batch_idx < (stopBeam-beamNumber); ++batch_idx) {
		divide_inplace<<<(psi_small_size-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>
				(psi_small_d + batch_idx*psi_small_size, PRISM_MAKE_CU_COMPLEX(psi_small_size, 0),psi_small_size); // normalize the FFT
			}

		// copy the result
		for (auto batch_idx = 0; batch_idx < (stopBeam-beamNumber); ++batch_idx) {
		cudaErrchk(cudaMemcpyAsync(Scompact_slice_ph,&psi_small_d[batch_idx*psi_small_size],psi_small_size * sizeof(PRISM_CUDA_COMPLEX_FLOAT),cudaMemcpyDeviceToHost,stream));
		cudaStreamSynchronize(stream);
		memcpy(&pars.Scompact[beamNumber * pars.Scompact.get_dimj() * pars.Scompact.get_dimi()], &Scompact_slice_ph[0], psi_small_size * sizeof(PRISM_CUDA_COMPLEX_FLOAT));
			}
	}

	inline void setupArrays2(Parameters<PRISM_FLOAT_PRECISION>& pars){

		// setup some needed arrays
		const PRISM_FLOAT_PRECISION pi = acos(-1);
		const std::complex<PRISM_FLOAT_PRECISION> i(0, 1);
		pars.Scompact = zeros_ND<3, complex<PRISM_FLOAT_PRECISION> >(
				{{pars.numberBeams, pars.imageSize[0] / 2, pars.imageSize[1] / 2}});
		pars.transmission = zeros_ND<3, complex<PRISM_FLOAT_PRECISION> >(
				{{pars.pot.get_dimk(), pars.pot.get_dimj(), pars.pot.get_dimi()}});
		{
			auto p = pars.pot.begin();
			for (auto &j:pars.transmission)j = exp(i * pars.sigma * (*p++));
		}
	}

	void fill_Scompact_GPU_singlexfer(Parameters <PRISM_FLOAT_PRECISION> &pars) {

		// This version transfers the entire transmission matrix a single time, which results in faster execution but requires more memory
#ifdef PRISM_BUILDING_GUI
		pars.progressbar->signalDescriptionMessage("Computing compact S-matrix");
		pars.progressbar->signalScompactUpdate(-1, pars.numberBeams);
#endif
		CudaParameters<PRISM_FLOAT_PRECISION> cuda_pars;

		// determine the batch size to use
        pars.meta.batch_size_GPU = min(pars.meta.batch_size_target_GPU, max((size_t)1, pars.numberBeams / max((size_t)1,(pars.meta.NUM_STREAMS_PER_GPU*pars.meta.NUM_GPUS))));

		// setup some arrays
		setupArrays2(pars);

		// create CUDA streams
		createStreamsAndPlans2(pars, cuda_pars);

		// create page-locked (pinned) host memory buffers
		allocatePinnedHostMemory_singlexfer2(pars, cuda_pars);

		// copy to pinned memory
		copyToPinnedMemory_singlexfer2(pars, cuda_pars);

		// allocate memory on the GPUs
		allocateDeviceMemory_singlexfer2(pars, cuda_pars);

		// copy to GPUs
		copyToDeviceMemory_singlexfer2(pars, cuda_pars);

		// launch workers
		launchWorkers_singlexfer2(pars, cuda_pars);

		// free memory on the host/device
		cleanupMemory2(pars, cuda_pars);
	}

	void fill_Scompact_GPU_streaming(Parameters <PRISM_FLOAT_PRECISION> &pars) {

#ifdef PRISM_BUILDING_GUI
		pars.progressbar->signalDescriptionMessage("Computing compact S-matrix");
		pars.progressbar->signalScompactUpdate(-1, pars.numberBeams);
#endif
		// This version streams each slice of the transmission matrix, which is less efficient but can tolerate very large arrays
		//initialize data
		CudaParameters<PRISM_FLOAT_PRECISION> cuda_pars;

		// determine the batch size to use
		pars.meta.batch_size_GPU = min(pars.meta.batch_size_target_GPU, max((size_t)1, pars.numberBeams / max((size_t)1,(pars.meta.NUM_STREAMS_PER_GPU*pars.meta.NUM_GPUS))));

		// setup some arrays
		setupArrays2(pars);

		// create CUDA streams and cuFFT plans
		createStreamsAndPlans2(pars, cuda_pars);

		// create page-locked (pinned) host memory buffers
		allocatePinnedHostMemory_streaming2(pars, cuda_pars);

		// copy to pinned memory
		copyToPinnedMemory_streaming2(pars, cuda_pars);

		// allocate memory on the GPUs
		allocateDeviceMemory_streaming2(pars, cuda_pars);

		// copy to GPUs
		copyToDeviceMemory_streaming2(pars, cuda_pars);

		// launch workers
		launchWorkers_streaming(pars, cuda_pars);

		// free memory on the host/device
		cleanupMemory2(pars, cuda_pars);
	}
}