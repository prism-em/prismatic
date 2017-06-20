// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

// / Calculate result of Multislice simulation using GPU and (potentially) CPU. Multiple GPU threads are launched, each with
// their own memory buffers. Page-locked host memory is allocated so that memory transfers to the GPU can occur asynchronously,
// and memory allocation for the GPU occurs only once, as each call to cudaMalloc will potentially interrupt concurrent execution.
// Each GPU/CPU worker thread repeatedly calls getWorkID to be assigned probe positions to compute. This queue mechanism
// ensures that both the CPU and GPU are kept busy.

// For variable naming, the suffixes are "_d" for "device" (1 copy per GPU), "_ds" for "device stream (1 copy per stream), "_ph" for "pinned host"

#include "Multislice_calcOutput.cuh"
#include "Multislice_calcOutput.h"
#include "cuComplex.h"
#include "cufft.h"
#include "utility.cuh"
#include "params.cuh"

namespace Prismatic{
	extern std::mutex fftw_plan_lock;

	inline void createPlansAndStreamsM(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
									   CudaParameters<PRISMATIC_FLOAT_PRECISION> &cuda_pars){
		// create CUDA streams
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;
		cuda_pars.streams   = new cudaStream_t[total_num_streams];
		cuda_pars.cufft_plans = new cufftHandle[total_num_streams];

		// batch parameters for cuFFT
		const int rank = 2;
		int n[] = {(int)pars.psiProbeInit.get_dimj(), (int)pars.psiProbeInit.get_dimi()};
		const int howmany = pars.meta.batch_size_GPU;
		int idist = n[0]*n[1];
		int odist = n[0]*n[1];
		int istride = 1;
		int ostride = 1;
		int *inembed = n;
		int *onembed = n;

		for (auto j = 0; j < total_num_streams; ++j){
			cudaSetDevice(j % pars.meta.NUM_GPUS);
			cudaErrchk(cudaStreamCreate(&cuda_pars.streams[j]));
			cufftErrchk(cufftPlanMany(&cuda_pars.cufft_plans[j], rank, n, inembed, istride, idist, onembed, ostride, odist, PRISMATIC_CUFFT_PLAN_TYPE, howmany));
			cufftErrchk(cufftSetStream(cuda_pars.cufft_plans[j], cuda_pars.streams[j]));
		}
	}

	inline void allocatePinnedHostMemory_M(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
	                                       CudaParameters<PRISMATIC_FLOAT_PRECISION> &cuda_pars){
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;
		cuda_pars.output_ph = new PRISMATIC_FLOAT_PRECISION*[total_num_streams];
		// allocate pinned memory
		cudaErrchk(cudaMallocHost((void **)&cuda_pars.PsiProbeInit_ph, pars.psiProbeInit.size()*sizeof(std::complex<PRISMATIC_FLOAT_PRECISION>)));
		cudaErrchk(cudaMallocHost((void **)&cuda_pars.trans_ph,        pars.transmission.size()*sizeof(std::complex<PRISMATIC_FLOAT_PRECISION>)));
		cudaErrchk(cudaMallocHost((void **)&cuda_pars.prop_ph,         pars.prop.size()*sizeof(std::complex<PRISMATIC_FLOAT_PRECISION>)));
		cudaErrchk(cudaMallocHost((void **)&cuda_pars.qxa_ph,          pars.qxa.size()*sizeof(PRISMATIC_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **)&cuda_pars.qya_ph,          pars.qya.size()*sizeof(PRISMATIC_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **)&cuda_pars.alphaInd_ph,     pars.alphaInd.size()*sizeof(PRISMATIC_FLOAT_PRECISION)));
		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaMallocHost((void **) &cuda_pars.output_ph[s], pars.output.get_dimi() * sizeof(PRISMATIC_FLOAT_PRECISION)));
		}
	}

	inline void copyToPinnedMemory_M(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
	                                          CudaParameters<PRISMATIC_FLOAT_PRECISION> &cuda_pars){
		// copy host memory to pinned
		memcpy(cuda_pars.PsiProbeInit_ph, &pars.psiProbeInit[0], pars.psiProbeInit.size() * sizeof(std::complex<PRISMATIC_FLOAT_PRECISION>));
		memcpy(cuda_pars.trans_ph,        &pars.transmission[0], pars.transmission.size() * sizeof(std::complex<PRISMATIC_FLOAT_PRECISION>));
		memcpy(cuda_pars.prop_ph,         &pars.prop[0],         pars.prop.size()         * sizeof(std::complex<PRISMATIC_FLOAT_PRECISION>));
		memcpy(cuda_pars.qxa_ph,          &pars.qxa[0],          pars.qxa.size()          * sizeof(PRISMATIC_FLOAT_PRECISION));
		memcpy(cuda_pars.qya_ph,          &pars.qya[0],          pars.qya.size()          * sizeof(PRISMATIC_FLOAT_PRECISION));
		memcpy(cuda_pars.alphaInd_ph,     &pars.alphaInd[0],     pars.alphaInd.size()     * sizeof(PRISMATIC_FLOAT_PRECISION));
	}

	inline void allocateDeviceMemory_singlexferM(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
	                                             CudaParameters<PRISMATIC_FLOAT_PRECISION> &cuda_pars){
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;

		// pointers to read-only GPU memory (one copy per GPU)
		cuda_pars.PsiProbeInit_d  = new PRISMATIC_CUDA_COMPLEX_FLOAT*[pars.meta.NUM_GPUS];
		cuda_pars.trans_d		  = new PRISMATIC_CUDA_COMPLEX_FLOAT*[pars.meta.NUM_GPUS];
		cuda_pars.prop_d 		  = new PRISMATIC_CUDA_COMPLEX_FLOAT*[pars.meta.NUM_GPUS];
		cuda_pars.qxa_d 		  = new PRISMATIC_FLOAT_PRECISION*[pars.meta.NUM_GPUS];
		cuda_pars.qya_d 		  = new PRISMATIC_FLOAT_PRECISION*[pars.meta.NUM_GPUS];
		cuda_pars.alphaInd_d      = new PRISMATIC_FLOAT_PRECISION*[pars.meta.NUM_GPUS];

		// pointers to read/write GPU memory (one per stream)
		cuda_pars.psi_ds 			  = new PRISMATIC_CUDA_COMPLEX_FLOAT*[total_num_streams];
		cuda_pars.psi_intensity_ds    = new PRISMATIC_FLOAT_PRECISION*[total_num_streams];
		cuda_pars.integratedOutput_ds = new PRISMATIC_FLOAT_PRECISION*[total_num_streams];

		// allocate memory on each GPU
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.PsiProbeInit_d[g],     pars.psiProbeInit.size()   * sizeof(std::complex<PRISMATIC_FLOAT_PRECISION>)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.trans_d[g],            pars.transmission.size()   * sizeof(std::complex<PRISMATIC_FLOAT_PRECISION>)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.prop_d[g],             pars.prop.size()           * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.qxa_d[g],              pars.qxa.size()            * sizeof(PRISMATIC_FLOAT_PRECISION)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.qya_d[g],              pars.qya.size()            * sizeof(PRISMATIC_FLOAT_PRECISION)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.alphaInd_d[g],         pars.alphaInd.size()       * sizeof(PRISMATIC_FLOAT_PRECISION)));
		}

		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaSetDevice(s % pars.meta.NUM_GPUS));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.psi_ds[s],              pars.meta.batch_size_GPU*pars.psiProbeInit.size() * sizeof(std::complex<PRISMATIC_FLOAT_PRECISION>)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.psi_intensity_ds[s],    pars.meta.batch_size_GPU*pars.psiProbeInit.size() * sizeof(PRISMATIC_FLOAT_PRECISION)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.integratedOutput_ds[s], pars.detectorAngles.size()                        * sizeof(PRISMATIC_FLOAT_PRECISION)));
			cudaErrchk(cudaMemset(cuda_pars.psi_ds[s], 0,                      pars.meta.batch_size_GPU*pars.psiProbeInit.size() * sizeof(std::complex<PRISMATIC_FLOAT_PRECISION>)));
			cudaErrchk(cudaMemset(cuda_pars.psi_intensity_ds[s], 0,            pars.meta.batch_size_GPU*pars.psiProbeInit.size() * sizeof(PRISMATIC_FLOAT_PRECISION)));
			cudaErrchk(cudaMemset(cuda_pars.integratedOutput_ds[s], 0,         pars.detectorAngles.size()                        * sizeof(PRISMATIC_FLOAT_PRECISION)));
		}
	}

	inline void allocateDeviceMemory_streamingM(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
	                                            CudaParameters<PRISMATIC_FLOAT_PRECISION> &cuda_pars){
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;

		// pointers to read-only GPU memory (one copy per GPU)
		cuda_pars.PsiProbeInit_d  = new PRISMATIC_CUDA_COMPLEX_FLOAT*[pars.meta.NUM_GPUS];
		cuda_pars.prop_d 	   	  = new PRISMATIC_CUDA_COMPLEX_FLOAT*[pars.meta.NUM_GPUS];
		cuda_pars.qxa_d 		  = new PRISMATIC_FLOAT_PRECISION*[pars.meta.NUM_GPUS];
		cuda_pars.qya_d 		  = new PRISMATIC_FLOAT_PRECISION*[pars.meta.NUM_GPUS];
		cuda_pars.alphaInd_d 	  = new PRISMATIC_FLOAT_PRECISION*[pars.meta.NUM_GPUS];

		// pointers to read/write GPU memory (one per stream)
		cuda_pars.trans_d 		      = new PRISMATIC_CUDA_COMPLEX_FLOAT*[total_num_streams];
		cuda_pars.psi_ds  		      = new PRISMATIC_CUDA_COMPLEX_FLOAT*[total_num_streams];
		cuda_pars.psi_intensity_ds    = new PRISMATIC_FLOAT_PRECISION*[total_num_streams];
		cuda_pars.integratedOutput_ds = new PRISMATIC_FLOAT_PRECISION*[total_num_streams];

		// allocate memory on each GPU
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.PsiProbeInit_d[g],     pars.psiProbeInit.size() * sizeof(pars.psiProbeInit[0])));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.prop_d[g],             pars.prop.size()         * sizeof(pars.prop[0])));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.qxa_d[g],              pars.qxa.size()          * sizeof(pars.qxa[0])));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.qya_d[g],              pars.qya.size()          * sizeof(pars.qya[0])));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.alphaInd_d[g],         pars.alphaInd.size()     * sizeof(pars.alphaInd[0])));
		}

		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaSetDevice(s % pars.meta.NUM_GPUS));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.trans_d[s],             pars.transmission.get_dimj() * pars.transmission.get_dimi() * sizeof(pars.transmission[0])));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.psi_ds[s],              pars.meta.batch_size_GPU*pars.psiProbeInit.size()           * sizeof(pars.psiProbeInit[0])));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.psi_intensity_ds[s],    pars.meta.batch_size_GPU*pars.psiProbeInit.size()           * sizeof(PRISMATIC_FLOAT_PRECISION)));
			cudaErrchk(cudaMalloc((void **) &cuda_pars.integratedOutput_ds[s], pars.detectorAngles.size()                                  * sizeof(PRISMATIC_FLOAT_PRECISION)));

			cudaErrchk(cudaMemset(cuda_pars.psi_ds[s],              0, pars.meta.batch_size_GPU*pars.psiProbeInit.size() * sizeof(pars.psiProbeInit[0])));
			cudaErrchk(cudaMemset(cuda_pars.psi_intensity_ds[s],    0, pars.meta.batch_size_GPU*pars.psiProbeInit.size() * sizeof(PRISMATIC_FLOAT_PRECISION)));
			cudaErrchk(cudaMemset(cuda_pars.integratedOutput_ds[s], 0, pars.detectorAngles.size()                        * sizeof(PRISMATIC_FLOAT_PRECISION)));
		}
	}

	inline void copyToGPUMemory_singlexferM(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
	                                        CudaParameters<PRISMATIC_FLOAT_PRECISION> &cuda_pars){

		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;

		// Copy memory to each GPU asynchronously from the pinned host memory spaces.
		// The streams are laid out so that consecutive streams represent different GPUs. If we
		// have more than one stream per GPU, then we want to interleave as much as possible
		int stream_id = 0;
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			stream_id = g;
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaMemcpyAsync(cuda_pars.PsiProbeInit_d[g], &cuda_pars.PsiProbeInit_ph[0], pars.psiProbeInit.size() * sizeof(pars.psiProbeInit[0]), cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));

			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(cuda_pars.trans_d[g], &cuda_pars.trans_ph[0], pars.transmission.size() * sizeof(pars.transmission[0]), cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));

			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(cuda_pars.prop_d[g], &cuda_pars.prop_ph[0], pars.prop.size() * sizeof(pars.prop[0]), cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));

			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(cuda_pars.qxa_d[g], &cuda_pars.qxa_ph[0], pars.qxa.size() * sizeof(pars.qxa[0]), cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));

			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(cuda_pars.qya_d[g], &cuda_pars.qya_ph[0], pars.qya.size() * sizeof(pars.qya[0]), cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));

			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(cuda_pars.alphaInd_d[g], &cuda_pars.alphaInd_ph[0], pars.alphaInd.size() * sizeof(pars.alphaInd[0]), cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));
		}

		// make sure transfers are complete
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g){
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaDeviceSynchronize());
		}
	}

	inline void copyToGPUMemory_streamingM(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
	                                       CudaParameters<PRISMATIC_FLOAT_PRECISION> &cuda_pars){
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;

		// Copy memory to each GPU asynchronously from the pinned host memory spaces.
		// The streams are laid out so that consecutive streams represent different GPUs. If we
		// have more than one stream per GPU, then we want to interleave as much as possible
		int stream_id = 0;
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			stream_id = g;
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaMemcpyAsync(cuda_pars.PsiProbeInit_d[g], &cuda_pars.PsiProbeInit_ph[0], pars.psiProbeInit.size() * sizeof(pars.psiProbeInit[0]), cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));

			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(cuda_pars.prop_d[g], &cuda_pars.prop_ph[0], pars.prop.size() * sizeof(pars.prop[0]), cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));

			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(cuda_pars.qxa_d[g], &cuda_pars.qxa_ph[0], pars.qxa.size() * sizeof(pars.qxa[0]), cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));

			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(cuda_pars.qya_d[g], &cuda_pars.qya_ph[0], pars.qya.size() * sizeof(pars.qya[0]), cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));

			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(cuda_pars.alphaInd_d[g], &cuda_pars.alphaInd_ph[0], pars.alphaInd.size() * sizeof(pars.alphaInd[0]), cudaMemcpyHostToDevice, cuda_pars.streams[stream_id]));
		}

		// make sure transfers are complete
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g){
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaDeviceSynchronize());
		}
	}

	inline void launchWorkers_singlexferM(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
	                                      CudaParameters<PRISMATIC_FLOAT_PRECISION> &cuda_pars){
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;

		vector<thread> workers_GPU;
		workers_GPU.reserve(total_num_streams); // prevents multiple reallocations
		size_t psi_size = pars.psiProbeInit.size();
		int stream_count = 0;
		const size_t PRISMATIC_PRINT_FREQUENCY_PROBES = max((size_t)1,pars.xp.size() * pars.yp.size() / 10); // for printing status
		WorkDispatcher dispatcher(0, pars.xp.size() * pars.yp.size());

		for (auto t = 0; t < total_num_streams; ++t){
			int GPU_num = stream_count % pars.meta.NUM_GPUS; // determine which GPU handles this job
			cudaStream_t& current_stream = cuda_pars.streams[stream_count];
			cout << "Launching GPU worker on stream #" << stream_count << " on GPU #" << GPU_num << '\n';

			// get pointers to the pre-copied arrays, making sure to get those on the current GPU
			PRISMATIC_CUDA_COMPLEX_FLOAT *current_PsiProbeInit_d = cuda_pars.PsiProbeInit_d[GPU_num];
			PRISMATIC_CUDA_COMPLEX_FLOAT *current_trans_d        = cuda_pars.trans_d[GPU_num];
			PRISMATIC_CUDA_COMPLEX_FLOAT *current_prop_d         = cuda_pars.prop_d[GPU_num];
			PRISMATIC_FLOAT_PRECISION *current_qxa_d             = cuda_pars.qxa_d[GPU_num];
			PRISMATIC_FLOAT_PRECISION *current_qya_d             = cuda_pars.qya_d[GPU_num];
			PRISMATIC_FLOAT_PRECISION *current_alphaInd_d        = cuda_pars.alphaInd_d[GPU_num];

			// get pointers to per-stream arrays
			PRISMATIC_CUDA_COMPLEX_FLOAT *current_psi_ds           = cuda_pars.psi_ds[stream_count];
			PRISMATIC_FLOAT_PRECISION *current_psi_intensity_ds    = cuda_pars.psi_intensity_ds[stream_count];
			PRISMATIC_FLOAT_PRECISION *current_integratedOutput_ds = cuda_pars.integratedOutput_ds[stream_count];
			PRISMATIC_FLOAT_PRECISION *current_output_ph           = cuda_pars.output_ph[stream_count];
			cufftHandle & current_cufft_plan                   = cuda_pars.cufft_plans[stream_count];

			// launch a new thread
			workers_GPU.push_back(thread([&pars, current_trans_d, current_PsiProbeInit_d, current_alphaInd_d, &dispatcher,
					                             current_psi_ds, current_psi_intensity_ds, current_integratedOutput_ds,
					                             GPU_num, current_qya_d, current_qxa_d, current_output_ph, &current_cufft_plan,
					                             current_prop_d, &current_stream, &psi_size, stream_count, &PRISMATIC_PRINT_FREQUENCY_PROBES, &cuda_pars]() {

				// set the GPU context
				cudaErrchk(cudaSetDevice(GPU_num)); // set current GPU

#ifndef NDEBUG
				{
//					 check memory usage on the GPU
					std::lock_guard<mutex> lock(Prismatic::mem_lock);
					size_t free_mem, total_mem;
					free_mem=total_mem=0;
					cudaErrchk(cudaMemGetInfo(&free_mem, &total_mem));
					pars.max_mem = std::max(total_mem - free_mem, pars.max_mem);
				}
#endif // NDEBUG
				size_t Nstart, Nstop;
				Nstart=Nstop=0;
				while (dispatcher.getWork(Nstart, Nstop, pars.meta.batch_size_GPU)){ // synchronously get work assignment
					while (Nstart < Nstop){
						if (Nstart % PRISMATIC_PRINT_FREQUENCY_PROBES < pars.meta.batch_size_GPU | Nstart == 100){
							cout << "Computing Probe Position #" << Nstart << "/" << pars.xp.size() * pars.yp.size() << '\n';
						}
//							getMultisliceProbe_GPU_singlexfer(pars, current_trans_d, current_PsiProbeInit_d, current_psi_ds, current_output_ph,
//							                                  current_psi_intensity_ds,
//							                                  current_integratedOutput_ds, current_qya_d, current_qxa_d,
//							                                  current_prop_d, ay, ax, pars.psiProbeInit.get_dimj(), pars.psiProbeInit.get_dimi(),
//							                                  current_alphaInd_d, current_cufft_plan, current_stream);
						getMultisliceProbe_GPU_singlexfer_batch(pars, current_trans_d, current_PsiProbeInit_d, current_psi_ds, current_output_ph,
						                                        current_psi_intensity_ds,
						                                        current_integratedOutput_ds, current_qya_d, current_qxa_d,
						                                        current_prop_d, Nstart, Nstop, pars.psiProbeInit.get_dimj(), pars.psiProbeInit.get_dimi(),
						                                        current_alphaInd_d, current_cufft_plan, current_stream);
#ifdef PRISMATIC_BUILDING_GUI
						pars.progressbar->signalOutputUpdate(Nstart, pars.xp.size() * pars.yp.size());
#endif
						Nstart=Nstop;
					}
				}
				cout << "GPU worker on stream #" << stream_count << " of GPU #" << GPU_num << "finished\n";
			}));
			++stream_count;
		}
		// now launch CPU work
		if (pars.meta.also_do_CPU_work){
			PRISMATIC_FFTW_INIT_THREADS();
			PRISMATIC_FFTW_PLAN_WITH_NTHREADS(pars.meta.NUM_THREADS);vector<thread> workers_CPU;
			workers_CPU.reserve(pars.meta.NUM_THREADS); // prevents multiple reallocations
			// If the batch size is too big, the work won't be spread over the threads, which will usually hurt more than the benefit
			// of batch FFT
			pars.meta.batch_size_CPU = min(pars.meta.batch_size_target_CPU, max((size_t)1, pars.xp.size() * pars.yp.size() / pars.meta.NUM_THREADS));
			cout << "multislice pars.meta.batch_size_GPU = " << pars.meta.batch_size_GPU << endl;
			for (auto t = 0; t < pars.meta.NUM_THREADS; ++t) {
				cout << "Launching CPU worker #" << t << endl;
				// push_back is better whenever constructing a new object
				workers_CPU.push_back(thread([&pars, &dispatcher, t, &PRISMATIC_PRINT_FREQUENCY_PROBES]() {
					size_t Nstart, Nstop, early_CPU_stop;
					Nstart=Nstop=0;
					// stop the CPU workers earlier than the GPU ones to prevent slower workers taking the last jobs and having to
					// wait longer for everything to complete
					if (pars.meta.NUM_GPUS > 0){
						// if there are no GPUs, make sure to do all work on CPU
						early_CPU_stop = (size_t)std::max((PRISMATIC_FLOAT_PRECISION)0.0, pars.xp.size() * pars.yp.size() - pars.meta.gpu_cpu_ratio);
					} else {
						early_CPU_stop = pars.xp.size() * pars.yp.size();
					}
					if (dispatcher.getWork(Nstart, Nstop, pars.meta.batch_size_CPU, early_CPU_stop)) { // synchronously get work assignment
						Array1D<std::complex<PRISMATIC_FLOAT_PRECISION> > psi_stack = zeros_ND<1, complex<PRISMATIC_FLOAT_PRECISION> >({{pars.psiProbeInit.size() * pars.meta.batch_size_CPU}});

						// setup batch FFTW parameters
						const int rank    = 2;
						int n[]           = {(int)pars.psiProbeInit.get_dimj(), (int)pars.psiProbeInit.get_dimi()};
						const int howmany = pars.meta.batch_size_CPU;
						int idist         = n[0]*n[1];
						int odist         = n[0]*n[1];
						int istride       = 1;
						int ostride       = 1;
						int *inembed      = n;
						int *onembed      = n;
						unique_lock<mutex> gatekeeper(fftw_plan_lock);
						PRISMATIC_FFTW_PLAN plan_forward = PRISMATIC_FFTW_PLAN_DFT_BATCH(rank, n, howmany,
						                                                         reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&psi_stack[0]), inembed,
						                                                         istride, idist,
						                                                         reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&psi_stack[0]), onembed,
						                                                         ostride, odist,
						                                                         FFTW_FORWARD, FFTW_MEASURE);
						PRISMATIC_FFTW_PLAN plan_inverse = PRISMATIC_FFTW_PLAN_DFT_BATCH(rank, n, howmany,
						                                                         reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&psi_stack[0]), inembed,
						                                                         istride, idist,
						                                                         reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&psi_stack[0]), onembed,
						                                                         ostride, odist,
						                                                         FFTW_BACKWARD, FFTW_MEASURE);
						gatekeeper.unlock();

						// main work loop
						do {
							while (Nstart < Nstop) {
								if (Nstart % PRISMATIC_PRINT_FREQUENCY_PROBES  < pars.meta.batch_size_CPU | Nstart == 100){
									cout << "Computing Probe Position #" << Nstart << "/" << pars.xp.size() * pars.yp.size() << endl;
								}
//							getMultisliceProbe_CPU(pars, ay, ax, plan_forward, plan_inverse, psi);
								getMultisliceProbe_CPU_batch(pars, Nstart, Nstop, plan_forward, plan_inverse, psi_stack);
#ifdef PRISMATIC_BUILDING_GUI
								pars.progressbar->signalOutputUpdate(Nstart, pars.xp.size() * pars.yp.size());
#endif
								Nstart=Nstop;
							}
							if (Nstop >= early_CPU_stop) break;
						} while(dispatcher.getWork(Nstart, Nstop, pars.meta.batch_size_CPU, early_CPU_stop));
						gatekeeper.lock();
						PRISMATIC_FFTW_DESTROY_PLAN(plan_forward);
						PRISMATIC_FFTW_DESTROY_PLAN(plan_inverse);
						gatekeeper.unlock();
					}
				}));
			}
			cout << "Waiting on CPU threads..." << endl;
			for (auto& t:workers_CPU)t.join();
			PRISMATIC_FFTW_CLEANUP_THREADS();
		}
		// synchronize threads
		cout << "Waiting on GPU threads..." << endl;
		for (auto& t:workers_GPU)t.join();

		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g){
			cudaSetDevice(g);
			cudaDeviceSynchronize();
		}
	}

	inline void launchWorkers_streamingM(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
	                                     CudaParameters<PRISMATIC_FLOAT_PRECISION> &cuda_pars){
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;
		vector<thread> workers_GPU;
		workers_GPU.reserve(total_num_streams); // prevents multiple reallocations
		size_t psi_size = pars.psiProbeInit.size();
		int stream_count = 0;
		const size_t PRISMATIC_PRINT_FREQUENCY_PROBES = max((size_t)1,pars.xp.size() * pars.yp.size() / 10); // for printing status
		WorkDispatcher dispatcher(0, pars.xp.size() * pars.yp.size());
		// If the batch size is too big, the work won't be spread over the threads, which will usually hurt more than the benefit
		// of batch FFT

		for (auto t = 0; t < total_num_streams; ++t){
			int GPU_num = stream_count % pars.meta.NUM_GPUS; // determine which GPU handles this job
			cudaStream_t& current_stream = cuda_pars.streams[stream_count];
			cout << "Launching GPU worker on stream #" << stream_count << " of GPU #" << GPU_num << endl;

			// get pointers to the pre-copied arrays, making sure to get those on the current GPU
			PRISMATIC_CUDA_COMPLEX_FLOAT *current_PsiProbeInit_d   = cuda_pars.PsiProbeInit_d[GPU_num];
			PRISMATIC_CUDA_COMPLEX_FLOAT *current_prop_d           = cuda_pars.prop_d[GPU_num];
			PRISMATIC_FLOAT_PRECISION *current_qxa_d               = cuda_pars.qxa_d[GPU_num];
			PRISMATIC_FLOAT_PRECISION *current_qya_d               = cuda_pars.qya_d[GPU_num];
			PRISMATIC_FLOAT_PRECISION *current_alphaInd_d          = cuda_pars.alphaInd_d[GPU_num];

			// get pointers to per-stream arrays
			PRISMATIC_CUDA_COMPLEX_FLOAT *current_trans_ds         = cuda_pars.trans_d[stream_count];
			PRISMATIC_CUDA_COMPLEX_FLOAT *current_psi_ds           = cuda_pars.psi_ds[stream_count];
			PRISMATIC_FLOAT_PRECISION *current_psi_intensity_ds    = cuda_pars.psi_intensity_ds[stream_count];
			PRISMATIC_FLOAT_PRECISION *current_integratedOutput_ds = cuda_pars.integratedOutput_ds[stream_count];
			PRISMATIC_FLOAT_PRECISION *current_output_ph           = cuda_pars.output_ph[stream_count];
			cufftHandle & current_cufft_plan                   = cuda_pars.cufft_plans[stream_count];
			// launch a new thread
			// push_back is better whenever constructing a new object
			workers_GPU.push_back(thread([&pars, current_trans_ds, current_PsiProbeInit_d, current_alphaInd_d, &dispatcher,
					                             current_psi_ds, current_psi_intensity_ds, current_integratedOutput_ds,
					                             GPU_num, current_qya_d, current_qxa_d, current_output_ph, current_cufft_plan,
					                             current_prop_d, &current_stream, &psi_size, stream_count, &PRISMATIC_PRINT_FREQUENCY_PROBES, &cuda_pars]()  {

				// set the GPU context
				cudaErrchk(cudaSetDevice(GPU_num)); // set current GPU


#ifndef NDEBUG
				{
//					 check memory usage on the GPU
					std::lock_guard<mutex> lock(Prismatic::mem_lock);
					size_t free_mem, total_mem;
					free_mem=total_mem=0;
					cudaErrchk(cudaMemGetInfo(&free_mem, &total_mem));
					pars.max_mem = std::max(total_mem - free_mem, pars.max_mem);
				}
#endif // NDEBUG

				size_t Nstart,Nstop;
				Nstart=Nstop=0;
				while (dispatcher.getWork(Nstart, Nstop, pars.meta.batch_size_GPU)){ // synchronously get work assignment
					while (Nstart < Nstop){
						if (Nstart % PRISMATIC_PRINT_FREQUENCY_PROBES < pars.meta.batch_size_GPU | Nstart == 100){
							cout << "Computing Probe Position #" << Nstart << "/" << pars.xp.size() * pars.yp.size() << endl;
						}
						getMultisliceProbe_GPU_streaming_batch(pars, current_trans_ds, cuda_pars.trans_ph, current_PsiProbeInit_d, current_psi_ds,
						                                       current_output_ph, current_psi_intensity_ds,
						                                       current_integratedOutput_ds, current_qya_d, current_qxa_d,
						                                       current_prop_d, Nstart, Nstop, pars.psiProbeInit.get_dimj(), pars.psiProbeInit.get_dimi(),
						                                       current_alphaInd_d, current_cufft_plan, current_stream);
#ifdef PRISMATIC_BUILDING_GUI
						pars.progressbar->signalOutputUpdate(Nstart, pars.xp.size() * pars.yp.size());
#endif
						Nstart = Nstop;
					}
				}
				cout << "GPU worker on stream #" << stream_count << " of GPU #" << GPU_num << "finished\n";
			}));
			++stream_count;
		}

		// now launch CPU work
		if (pars.meta.also_do_CPU_work){
			PRISMATIC_FFTW_INIT_THREADS();
			PRISMATIC_FFTW_PLAN_WITH_NTHREADS(pars.meta.NUM_THREADS);vector<thread> workers_CPU;
			workers_CPU.reserve(pars.meta.NUM_THREADS); // prevents multiple reallocations
			for (auto t = 0; t < pars.meta.NUM_THREADS; ++t) {
				cout << "Launching CPU worker #" << t << endl;
				// push_back is better whenever constructing a new object
				workers_CPU.push_back(thread([&pars, &dispatcher, t, &PRISMATIC_PRINT_FREQUENCY_PROBES]() {
					size_t Nstart, Nstop, early_CPU_stop;
					Nstart=Nstop=0;
					// stop the CPU workers earlier than the GPU ones to prevent slower workers taking the last jobs and having to
					// wait longer for everything to complete
					if (pars.meta.NUM_GPUS > 0){
						// if there are no GPUs, make sure to do all work on CPU
						early_CPU_stop = (size_t)std::max((PRISMATIC_FLOAT_PRECISION)0.0, pars.xp.size() * pars.yp.size() - pars.meta.gpu_cpu_ratio);
					} else {
						early_CPU_stop = pars.xp.size() * pars.yp.size();
					}
					if (dispatcher.getWork(Nstart, Nstop, pars.meta.batch_size_CPU, early_CPU_stop)) { // synchronously get work assignment
						Array1D<std::complex<PRISMATIC_FLOAT_PRECISION> > psi_stack = zeros_ND<1, complex<PRISMATIC_FLOAT_PRECISION> >({{pars.psiProbeInit.size() * pars.meta.batch_size_CPU}});

						// setup batch FFTW parameters
						const int rank = 2;
						int n[] = {(int)pars.psiProbeInit.get_dimj(), (int)pars.psiProbeInit.get_dimi()};
						const int howmany = pars.meta.batch_size_CPU;
						int idist = n[0]*n[1];
						int odist = n[0]*n[1];
						int istride = 1;
						int ostride = 1;
						int *inembed = n;
						int *onembed = n;
						unique_lock<mutex> gatekeeper(fftw_plan_lock);
						PRISMATIC_FFTW_PLAN plan_forward = PRISMATIC_FFTW_PLAN_DFT_BATCH(rank, n, howmany,
						                                                         reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&psi_stack[0]), inembed,
						                                                         istride, idist,
						                                                         reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&psi_stack[0]), onembed,
						                                                         ostride, odist,
						                                                         FFTW_FORWARD, FFTW_MEASURE);
						PRISMATIC_FFTW_PLAN plan_inverse = PRISMATIC_FFTW_PLAN_DFT_BATCH(rank, n, howmany,
						                                                         reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&psi_stack[0]), inembed,
						                                                         istride, idist,
						                                                         reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&psi_stack[0]), onembed,
						                                                         ostride, odist,
						                                                         FFTW_BACKWARD, FFTW_MEASURE);
						gatekeeper.unlock();

						// main work loop
						do {
							while (Nstart < Nstop) {
								if (Nstart % PRISMATIC_PRINT_FREQUENCY_PROBES  < pars.meta.batch_size_CPU | Nstart == 100){
									cout << "Computing Probe Position #" << Nstart << "/" << pars.xp.size() * pars.yp.size() << endl;
								}
//							getMultisliceProbe_CPU(pars, ay, ax, plan_forward, plan_inverse, psi);
								getMultisliceProbe_CPU_batch(pars, Nstart, Nstop, plan_forward, plan_inverse, psi_stack);
#ifdef PRISMATIC_BUILDING_GUI
								pars.progressbar->signalOutputUpdate(Nstart, pars.xp.size() * pars.yp.size());
#endif
								Nstart=Nstop;
							}
							if (Nstop >= early_CPU_stop) break;
						} while(dispatcher.getWork(Nstart, Nstop, pars.meta.batch_size_CPU, early_CPU_stop));
						gatekeeper.lock();
						PRISMATIC_FFTW_DESTROY_PLAN(plan_forward);
						PRISMATIC_FFTW_DESTROY_PLAN(plan_inverse);
						gatekeeper.unlock();
					}
					cout << "CPU worker #" << t << " finished\n";
				}));
			}
			cout << "Waiting on GPU threads..." << endl;
			for (auto& t:workers_CPU)t.join();
			PRISMATIC_FFTW_CLEANUP_THREADS();
		}
		// synchronize threads
		cout << "Waiting on GPU threads..." << endl;
		for (auto& t:workers_GPU)t.join();

		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g){
			cudaSetDevice(g);
			cudaDeviceSynchronize();
		}
	}

	inline void cleanupMemoryM(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
	                           CudaParameters<PRISMATIC_FLOAT_PRECISION> &cuda_pars){
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;
		// synchronize GPUs and cleanup data
		for (auto j = 0; j < pars.meta.NUM_GPUS; ++j){
			cudaErrchk(cudaSetDevice(j));
			cudaErrchk(cudaFree(cuda_pars.PsiProbeInit_d[j]));
			cudaErrchk(cudaFree(cuda_pars.trans_d[j]));
			cudaErrchk(cudaFree(cuda_pars.qxa_d[j]));
			cudaErrchk(cudaFree(cuda_pars.qya_d[j]));
			cudaErrchk(cudaFree(cuda_pars.prop_d[j]));
			cudaErrchk(cudaFree(cuda_pars.alphaInd_d[j]));
		}

		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaSetDevice(s % pars.meta.NUM_GPUS));
			cudaErrchk(cudaFree(cuda_pars.psi_ds[s]));
			cudaErrchk(cudaFree(cuda_pars.psi_intensity_ds[s]));
			cudaErrchk(cudaFree(cuda_pars.integratedOutput_ds[s]));
			cufftErrchk(cufftDestroy(cuda_pars.cufft_plans[s]));
		}

		// free pinned memory
		cudaErrchk(cudaFreeHost(cuda_pars.PsiProbeInit_ph));
		cudaErrchk(cudaFreeHost(cuda_pars.trans_ph));
		cudaErrchk(cudaFreeHost(cuda_pars.prop_ph));
		cudaErrchk(cudaFreeHost(cuda_pars.qxa_ph));
		cudaErrchk(cudaFreeHost(cuda_pars.qya_ph));
		cudaErrchk(cudaFreeHost(cuda_pars.alphaInd_ph));
		for (auto s =0; s < total_num_streams; ++s){
			cudaErrchk(cudaFreeHost(cuda_pars.output_ph[s]));
		}

		// destroy CUDA streams
		for (auto j = 0; j < total_num_streams; ++j){
			cudaSetDevice(j % pars.meta.NUM_GPUS);
			cudaErrchk(cudaStreamDestroy(cuda_pars.streams[j]));
		}
		for (auto j = 0; j < pars.meta.NUM_GPUS; ++j) {
			cudaErrchk(cudaSetDevice(j));
			cudaErrchk(cudaDeviceReset());
		}
		delete[] cuda_pars.streams;
		delete[] cuda_pars.cufft_plans;
		delete[] cuda_pars.PsiProbeInit_d;
		delete[] cuda_pars.trans_d;
		delete[] cuda_pars.prop_d;
		delete[] cuda_pars.qxa_d;
		delete[] cuda_pars.qya_d;
		delete[] cuda_pars.alphaInd_d;
		delete[] cuda_pars.psi_ds;
		delete[] cuda_pars.psi_intensity_ds;
		delete[] cuda_pars.integratedOutput_ds;
		delete[] cuda_pars.output_ph;
	}

	// computes the result of probe position ay,ax using the GPU. The effect of this function is the same as getMultisliceProbe_CPU
	__host__ void getMultisliceProbe_GPU_singlexfer(Parameters<PRISMATIC_FLOAT_PRECISION>& pars,
	                                                PRISMATIC_CUDA_COMPLEX_FLOAT* trans_d,
	                                                PRISMATIC_CUDA_COMPLEX_FLOAT* PsiProbeInit_d,
	                                                PRISMATIC_CUDA_COMPLEX_FLOAT* psi_ds,
	                                                PRISMATIC_FLOAT_PRECISION* output_ph,
	                                                PRISMATIC_FLOAT_PRECISION* psi_intensity_ds,
	                                                PRISMATIC_FLOAT_PRECISION* integratedOutput_ds,
									                const PRISMATIC_FLOAT_PRECISION* qya_d,
									                const PRISMATIC_FLOAT_PRECISION* qxa_d,
									                const PRISMATIC_CUDA_COMPLEX_FLOAT* prop_d,
									                const size_t ay,
									                const size_t ax,
									                const size_t dimj,
									                const size_t dimi,
									                const PRISMATIC_FLOAT_PRECISION* alphaInd_d,
									                const cufftHandle& plan,
									                cudaStream_t& stream){

		// initialize psi
		PRISMATIC_FLOAT_PRECISION yp = pars.yp[ay];
		PRISMATIC_FLOAT_PRECISION xp = pars.xp[ax];
		const size_t psi_size = dimj*dimi;
		initializePsi<<<(psi_size - 1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_ds, PsiProbeInit_d, qya_d, qxa_d, psi_size, yp, xp);
		for (auto planeNum = 0; planeNum < pars.numPlanes; ++planeNum) {
			cufftErrchk(PRISMATIC_CUFFT_EXECUTE(plan, &psi_ds[0], &psi_ds[0], CUFFT_INVERSE));
			multiply_inplace<<<(psi_size - 1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_ds, &trans_d[planeNum*psi_size], psi_size);
			cufftErrchk(PRISMATIC_CUFFT_EXECUTE(plan, &psi_ds[0], &psi_ds[0], CUFFT_FORWARD));
			multiply_inplace<<<(psi_size - 1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_ds, prop_d, psi_size);
			divide_inplace<<<(psi_size - 1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_ds, PRISMATIC_MAKE_CU_COMPLEX(psi_size, 0), psi_size);
		}
		abs_squared<<<(psi_size - 1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_intensity_ds, psi_ds, psi_size);
		formatOutput_GPU_integrate(pars, psi_intensity_ds, alphaInd_d, output_ph, integratedOutput_ds, ay, ax, dimj, dimi, stream);
}

	__host__ void getMultisliceProbe_GPU_singlexfer_batch(Parameters<PRISMATIC_FLOAT_PRECISION>& pars,
	                                                      PRISMATIC_CUDA_COMPLEX_FLOAT* trans_d,
	                                                      PRISMATIC_CUDA_COMPLEX_FLOAT* PsiProbeInit_d,
	                                                      PRISMATIC_CUDA_COMPLEX_FLOAT* psi_ds,
	                                                      PRISMATIC_FLOAT_PRECISION* output_ph,
	                                                      PRISMATIC_FLOAT_PRECISION* psi_intensity_ds,
	                                                      PRISMATIC_FLOAT_PRECISION* integratedOutput_ds,
	                                                      const PRISMATIC_FLOAT_PRECISION* qya_d,
	                                                      const PRISMATIC_FLOAT_PRECISION* qxa_d,
	                                                      const PRISMATIC_CUDA_COMPLEX_FLOAT* prop_d,
	                                                      const size_t Nstart,
	                                                      const size_t Nstop,
	                                                      const size_t dimj,
	                                                      const size_t dimi,
	                                                      const PRISMATIC_FLOAT_PRECISION* alphaInd_d,
	                                                      const cufftHandle& plan,
	                                                      cudaStream_t& stream){
		const size_t psi_size = dimj*dimi;
		for (auto batch_idx = 0; batch_idx < (Nstop-Nstart); ++batch_idx) {
			const size_t ay = (Nstart + batch_idx) / pars.xp.size();
			const size_t ax = (Nstart + batch_idx) % pars.xp.size();

			// initialize psi
			PRISMATIC_FLOAT_PRECISION yp = pars.yp[ay];
			PRISMATIC_FLOAT_PRECISION xp = pars.xp[ax];

			initializePsi << < (psi_size - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >> >
                (psi_ds + (batch_idx * psi_size), PsiProbeInit_d, qya_d, qxa_d, psi_size, yp, xp);
		}
		for (auto planeNum = 0; planeNum < pars.numPlanes; ++planeNum) {
			cufftErrchk(PRISMATIC_CUFFT_EXECUTE(plan, &psi_ds[0], &psi_ds[0], CUFFT_INVERSE));
			for (auto batch_idx = 0; batch_idx < (Nstop-Nstart); ++batch_idx) {
				multiply_inplace << < (psi_size - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >> >
						(psi_ds + (batch_idx * psi_size), &trans_d[planeNum * psi_size], psi_size);
			}
			cufftErrchk(PRISMATIC_CUFFT_EXECUTE(plan, &psi_ds[0], &psi_ds[0], CUFFT_FORWARD));
			for (auto batch_idx = 0; batch_idx < (Nstop-Nstart); ++batch_idx) {
				multiply_inplace << < (psi_size - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >> >
						(psi_ds + (batch_idx * psi_size), prop_d, psi_size);
				divide_inplace << < (psi_size - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >> >
						(psi_ds + (batch_idx * psi_size), PRISMATIC_MAKE_CU_COMPLEX(psi_size, 0), psi_size);
			}
		}
		abs_squared << < ( psi_size*(Nstop-Nstart) - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >> > (psi_intensity_ds, psi_ds, psi_size*(Nstop-Nstart));
		for (auto batch_idx = 0; batch_idx < (Nstop-Nstart); ++batch_idx) {
			const size_t ay = (Nstart + batch_idx) / pars.xp.size();
			const size_t ax = (Nstart + batch_idx) % pars.xp.size();
			formatOutput_GPU_integrate(pars, psi_intensity_ds + (batch_idx * psi_size),
			                           alphaInd_d, output_ph, integratedOutput_ds, ay, ax, dimj, dimi, stream);
		}
	}

	__host__ void getMultisliceProbe_GPU_streaming(Parameters<PRISMATIC_FLOAT_PRECISION>& pars,
	                                               PRISMATIC_CUDA_COMPLEX_FLOAT* trans_d,
	                                               const complex<PRISMATIC_FLOAT_PRECISION>* trans_ph,
	                                               PRISMATIC_CUDA_COMPLEX_FLOAT* PsiProbeInit_d,
	                                               PRISMATIC_CUDA_COMPLEX_FLOAT* psi_ds,
	                                               PRISMATIC_FLOAT_PRECISION* output_ph,
	                                               PRISMATIC_FLOAT_PRECISION* psi_intensity_ds,
	                                               PRISMATIC_FLOAT_PRECISION* integratedOutput_ds,
	                                               const PRISMATIC_FLOAT_PRECISION* qya_d,
	                                               const PRISMATIC_FLOAT_PRECISION* qxa_d,
	                                               const PRISMATIC_CUDA_COMPLEX_FLOAT* prop_d,
	                                               const size_t& ay,
	                                               const size_t& ax,
	                                               const size_t dimj,
	                                               const size_t dimi,
	                                               const PRISMATIC_FLOAT_PRECISION* alphaInd_d,
	                                               const cufftHandle& plan,
	                                               cudaStream_t& stream){
		// initialize psi
		PRISMATIC_FLOAT_PRECISION yp = pars.yp[ay];
		PRISMATIC_FLOAT_PRECISION xp = pars.xp[ax];
		const size_t psi_size = dimj*dimi;
		initializePsi<<<(psi_size - 1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_ds, PsiProbeInit_d, qya_d, qxa_d, dimj*dimi, yp, xp);


		for (auto planeNum = 0; planeNum < pars.numPlanes; ++planeNum) {
			cudaErrchk(cudaMemcpyAsync(trans_d, &trans_ph[planeNum*psi_size], psi_size * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT), cudaMemcpyHostToDevice, stream));
			cufftErrchk(PRISMATIC_CUFFT_EXECUTE(plan, &psi_ds[0], &psi_ds[0], CUFFT_INVERSE));
			multiply_inplace<<<(psi_size - 1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_ds, trans_d, psi_size);
			cufftErrchk(PRISMATIC_CUFFT_EXECUTE(plan, &psi_ds[0], &psi_ds[0], CUFFT_FORWARD));
			multiply_inplace<<<(psi_size - 1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_ds, prop_d, psi_size);
			divide_inplace<<<(psi_size - 1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_ds, PRISMATIC_MAKE_CU_COMPLEX(psi_size, 0), psi_size);
		}
		abs_squared<<<(psi_size - 1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_intensity_ds, psi_ds, psi_size);
		formatOutput_GPU_integrate(pars, psi_intensity_ds, alphaInd_d, output_ph, integratedOutput_ds, ay, ax, dimj, dimi,stream);
	}


	__host__ void getMultisliceProbe_GPU_streaming_batch(Parameters<PRISMATIC_FLOAT_PRECISION>& pars,
	                                                     PRISMATIC_CUDA_COMPLEX_FLOAT* trans_d,
	                                                     const complex<PRISMATIC_FLOAT_PRECISION>* trans_ph,
	                                                     PRISMATIC_CUDA_COMPLEX_FLOAT* PsiProbeInit_d,
	                                                     PRISMATIC_CUDA_COMPLEX_FLOAT* psi_ds,
	                                                     PRISMATIC_FLOAT_PRECISION* output_ph,
	                                                     PRISMATIC_FLOAT_PRECISION* psi_intensity_ds,
	                                                     PRISMATIC_FLOAT_PRECISION* integratedOutput_ds,
	                                                     const PRISMATIC_FLOAT_PRECISION* qya_d,
	                                                     const PRISMATIC_FLOAT_PRECISION* qxa_d,
	                                                     const PRISMATIC_CUDA_COMPLEX_FLOAT* prop_d,
	                                                     const size_t Nstart,
	                                                     const size_t Nstop,
	                                                     const size_t dimj,
	                                                     const size_t dimi,
	                                                     const PRISMATIC_FLOAT_PRECISION* alphaInd_d,
	                                                     const cufftHandle& plan,
	                                                     cudaStream_t& stream){

		// initialize psi
		const size_t psi_size = dimj*dimi;
		for (auto batch_idx = 0; batch_idx < (Nstop-Nstart); ++batch_idx) {
			const size_t ay = (Nstart + batch_idx) / pars.xp.size();
			const size_t ax = (Nstart + batch_idx) % pars.xp.size();
			PRISMATIC_FLOAT_PRECISION yp = pars.yp[ay];
			PRISMATIC_FLOAT_PRECISION xp = pars.xp[ax];
			initializePsi << < (psi_size - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >> >
			                                                (psi_ds + (batch_idx * psi_size), PsiProbeInit_d, qya_d, qxa_d, psi_size, yp, xp);
		}

		for (auto planeNum = 0; planeNum < pars.numPlanes; ++planeNum) {

			cudaErrchk(cudaMemcpyAsync(trans_d, &trans_ph[planeNum*psi_size], psi_size * sizeof(PRISMATIC_CUDA_COMPLEX_FLOAT), cudaMemcpyHostToDevice, stream));
			cufftErrchk(PRISMATIC_CUFFT_EXECUTE(plan, &psi_ds[0], &psi_ds[0], CUFFT_INVERSE));
			for (auto batch_idx = 0; batch_idx < (Nstop-Nstart); ++batch_idx) {
				multiply_inplace << < (psi_size - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >> >
				                                                   (psi_ds + (batch_idx * psi_size), trans_d, psi_size);
			}
			cufftErrchk(PRISMATIC_CUFFT_EXECUTE(plan, &psi_ds[0], &psi_ds[0], CUFFT_FORWARD));
			for (auto batch_idx = 0; batch_idx < (Nstop-Nstart); ++batch_idx) {
				multiply_inplace << < (psi_size - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >> >
				                                                   (psi_ds + (batch_idx * psi_size), prop_d, psi_size);
				divide_inplace << < (psi_size - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >> >
				                                                 (psi_ds + (batch_idx * psi_size), PRISMATIC_MAKE_CU_COMPLEX(psi_size, 0), psi_size);
			}
		}
		abs_squared << < (psi_size*(Nstop-Nstart) - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >> > (psi_intensity_ds, psi_ds, psi_size*(Nstop-Nstart));
		for (auto batch_idx = 0; batch_idx < (Nstop-Nstart); ++batch_idx) {
			const size_t ay = (Nstart + batch_idx) / pars.xp.size();
			const size_t ax = (Nstart + batch_idx) % pars.xp.size();
			formatOutput_GPU_integrate(pars, psi_intensity_ds + (batch_idx * psi_size),
			                           alphaInd_d, output_ph, integratedOutput_ds, ay, ax, dimj, dimi, stream);
		}
	}
    __host__ void buildMultisliceOutput_GPU_singlexfer(Parameters <PRISMATIC_FLOAT_PRECISION> &pars){
		using namespace std;
#ifdef PRISMATIC_BUILDING_GUI
	    pars.progressbar->signalDescriptionMessage("Computing final output (Multislice)");
#endif
		CudaParameters<PRISMATIC_FLOAT_PRECISION> cuda_pars;

		// determine the batch size to use
	    pars.meta.batch_size_GPU = min(pars.meta.batch_size_target_GPU, max((size_t)1, pars.xp.size() * pars.yp.size() / max((size_t)1, (pars.meta.NUM_STREAMS_PER_GPU*pars.meta.NUM_GPUS))));

		// populate the Multislice output stack dividing the work between GPUs and CPU cores.
		// this version assumes the full trans array fits into DRAM on each GPU

		createPlansAndStreamsM(pars, cuda_pars);

	    // create page-locked (pinned) host memory buffers
	    allocatePinnedHostMemory_M(pars, cuda_pars);

	    copyToPinnedMemory_M(pars, cuda_pars);

	    // allocate memory on the GPUs
	    allocateDeviceMemory_singlexferM(pars, cuda_pars);

	    // copy memory to GPUs
	    copyToGPUMemory_singlexferM(pars, cuda_pars);

	    // launch GPU and CPU workers
	    launchWorkers_singlexferM(pars, cuda_pars);

	    // free memory
	    cleanupMemoryM(pars, cuda_pars);
	}

	__host__ void buildMultisliceOutput_GPU_streaming(Parameters <PRISMATIC_FLOAT_PRECISION> &pars){
#ifdef PRISMATIC_BUILDING_GUI
		pars.progressbar->signalDescriptionMessage("Computing final output (Multislice)");
#endif
		// populate the Multislice output stack dividing the work between GPUs and CPU cores.
		// this version assumes the full trans array fits into DRAM on each GPU

		using namespace std;

		CudaParameters<PRISMATIC_FLOAT_PRECISION> cuda_pars;

		// determine the batch size to use
		pars.meta.batch_size_GPU = min(pars.meta.batch_size_target_GPU, max((size_t)1, pars.xp.size() * pars.yp.size() / max((size_t)1, (pars.meta.NUM_STREAMS_PER_GPU*pars.meta.NUM_GPUS))));

		// populate the Multislice output stack dividing the work between GPUs and CPU cores.
		// this version assumes the full trans array fits into DRAM on each GPU

		createPlansAndStreamsM(pars, cuda_pars);

		allocatePinnedHostMemory_M(pars, cuda_pars);

		copyToPinnedMemory_M(pars,cuda_pars);

		// allocate memory on the GPUs
		allocateDeviceMemory_streamingM(pars, cuda_pars);

		// copy memory to GPUs
		copyToGPUMemory_streamingM(pars, cuda_pars);

		// launch GPU and CPU workers
		launchWorkers_streamingM(pars, cuda_pars);

		// free memory
		cleanupMemoryM(pars, cuda_pars);
	}
}