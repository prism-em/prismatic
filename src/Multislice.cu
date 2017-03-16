// Calculate result of Multislice simulation using GPU and (potentially) CPU. Multiple GPU threads are launched, each with
// their own memory buffers. Page-locked host memory is allocated so that memory transfers to the GPU can occur asynchronously,
// and memory allocation for the GPU occurs only once, as each call to cudaMalloc will potentially interrupt concurrent execution.
// Each GPU/CPU worker thread repeatedly calls getWorkID to be assigned probe positions to compute. This queue mechanism
// ensures that both the CPU and GPU are kept busy.

// For variable naming, the suffixes are "_d" for "device" (1 copy per GPU), "_ds" for "device stream (1 copy per stream), "_ph" for "pinned host"

#include "Multislice.cuh"
#include "Multislice.h"
#include "cuComplex.h"
#include "cufft.h"
#include "getWorkID.h"
#include <iostream>
#include "fftw3.h"
#include "utility.h"
#include "utility.cuh"

namespace PRISM{

	// computes the result of probe position ay,ax using the GPU. The effect of this function is the same as getMultisliceProbe_CPU
	__host__ void getMultisliceProbe_GPU(Parameters<PRISM_FLOAT_PRECISION>& pars,
									    PRISM_CUDA_COMPLEX_FLOAT* trans_d,
									    PRISM_CUDA_COMPLEX_FLOAT* PsiProbeInit_d,
									    PRISM_CUDA_COMPLEX_FLOAT* psi_ds,
									    PRISM_FLOAT_PRECISION* output_ph,
									    PRISM_FLOAT_PRECISION* psi_intensity_ds,
									    PRISM_FLOAT_PRECISION* integratedOutput_ds,
									    const PRISM_FLOAT_PRECISION* qya_d,
									    const PRISM_FLOAT_PRECISION* qxa_d,
									    const PRISM_CUDA_COMPLEX_FLOAT* prop_d,
									    const size_t& ay,
									    const size_t& ax,
									    const size_t dimj,
									    const size_t dimi,
									    const PRISM_FLOAT_PRECISION* alphaInd_d,
									    const cufftHandle& plan,
									    cudaStream_t& stream){

		// initialize psi
		PRISM_FLOAT_PRECISION yp = pars.yp[ay];
		PRISM_FLOAT_PRECISION xp = pars.xp[ax];
		const size_t N = dimj*dimi;
		initializePsi<<<(N-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_ds, PsiProbeInit_d, qya_d, qxa_d, dimj*dimi, yp, xp);

		for (auto planeNum = 0; planeNum < pars.numPlanes; ++planeNum) {
			cufftErrchk(PRISM_CUFFT_EXECUTE(plan, &psi_ds[0], &psi_ds[0], CUFFT_INVERSE));
			multiply_inplace<<<(N-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_ds, &trans_d[planeNum*N], N);
			cufftErrchk(PRISM_CUFFT_EXECUTE(plan, &psi_ds[0], &psi_ds[0], CUFFT_FORWARD));
			multiply_inplace<<<(N-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_ds, prop_d, N);
			divide_inplace<<<(N-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_ds, PRISM_MAKE_CU_COMPLEX(N, 0), N);
		}
		abs_squared<<<(N-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_intensity_ds, psi_ds, N);
		formatOutput_GPU_integrate(pars, psi_intensity_ds, alphaInd_d, output_ph, integratedOutput_ds, ay, ax, dimj, dimi,stream);
}
    __host__ void buildMultisliceOutput_GPU(Parameters <PRISM_FLOAT_PRECISION> &pars,
                                            Array3D <std::complex<PRISM_FLOAT_PRECISION>> &trans,
                                            Array2D <std::complex<PRISM_FLOAT_PRECISION>> &PsiProbeInit,
                                            Array2D <PRISM_FLOAT_PRECISION> &alphaInd) {


		cudaErrchk(cudaDeviceReset());
		cudaErrchk(cudaFree(0));
		cout << "debug pars.prop" << endl;
		cout << "Psi dim y = " << PsiProbeInit.get_dimj() << endl;
		cout << "Psi dim x = " << PsiProbeInit.get_dimi() << endl;
		cout << "number of planes = " << pars.numPlanes << endl;
		for (auto i =0; i < 10; ++i)cout << pars.prop[i] << endl;
		// populate the Multislice output stack dividing the work between GPUs and CPU cores.
		// this version assumes the full trans array fits into DRAM on each GPU

		using namespace std;

		// create CUDA streams
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;
		cudaStream_t streams[total_num_streams];
		cufftHandle cufft_plan[total_num_streams];
		cout <<"total_num_streams = " << total_num_streams<< endl;
		for (auto j = 0; j < total_num_streams; ++j){
			cudaSetDevice(j % pars.meta.NUM_GPUS);
			cudaErrchk(cudaStreamCreate(&streams[j]));
			cufftErrchk(cufftPlan2d(&cufft_plan[j], PsiProbeInit.get_dimi(), PsiProbeInit.get_dimj(), CUFFT_C2C));
			cufftErrchk(cufftSetStream(cufft_plan[j], streams[j]));
		}


		vector<thread> workers_GPU;
		workers_GPU.reserve(total_num_streams); // prevents multiple reallocations


		// pointers to pinned host memory for async transfers
		std::complex<PRISM_FLOAT_PRECISION>  *PsiProbeInit_ph;
		std::complex<PRISM_FLOAT_PRECISION>  *trans_ph;
		std::complex<PRISM_FLOAT_PRECISION>  *prop_ph;
		PRISM_FLOAT_PRECISION                *qxa_ph;
		PRISM_FLOAT_PRECISION                *qya_ph;
		PRISM_FLOAT_PRECISION                *alphaInd_ph;
		PRISM_FLOAT_PRECISION                *output_ph[total_num_streams];

		// allocate pinned memory
		cudaErrchk(cudaMallocHost((void **)&PsiProbeInit_ph, PsiProbeInit.size()*sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
		cudaErrchk(cudaMallocHost((void **)&trans_ph,        trans.size()*sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
		cudaErrchk(cudaMallocHost((void **)&prop_ph,         pars.prop.size()*sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
		cudaErrchk(cudaMallocHost((void **)&qxa_ph,          pars.qxa.size()*sizeof(PRISM_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **)&qya_ph,          pars.qya.size()*sizeof(PRISM_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **)&alphaInd_ph,     alphaInd.size()*sizeof(PRISM_FLOAT_PRECISION)));
		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaMallocHost((void **) &output_ph[s], pars.stack.get_dimj() * pars.stack.get_dimi() * sizeof(PRISM_FLOAT_PRECISION)));
		}
		// copy host memory to pinned
		memcpy(PsiProbeInit_ph, &PsiProbeInit[0], PsiProbeInit.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>));
		memcpy(trans_ph,        &trans[0],        trans.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>));
		memcpy(prop_ph,         &pars.prop[0],    pars.prop.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>));
		memcpy(qxa_ph,          &pars.qxa[0],     pars.qxa.size() * sizeof(PRISM_FLOAT_PRECISION));
		memcpy(qya_ph,          &pars.qya[0],     pars.qya.size() * sizeof(PRISM_FLOAT_PRECISION));
		memcpy(alphaInd_ph,     &alphaInd[0],     alphaInd.size() * sizeof(PRISM_FLOAT_PRECISION));


		// pointers to read-only GPU memory (one copy per GPU)
		PRISM_CUDA_COMPLEX_FLOAT *PsiProbeInit_d[pars.meta.NUM_GPUS];
		PRISM_CUDA_COMPLEX_FLOAT *trans_d[pars.meta.NUM_GPUS];
		PRISM_CUDA_COMPLEX_FLOAT *prop_d[pars.meta.NUM_GPUS];
		PRISM_FLOAT_PRECISION    *qxa_d[pars.meta.NUM_GPUS];
		PRISM_FLOAT_PRECISION    *qya_d[pars.meta.NUM_GPUS];
	    PRISM_FLOAT_PRECISION    *alphaInd_d[pars.meta.NUM_GPUS];

		// pointers to read/write GPU memory (one per stream)
		PRISM_CUDA_COMPLEX_FLOAT *psi_ds[total_num_streams];
		PRISM_FLOAT_PRECISION    *psi_intensity_ds[total_num_streams];
		PRISM_FLOAT_PRECISION    *integratedOutput_ds[total_num_streams];

		// allocate memory on each GPU
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaMalloc((void **) &PsiProbeInit_d[g],     PsiProbeInit.size()        * sizeof(PsiProbeInit[0])));
			cudaErrchk(cudaMalloc((void **) &trans_d[g],            trans.size()               * sizeof(trans[0])));
			cudaErrchk(cudaMalloc((void **) &prop_d[g],             pars.prop.size()           * sizeof(pars.prop[0])));
			cudaErrchk(cudaMalloc((void **) &qxa_d[g],              pars.qxa.size()            * sizeof(pars.qxa[0])));
			cudaErrchk(cudaMalloc((void **) &qya_d[g],              pars.qya.size()            * sizeof(pars.qya[0])));
			cudaErrchk(cudaMalloc((void **) &alphaInd_d[g],         alphaInd.size()            * sizeof(alphaInd[0])));
		}

		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaSetDevice(s % pars.meta.NUM_GPUS));
			cudaErrchk(cudaMalloc((void **) &psi_ds[s],              PsiProbeInit.size()        * sizeof(PsiProbeInit[0])));
			cudaErrchk(cudaMalloc((void **) &psi_intensity_ds[s],    PsiProbeInit.size()        * sizeof(PRISM_FLOAT_PRECISION)));
			cudaErrchk(cudaMalloc((void **) &integratedOutput_ds[s], pars.detectorAngles.size() * sizeof(PRISM_FLOAT_PRECISION)));
			cudaErrchk(cudaMemset(psi_ds[s], 0, PsiProbeInit.size()        * sizeof(PsiProbeInit[0])));
			cudaErrchk(cudaMemset(psi_intensity_ds[s], 0, PsiProbeInit.size()        * sizeof(PRISM_FLOAT_PRECISION)));
			cudaErrchk(cudaMemset(integratedOutput_ds[s], 0, pars.detectorAngles.size() * sizeof(PRISM_FLOAT_PRECISION)));
		}


		// Copy memory to each GPU asynchronously from the pinned host memory spaces.
		// The streams are laid out so that consecutive streams represent different GPUs. If we
		// have more than one stream per GPU, then we want to interleave as much as possible
		int stream_id = 0;
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			stream_id = g;
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaMemcpyAsync(PsiProbeInit_d[g], &PsiProbeInit_ph[0],
			                      PsiProbeInit.size() * sizeof(PsiProbeInit[0]), cudaMemcpyHostToDevice, streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cout << "stream_id = " << stream_id << endl;
			cudaErrchk(cudaMemcpyAsync(trans_d[g], &trans_ph[0],
			                      trans.size() * sizeof(trans[0]), cudaMemcpyHostToDevice, streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cout << "stream_id = " << stream_id << endl;
			cudaErrchk(cudaMemcpyAsync(prop_d[g], &prop_ph[0],
			                      pars.prop.size() * sizeof(pars.prop[0]), cudaMemcpyHostToDevice, streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cout << "stream_id = " << stream_id << endl;
			cudaErrchk(cudaMemcpyAsync(qxa_d[g], &qxa_ph[0],
			                      pars.qxa.size() * sizeof(pars.qxa[0]), cudaMemcpyHostToDevice, streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cout << "stream_id = " << stream_id << endl;
			cudaErrchk(cudaMemcpyAsync(qya_d[g], &qya_ph[0],
			                      pars.qya.size() * sizeof(pars.qya[0]), cudaMemcpyHostToDevice, streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cout << "stream_id = " << stream_id << endl;
			cudaErrchk(cudaMemcpyAsync(alphaInd_d[g], &alphaInd_ph[0],
			                      alphaInd.size() * sizeof(alphaInd[0]), cudaMemcpyHostToDevice, streams[stream_id]));
		}

		// make sure transfers are complete
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g){
			cudaSetDevice(g);
			cudaDeviceSynchronize();
		}

		size_t psi_size = PsiProbeInit.size();
		int stream_count = 0;
		setWorkStartStop(0, pars.xp.size() * pars.yp.size());
		for (auto t = 0; t < total_num_streams; ++t){
			int GPU_num = stream_count % pars.meta.NUM_GPUS; // determine which GPU handles this job
			cudaStream_t& current_stream = streams[stream_count];
			cout << "Launching GPU worker on stream #" << stream_count << " of GPU #" << GPU_num << '\n';

			// get pointers to the pre-copied arrays, making sure to get those on the current GPU
			PRISM_CUDA_COMPLEX_FLOAT *current_PsiProbeInit_d = PsiProbeInit_d[GPU_num];
			PRISM_CUDA_COMPLEX_FLOAT *current_trans_d = trans_d[GPU_num];
			PRISM_CUDA_COMPLEX_FLOAT *current_prop_d  = prop_d[GPU_num];
			PRISM_FLOAT_PRECISION *current_qxa_d      = qxa_d[GPU_num];
			PRISM_FLOAT_PRECISION *current_qya_d      = qya_d[GPU_num];
			PRISM_FLOAT_PRECISION *current_alphaInd_d = alphaInd_d[GPU_num];

			// get pointers to per-stream arrays
			PRISM_CUDA_COMPLEX_FLOAT *current_psi_ds           = psi_ds[stream_count];
			PRISM_FLOAT_PRECISION *current_psi_intensity_ds    = psi_intensity_ds[stream_count];
			PRISM_FLOAT_PRECISION *current_integratedOutput_ds = integratedOutput_ds[stream_count];
			PRISM_FLOAT_PRECISION *current_output_ph           = output_ph[stream_count];
			cufftHandle & current_cufft_plan = cufft_plan[stream_count];
			// launch a new thread
			// emplace_back is better whenever constructing a new object
			workers_GPU.emplace_back(thread([&pars, current_trans_d, current_PsiProbeInit_d, &alphaInd, current_alphaInd_d,
					                                current_psi_ds, current_psi_intensity_ds, current_integratedOutput_ds,
					                                GPU_num, current_qya_d, current_qxa_d, current_output_ph, current_cufft_plan,
					                                current_prop_d, &current_stream, &psi_size, &PsiProbeInit, stream_count]() {

				// set the GPU context
				cudaErrchk(cudaSetDevice(GPU_num)); // set current GPU
				size_t Nstart, Nstop, ay, ax;
				while (getWorkID(pars, Nstart, Nstop)){ // synchronously get work assignment
					while (Nstart != Nstop){
						ay = Nstart / pars.xp.size();
						ax = Nstart % pars.xp.size();
						if (ax==0 & ay==0) {
							cout << "ON GPU" << endl;
							cout << "pars.stack.at(0,0,0) =  " << pars.stack.at(0, 0, 0) << endl;

							complex<float> deb2;
							for (auto jj = 0; jj < 100; ++jj) {
								cudaMemcpy(&deb2, current_psi_ds + jj, sizeof(PRISM_FLOAT_PRECISION),
								           cudaMemcpyDeviceToHost);
								cout << "psi_ds before = " << deb2 << endl;
							}
						}
							getMultisliceProbe_GPU(pars, current_trans_d, current_PsiProbeInit_d, current_psi_ds, current_output_ph,
						                       current_psi_intensity_ds,
						                       current_integratedOutput_ds, current_qya_d, current_qxa_d,
						                       current_prop_d, ay, ax, PsiProbeInit.get_dimj(), PsiProbeInit.get_dimi(),
						                       current_alphaInd_d, current_cufft_plan, current_stream);
						if (ax==0 & ay==0){
							cout <<"ON GPU" << endl;
							cout << "pars.stack.at(0,0,0,0) =  " << pars.stack.at(0,0,0,0) << endl;
							float s = 0;
							for (auto jj = 0; jj < pars.stack.get_dimj(); ++jj)s+=pars.stack.at(0,0,jj,0);
							cout << "sum = " << s << endl;
							if (std::isnan(s) | s < 0.01)
							{
								cout << "MESSED UP!!!\n";
								cout << "MESSED UP!!!\n";
								cout << "MESSED UP!!!\n";
								cout << "MESSED UP!!!\n";
								float ans;
								for (auto jj = 0; jj < pars.detectorAngles.size(); ++jj){
									cudaMemcpy(&ans, current_integratedOutput_ds + jj, sizeof(PRISM_FLOAT_PRECISION), cudaMemcpyDeviceToHost);
									cout << "current_integratedOutput_ds = " << ans << endl;
								}
								for (auto jj = 0; jj < 100; ++jj){
									cudaMemcpy(&ans, current_alphaInd_d + jj, sizeof(PRISM_FLOAT_PRECISION), cudaMemcpyDeviceToHost);
									cout << "current_alphaInd = " << ans << endl;
								}
								for (auto jj = 0; jj < 100; ++jj){
									cudaMemcpy(&ans, current_qya_d + jj, sizeof(PRISM_FLOAT_PRECISION), cudaMemcpyDeviceToHost);
									cout << "qya_d = " << ans << endl;
								}
								for (auto jj = 0; jj < 100; ++jj){
									cudaMemcpy(&ans, current_qxa_d + jj, sizeof(PRISM_FLOAT_PRECISION), cudaMemcpyDeviceToHost);
									cout << "qxa_d = " << ans << endl;
								}
								complex<float> deb;
								for (auto jj = 0; jj < 100; ++jj){
									cudaMemcpy(&deb, current_psi_ds + jj, sizeof(PRISM_FLOAT_PRECISION), cudaMemcpyDeviceToHost);
									cout << "psi_ds = " << deb << endl;
								}
								for (auto jj = 0; jj < 100; ++jj){
									cudaMemcpy(&deb, current_PsiProbeInit_d + jj, sizeof(PRISM_FLOAT_PRECISION), cudaMemcpyDeviceToHost);
									cout << "PsiProbeInit_d = " << deb << endl;
								}



							} //else{
//								float ans;
//								for (auto jj = 0; jj < pars.detectorAngles.size(); ++jj){
//									cudaMemcpy(&ans, current_integratedOutput_ds + jj, sizeof(PRISM_FLOAT_PRECISION), cudaMemcpyDeviceToHost);
//									cout << "current_integratedOutput_ds = " << ans << endl;
//								}
//							}
						}
						++Nstart;
					}
				}
				cout << "GPU worker on stream #" << stream_count << " of GPU #" << GPU_num << "finished\n";
			}));

			++stream_count;
		}


		// now launch CPU work
		PRISM_FFTW_INIT_THREADS();
		PRISM_FFTW_PLAN_WITH_NTHREADS(pars.meta.NUM_THREADS);
		if (pars.meta.also_do_CPU_work){
			vector<thread> workers_CPU;
			workers_CPU.reserve(pars.meta.NUM_THREADS); // prevents multiple reallocations
			for (auto t = 0; t < pars.meta.NUM_THREADS; ++t) {
				cout << "Launching CPU worker #" << t << '\n';
				// emplace_back is better whenever constructing a new object
				workers_CPU.emplace_back(thread([&pars, t, &trans, &alphaInd, &PsiProbeInit]() {
					size_t Nstart, Nstop, ay, ax;
					while (getWorkID(pars, Nstart, Nstop)) { // synchronously get work assignment
						while (Nstart != Nstop) {
							ay = Nstart / pars.xp.size();
							ax = Nstart % pars.xp.size();
							getMultisliceProbe_CPU(pars, trans, PsiProbeInit, ay, ax, alphaInd);
							++Nstart;
						}
					}

					cout << "CPU worker #" << t << " finished\n";
				}));
			}
			cout << "Waiting on GPU threads..." << endl;
			for (auto& t:workers_CPU)t.join();
		}
		// synchronize threads
		cout << "Waiting on GPU threads..." << endl;
		for (auto& t:workers_GPU)t.join();


		PRISM_FFTW_CLEANUP_THREADS();
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g){
			cudaSetDevice(g);
			cudaDeviceSynchronize();
		}

		// copy the results of the GPU, which are in pinned memory, back to the actual stack. The CPU work populates the
		// beginning, so make sure to copy from the offset of where the GPU started. Launch this copy on a background thread
		// while we cleanup the GPU
//		const size_t GPU_start_offset = (size_t)CPU_stop*pars.stack.get_dimk()*pars.stack.get_dimj()*pars.stack.get_dimi();
//		std::thread copy_t([&GPU_start_offset, &pars, &stack_ph](){
//			memcpy(&pars.stack[GPU_start_offset],
//			       &stack_ph[GPU_start_offset],
//			       (pars.stack.size()-GPU_start_offset) * sizeof(PRISM_FLOAT_PRECISION));
//		});

		// synchronize GPUs and cleanup data
		for (auto j = 0; j < pars.meta.NUM_GPUS; ++j){
			cudaErrchk(cudaSetDevice(j));
//			cudaErrchk(cudaDeviceSynchronize());
			cudaErrchk(cudaFree(PsiProbeInit_d[j]));
			cudaErrchk(cudaFree(trans_d[j]));
			cudaErrchk(cudaFree(qxa_d[j]));
			cudaErrchk(cudaFree(qya_d[j]));
			cudaErrchk(cudaFree(prop_d[j]));
			cudaErrchk(cudaFree(alphaInd_d[j]));
//			cudaErrchk(cudaFree(integratedOutput_d[j]));
		}

		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaSetDevice(s % pars.meta.NUM_GPUS));
			cudaErrchk(cudaFree(psi_ds[s]));
			cudaErrchk(cudaFree(psi_intensity_ds[s]));
			cudaErrchk(cudaFree(integratedOutput_ds[s]));
			cufftErrchk(cufftDestroy(cufft_plan[s]));
		}


		// free pinned memory
		cudaErrchk(cudaFreeHost(PsiProbeInit_ph));
		cudaErrchk(cudaFreeHost(trans_ph));
		cudaErrchk(cudaFreeHost(prop_ph));
		cudaErrchk(cudaFreeHost(qxa_ph));
		cudaErrchk(cudaFreeHost(qya_ph));
		cudaErrchk(cudaFreeHost(alphaInd_ph));
		for (auto s =0; s < total_num_streams; ++s){
			cudaErrchk(cudaFreeHost(output_ph[s]));
		}

		// destroy CUDA streams
		for (auto j = 0; j < total_num_streams; ++j){
			cudaSetDevice(j % pars.meta.NUM_GPUS);
			cudaErrchk(cudaStreamDestroy(streams[j]));
		}
		for (auto j = 0; j < pars.meta.NUM_GPUS; ++j) {
			cudaErrchk(cudaSetDevice(j));
			cudaErrchk(cudaDeviceReset());
		}

//		// make sure the copy is finished
//		copy_t.join();
	}
}
