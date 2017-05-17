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

#include "Multislice.cuh"
#include "Multislice.h"
#include "cuComplex.h"
#include "cufft.h"
#include "WorkDispatcher.h"
#include <iostream>
#include "fftw3.h"
#include "utility.h"
#include "utility.cuh"

namespace PRISM{
	extern std::mutex fftw_plan_lock;
	// computes the result of probe position ay,ax using the GPU. The effect of this function is the same as getMultisliceProbe_CPU
	__host__ void getMultisliceProbe_GPU_singlexfer(Parameters<PRISM_FLOAT_PRECISION>& pars,
	                                                PRISM_CUDA_COMPLEX_FLOAT* trans_d,
	                                                PRISM_CUDA_COMPLEX_FLOAT* PsiProbeInit_d,
	                                                PRISM_CUDA_COMPLEX_FLOAT* psi_ds,
	                                                PRISM_FLOAT_PRECISION* output_ph,
	                                                PRISM_FLOAT_PRECISION* psi_intensity_ds,
	                                                PRISM_FLOAT_PRECISION* integratedOutput_ds,
									                const PRISM_FLOAT_PRECISION* qya_d,
									                const PRISM_FLOAT_PRECISION* qxa_d,
									                const PRISM_CUDA_COMPLEX_FLOAT* prop_d,
									                const size_t ay,
									                const size_t ax,
									                const size_t dimj,
									                const size_t dimi,
									                const PRISM_FLOAT_PRECISION* alphaInd_d,
									                const cufftHandle& plan,
									                cudaStream_t& stream){

		// initialize psi
		PRISM_FLOAT_PRECISION yp = pars.yp[ay];
		PRISM_FLOAT_PRECISION xp = pars.xp[ax];
		const size_t N = dimj*dimi;
		const size_t num_blocks = std::min(pars.target_num_blocks, (N - 1) / BLOCK_SIZE1D + 1);
		initializePsi<<<num_blocks,BLOCK_SIZE1D, 0, stream>>>(psi_ds, PsiProbeInit_d, qya_d, qxa_d, N, yp, xp);
		for (auto planeNum = 0; planeNum < pars.numPlanes; ++planeNum) {
			cufftErrchk(PRISM_CUFFT_EXECUTE(plan, &psi_ds[0], &psi_ds[0], CUFFT_INVERSE));
			multiply_inplace<<<num_blocks,BLOCK_SIZE1D, 0, stream>>>(psi_ds, &trans_d[planeNum*N], N);
			cufftErrchk(PRISM_CUFFT_EXECUTE(plan, &psi_ds[0], &psi_ds[0], CUFFT_FORWARD));
			multiply_inplace<<<num_blocks,BLOCK_SIZE1D, 0, stream>>>(psi_ds, prop_d, N);
			divide_inplace<<<num_blocks,BLOCK_SIZE1D, 0, stream>>>(psi_ds, PRISM_MAKE_CU_COMPLEX(N, 0), N);
		}
		abs_squared<<<num_blocks,BLOCK_SIZE1D, 0, stream>>>(psi_intensity_ds, psi_ds, N);
		formatOutput_GPU_integrate(pars, psi_intensity_ds, alphaInd_d, output_ph, integratedOutput_ds, ay, ax, dimj, dimi,stream);
}

	__host__ void getMultisliceProbe_GPU_singlexfer_batch(Parameters<PRISM_FLOAT_PRECISION>& pars,
	                                                      PRISM_CUDA_COMPLEX_FLOAT* trans_d,
	                                                      PRISM_CUDA_COMPLEX_FLOAT* PsiProbeInit_d,
	                                                      PRISM_CUDA_COMPLEX_FLOAT* psi_ds,
	                                                      PRISM_FLOAT_PRECISION* output_ph,
	                                                      PRISM_FLOAT_PRECISION* psi_intensity_ds,
	                                                      PRISM_FLOAT_PRECISION* integratedOutput_ds,
	                                                      const PRISM_FLOAT_PRECISION* qya_d,
	                                                      const PRISM_FLOAT_PRECISION* qxa_d,
	                                                      const PRISM_CUDA_COMPLEX_FLOAT* prop_d,
	                                                      const size_t Nstart,
	                                                      const size_t Nstop,
	                                                      const size_t dimj,
	                                                      const size_t dimi,
	                                                      const PRISM_FLOAT_PRECISION* alphaInd_d,
	                                                      const cufftHandle& plan,
	                                                      cudaStream_t& stream){
		const size_t N = dimj*dimi;
		const size_t num_blocks = std::min(pars.target_num_blocks, (N - 1) / BLOCK_SIZE1D + 1);
		for (auto batch_idx = 0; batch_idx < (Nstop-Nstart); ++batch_idx) {
			const size_t ay = (Nstart + batch_idx) / pars.xp.size();
			const size_t ax = (Nstart + batch_idx) % pars.xp.size();

			// initialize psi
			PRISM_FLOAT_PRECISION yp = pars.yp[ay];
			PRISM_FLOAT_PRECISION xp = pars.xp[ax];

//			initializePsi << < (N - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >> >
			initializePsi << < num_blocks, BLOCK_SIZE1D, 0, stream >> >
                (psi_ds + (batch_idx * N), PsiProbeInit_d, qya_d, qxa_d, N, yp, xp);
		}
		for (auto planeNum = 0; planeNum < pars.numPlanes; ++planeNum) {
			cufftErrchk(PRISM_CUFFT_EXECUTE(plan, &psi_ds[0], &psi_ds[0], CUFFT_INVERSE));
			for (auto batch_idx = 0; batch_idx < (Nstop-Nstart); ++batch_idx) {
				multiply_inplace << < num_blocks, BLOCK_SIZE1D, 0, stream >> >
						(psi_ds + (batch_idx * N), &trans_d[planeNum * N], N);
			}
			cufftErrchk(PRISM_CUFFT_EXECUTE(plan, &psi_ds[0], &psi_ds[0], CUFFT_FORWARD));
			for (auto batch_idx = 0; batch_idx < (Nstop-Nstart); ++batch_idx) {
				multiply_inplace << < num_blocks, BLOCK_SIZE1D, 0, stream >> >
						(psi_ds + (batch_idx * N), prop_d, N);
				divide_inplace << < num_blocks, BLOCK_SIZE1D, 0, stream >> >
						(psi_ds + (batch_idx * N), PRISM_MAKE_CU_COMPLEX(N, 0), N);
			}
		}
		abs_squared << < num_blocks, BLOCK_SIZE1D, 0, stream >> > (psi_intensity_ds, psi_ds, N*(Nstop-Nstart));
//		if (Nstart==0 ) {
//			cout << "after propagation" << endl;
//			{
////				std::complex<PRISM_FLOAT_PRECISION> ans_cx;
//				PRISM_FLOAT_PRECISION ans;
//				PRISM_FLOAT_PRECISION psi_intensity_sum = 0;
//				for (auto batch_idx = 0; batch_idx < (Nstop-Nstart); ++batch_idx) {
//					psi_intensity_sum = 0;
//					cout << "BATCH ID = " << batch_idx << endl;
//					for (auto i = 0; i < 10; ++i) {
//						cudaErrchk(cudaMemcpy(&ans, psi_intensity_ds + i + ((batch_idx * N)), sizeof(ans), cudaMemcpyDeviceToHost));
//						cout << "psi_intensity_ds[" << i  + (batch_idx * N) << "] = " << ans << endl;
//					}
//					for (auto i = 0; i < pars.psiProbeInit.size(); ++i) {
//						cudaErrchk(cudaMemcpy(&ans, psi_intensity_ds + i, sizeof(ans), cudaMemcpyDeviceToHost));
//						psi_intensity_sum += ans;
//					}
//					Array2D<PRISM_FLOAT_PRECISION> debug = zeros_ND<2, PRISM_FLOAT_PRECISION >({{pars.psiProbeInit.get_dimj(), pars.psiProbeInit.get_dimi()}});
//					cudaErrchk(cudaMemcpy(&debug[0], psi_intensity_ds, pars.psiProbeInit.size() * sizeof(debug[0]), cudaMemcpyDeviceToHost));
//					debug.toMRC_f("p1.mrc");
//					cudaErrchk(cudaMemcpy(&debug[0], psi_intensity_ds + pars.psiProbeInit.size(), pars.psiProbeInit.size() * sizeof(debug[0]), cudaMemcpyDeviceToHost));
//					debug.toMRC_f("p2.mrc");
//					cout << "PSI_INTENSITY SUM = " << psi_intensity_sum << endl;
//				}
//			}
//		}

		for (auto batch_idx = 0; batch_idx < (Nstop-Nstart); ++batch_idx) {
			const size_t ay = (Nstart + batch_idx) / pars.xp.size();
			const size_t ax = (Nstart + batch_idx) % pars.xp.size();
			formatOutput_GPU_integrate(pars, psi_intensity_ds + (batch_idx * N),
			                           alphaInd_d, output_ph, integratedOutput_ds, ay, ax, dimj, dimi, stream);
		}
	}




	__host__ void getMultisliceProbe_GPU_streaming(Parameters<PRISM_FLOAT_PRECISION>& pars,
	                                               PRISM_CUDA_COMPLEX_FLOAT* trans_d,
	                                               const complex<PRISM_FLOAT_PRECISION>* trans_ph,
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
		const size_t num_blocks = std::min(pars.target_num_blocks, (N - 1) / BLOCK_SIZE1D + 1);
		initializePsi<<<num_blocks,BLOCK_SIZE1D, 0, stream>>>(psi_ds, PsiProbeInit_d, qya_d, qxa_d, dimj*dimi, yp, xp);


		for (auto planeNum = 0; planeNum < pars.numPlanes; ++planeNum) {
			cudaErrchk(cudaMemcpyAsync(trans_d, &trans_ph[planeNum*N], N * sizeof(PRISM_CUDA_COMPLEX_FLOAT), cudaMemcpyHostToDevice, stream));
			cufftErrchk(PRISM_CUFFT_EXECUTE(plan, &psi_ds[0], &psi_ds[0], CUFFT_INVERSE));
			multiply_inplace<<<num_blocks,BLOCK_SIZE1D, 0, stream>>>(psi_ds, trans_d, N);
			cufftErrchk(PRISM_CUFFT_EXECUTE(plan, &psi_ds[0], &psi_ds[0], CUFFT_FORWARD));
			multiply_inplace<<<num_blocks,BLOCK_SIZE1D, 0, stream>>>(psi_ds, prop_d, N);
			divide_inplace<<<num_blocks,BLOCK_SIZE1D, 0, stream>>>(psi_ds, PRISM_MAKE_CU_COMPLEX(N, 0), N);
		}
		abs_squared<<<num_blocks,BLOCK_SIZE1D, 0, stream>>>(psi_intensity_ds, psi_ds, N);
		formatOutput_GPU_integrate(pars, psi_intensity_ds, alphaInd_d, output_ph, integratedOutput_ds, ay, ax, dimj, dimi,stream);
	}


	__host__ void getMultisliceProbe_GPU_streaming_batch(Parameters<PRISM_FLOAT_PRECISION>& pars,
	                                                     PRISM_CUDA_COMPLEX_FLOAT* trans_d,
	                                                     const complex<PRISM_FLOAT_PRECISION>* trans_ph,
	                                                     PRISM_CUDA_COMPLEX_FLOAT* PsiProbeInit_d,
	                                                     PRISM_CUDA_COMPLEX_FLOAT* psi_ds,
	                                                     PRISM_FLOAT_PRECISION* output_ph,
	                                                     PRISM_FLOAT_PRECISION* psi_intensity_ds,
	                                                     PRISM_FLOAT_PRECISION* integratedOutput_ds,
	                                                     const PRISM_FLOAT_PRECISION* qya_d,
	                                                     const PRISM_FLOAT_PRECISION* qxa_d,
	                                                     const PRISM_CUDA_COMPLEX_FLOAT* prop_d,
	                                                     const size_t Nstart,
	                                                     const size_t Nstop,
	                                                     const size_t dimj,
	                                                     const size_t dimi,
	                                                     const PRISM_FLOAT_PRECISION* alphaInd_d,
	                                                     const cufftHandle& plan,
	                                                     cudaStream_t& stream){

		// initialize psi
		const size_t N = dimj*dimi;
		const size_t num_blocks = std::min(pars.target_num_blocks, (N - 1) / BLOCK_SIZE1D + 1);
		for (auto batch_idx = 0; batch_idx < (Nstop-Nstart); ++batch_idx) {
			const size_t ay = (Nstart + batch_idx) / pars.xp.size();
			const size_t ax = (Nstart + batch_idx) % pars.xp.size();
			PRISM_FLOAT_PRECISION yp = pars.yp[ay];
			PRISM_FLOAT_PRECISION xp = pars.xp[ax];
			initializePsi << < num_blocks, BLOCK_SIZE1D, 0, stream >> >
			                                                (psi_ds + (batch_idx * N), PsiProbeInit_d, qya_d, qxa_d, N, yp, xp);
		}

		for (auto planeNum = 0; planeNum < pars.numPlanes; ++planeNum) {

			cudaErrchk(cudaMemcpyAsync(trans_d, &trans_ph[planeNum*N], N * sizeof(PRISM_CUDA_COMPLEX_FLOAT), cudaMemcpyHostToDevice, stream));
			cufftErrchk(PRISM_CUFFT_EXECUTE(plan, &psi_ds[0], &psi_ds[0], CUFFT_INVERSE));
			for (auto batch_idx = 0; batch_idx < (Nstop-Nstart); ++batch_idx) {
				multiply_inplace << < num_blocks, BLOCK_SIZE1D, 0, stream >> >
				                                                   (psi_ds + (batch_idx * N), trans_d, N);
			}
			cufftErrchk(PRISM_CUFFT_EXECUTE(plan, &psi_ds[0], &psi_ds[0], CUFFT_FORWARD));
			for (auto batch_idx = 0; batch_idx < (Nstop-Nstart); ++batch_idx) {
				multiply_inplace << < num_blocks, BLOCK_SIZE1D, 0, stream >> >
				                                                   (psi_ds + (batch_idx * N), prop_d, N);
				divide_inplace << < num_blocks, BLOCK_SIZE1D, 0, stream >> >
				                                                 (psi_ds + (batch_idx * N), PRISM_MAKE_CU_COMPLEX(N, 0), N);
			}
		}
		abs_squared << < num_blocks, BLOCK_SIZE1D, 0, stream >> > (psi_intensity_ds, psi_ds, N*(Nstop-Nstart));
		for (auto batch_idx = 0; batch_idx < (Nstop-Nstart); ++batch_idx) {
			const size_t ay = (Nstart + batch_idx) / pars.xp.size();
			const size_t ax = (Nstart + batch_idx) % pars.xp.size();
			formatOutput_GPU_integrate(pars, psi_intensity_ds + (batch_idx * N),
			                           alphaInd_d, output_ph, integratedOutput_ds, ay, ax, dimj, dimi, stream);
		}
	}

    __host__ void buildMultisliceOutput_GPU_singlexfer(Parameters <PRISM_FLOAT_PRECISION> &pars){
#ifdef PRISM_BUILDING_GUI
	    pars.progressbar->signalDescriptionMessage("Computing final output");
#endif
		// determine the batch size to use
	    pars.meta.batch_size_GPU = min(pars.meta.batch_size_target_GPU, pars.xp.size() * pars.yp.size() / (pars.meta.NUM_STREAMS_PER_GPU*pars.meta.NUM_GPUS));
		// populate the Multislice output stack dividing the work between GPUs and CPU cores.
		// this version assumes the full trans array fits into DRAM on each GPU
		using namespace std;
		// create CUDA streams
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;
        cudaStream_t *streams   = new cudaStream_t[total_num_streams];
        cufftHandle *cufft_plan = new cufftHandle[total_num_streams];
//		cudaStream_t streams[total_num_streams];
//		cufftHandle cufft_plan[total_num_streams];


	    // batch parameters for cuFFT
	    const int rank = 2;
	    int n[] = {(int)pars.psiProbeInit.get_dimj(), (int)pars.psiProbeInit.get_dimi()};
	    const int howmany = pars.meta.batch_size_GPU;
	    cout <<"pars.meta.batch_size_GPU= " << pars.meta.batch_size_GPU<< endl;
	    int idist = n[0]*n[1];
	    int odist = n[0]*n[1];
	    int istride = 1;
	    int ostride = 1;
	    int *inembed = n;
	    int *onembed = n;

		for (auto j = 0; j < total_num_streams; ++j){
			cudaSetDevice(j % pars.meta.NUM_GPUS);
			cudaErrchk(cudaStreamCreate(&streams[j]));
			cufftErrchk(cufftPlanMany(&cufft_plan[j], rank, n, inembed, istride, idist, onembed, ostride, odist, PRISM_CUFFT_PLAN_TYPE, howmany));
//			cufftErrchk(cufftPlan2d(&cufft_plan[j], pars.psiProbeInit.get_dimj(), pars.psiProbeInit.get_dimi(), PRISM_CUFFT_PLAN_TYPE));
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
//		PRISM_FLOAT_PRECISION                *output_ph[total_num_streams];
		PRISM_FLOAT_PRECISION                **output_ph = new PRISM_FLOAT_PRECISION*[total_num_streams];
		// allocate pinned memory
		cudaErrchk(cudaMallocHost((void **)&PsiProbeInit_ph, pars.psiProbeInit.size()*sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
		cudaErrchk(cudaMallocHost((void **)&trans_ph,        pars.transmission.size()*sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
		cudaErrchk(cudaMallocHost((void **)&prop_ph,         pars.prop.size()*sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
		cudaErrchk(cudaMallocHost((void **)&qxa_ph,          pars.qxa.size()*sizeof(PRISM_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **)&qya_ph,          pars.qya.size()*sizeof(PRISM_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **)&alphaInd_ph,     pars.alphaInd.size()*sizeof(PRISM_FLOAT_PRECISION)));
		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaMallocHost((void **) &output_ph[s], pars.output.get_dimi() * sizeof(PRISM_FLOAT_PRECISION)));
		}
		// copy host memory to pinned
		memcpy(PsiProbeInit_ph, &pars.psiProbeInit[0], pars.psiProbeInit.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>));
		memcpy(trans_ph,        &pars.transmission[0],        pars.transmission.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>));
		memcpy(prop_ph,         &pars.prop[0],    pars.prop.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>));
		memcpy(qxa_ph,          &pars.qxa[0],     pars.qxa.size() * sizeof(PRISM_FLOAT_PRECISION));
		memcpy(qya_ph,          &pars.qya[0],     pars.qya.size() * sizeof(PRISM_FLOAT_PRECISION));
		memcpy(alphaInd_ph,     &pars.alphaInd[0],     pars.alphaInd.size() * sizeof(PRISM_FLOAT_PRECISION));


		// pointers to read-only GPU memory (one copy per GPU)
		PRISM_CUDA_COMPLEX_FLOAT **PsiProbeInit_d = new PRISM_CUDA_COMPLEX_FLOAT*[pars.meta.NUM_GPUS];
		PRISM_CUDA_COMPLEX_FLOAT **trans_d		  = new PRISM_CUDA_COMPLEX_FLOAT*[pars.meta.NUM_GPUS];
		PRISM_CUDA_COMPLEX_FLOAT **prop_d 		  = new PRISM_CUDA_COMPLEX_FLOAT*[pars.meta.NUM_GPUS];
		PRISM_FLOAT_PRECISION    **qxa_d 		  = new PRISM_FLOAT_PRECISION*[pars.meta.NUM_GPUS];
		PRISM_FLOAT_PRECISION    **qya_d 		  = new PRISM_FLOAT_PRECISION*[pars.meta.NUM_GPUS];
		PRISM_FLOAT_PRECISION    **alphaInd_d     = new PRISM_FLOAT_PRECISION*[pars.meta.NUM_GPUS];

		// pointers to read/write GPU memory (one per stream)
		PRISM_CUDA_COMPLEX_FLOAT **psi_ds 			   = new PRISM_CUDA_COMPLEX_FLOAT*[total_num_streams];
		PRISM_FLOAT_PRECISION    **psi_intensity_ds    = new PRISM_FLOAT_PRECISION*[total_num_streams];
		PRISM_FLOAT_PRECISION    **integratedOutput_ds = new PRISM_FLOAT_PRECISION*[total_num_streams];
//		PRISM_CUDA_COMPLEX_FLOAT *PsiProbeInit_d[pars.meta.NUM_GPUS];
//		PRISM_CUDA_COMPLEX_FLOAT *trans_d[pars.meta.NUM_GPUS];
//		PRISM_CUDA_COMPLEX_FLOAT *prop_d[pars.meta.NUM_GPUS];
//		PRISM_FLOAT_PRECISION    *qxa_d[pars.meta.NUM_GPUS];
//		PRISM_FLOAT_PRECISION    *qya_d[pars.meta.NUM_GPUS];
//	    PRISM_FLOAT_PRECISION    *alphaInd_d[pars.meta.NUM_GPUS];
//
//		// pointers to read/write GPU memory (one per stream)
//		PRISM_CUDA_COMPLEX_FLOAT *psi_ds[total_num_streams];
//		PRISM_FLOAT_PRECISION    *psi_intensity_ds[total_num_streams];
//		PRISM_FLOAT_PRECISION    *integratedOutput_ds[total_num_streams];

		// allocate memory on each GPU
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaMalloc((void **) &PsiProbeInit_d[g],     pars.psiProbeInit.size()   * sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
			cudaErrchk(cudaMalloc((void **) &trans_d[g],            pars.transmission.size()   * sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
			cudaErrchk(cudaMalloc((void **) &prop_d[g],             pars.prop.size()           * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMalloc((void **) &qxa_d[g],              pars.qxa.size()            * sizeof(PRISM_FLOAT_PRECISION)));
			cudaErrchk(cudaMalloc((void **) &qya_d[g],              pars.qya.size()            * sizeof(PRISM_FLOAT_PRECISION)));
			cudaErrchk(cudaMalloc((void **) &alphaInd_d[g],         pars.alphaInd.size()       * sizeof(PRISM_FLOAT_PRECISION)));
		}

	    cout << "pars.psiProbeInit.get_dimj()  = " << pars.psiProbeInit.get_dimj()<< endl;
	    cout << "pars.psiProbeInit.get_dimj() = " << pars.psiProbeInit.get_dimi()<< endl;
	    cout << "pars.psiProbeInit.size()  = " << pars.psiProbeInit.size()<< endl;
	    cout << "pars.meta.batch_size_GPU  = " << pars.meta.batch_size_GPU<< endl;
	    cout << "pars.meta.batch_size_GPU*pars.psiProbeInit.size() * sizeof(pars.psiProbeInit[0]) = " << pars.meta.batch_size_GPU*pars.psiProbeInit.size()        * sizeof(pars.psiProbeInit[0]) << endl;
	    cout << "pars.meta.batch_size_GPU*pars.psiProbeInit.size() *  sizeof(PRISM_FLOAT_PRECISION) = " << pars.meta.batch_size_GPU*pars.psiProbeInit.size()        * sizeof(PRISM_FLOAT_PRECISION)<< endl;

	    for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaSetDevice(s % pars.meta.NUM_GPUS));
			cudaErrchk(cudaMalloc((void **) &psi_ds[s],              pars.meta.batch_size_GPU*pars.psiProbeInit.size()        * sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
			cudaErrchk(cudaMalloc((void **) &psi_intensity_ds[s],    pars.meta.batch_size_GPU*pars.psiProbeInit.size()        * sizeof(PRISM_FLOAT_PRECISION)));
			cudaErrchk(cudaMalloc((void **) &integratedOutput_ds[s], pars.detectorAngles.size() * sizeof(PRISM_FLOAT_PRECISION)));
			cudaErrchk(cudaMemset(psi_ds[s], 0,                      pars.meta.batch_size_GPU*pars.psiProbeInit.size()  * sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
			cudaErrchk(cudaMemset(psi_intensity_ds[s], 0,            pars.meta.batch_size_GPU*pars.psiProbeInit.size()        * sizeof(PRISM_FLOAT_PRECISION)));
			cudaErrchk(cudaMemset(integratedOutput_ds[s], 0,         pars.detectorAngles.size() * sizeof(PRISM_FLOAT_PRECISION)));
		}


		// Copy memory to each GPU asynchronously from the pinned host memory spaces.
		// The streams are laid out so that consecutive streams represent different GPUs. If we
		// have more than one stream per GPU, then we want to interleave as much as possible
		int stream_id = 0;
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			stream_id = g;
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaMemcpyAsync(PsiProbeInit_d[g], &PsiProbeInit_ph[0],
			                      pars.psiProbeInit.size() * sizeof(pars.psiProbeInit[0]), cudaMemcpyHostToDevice, streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(trans_d[g], &trans_ph[0],
			                      pars.transmission.size() * sizeof(pars.transmission[0]), cudaMemcpyHostToDevice, streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(prop_d[g], &prop_ph[0],
			                      pars.prop.size() * sizeof(pars.prop[0]), cudaMemcpyHostToDevice, streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(qxa_d[g], &qxa_ph[0],
			                      pars.qxa.size() * sizeof(pars.qxa[0]), cudaMemcpyHostToDevice, streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(qya_d[g], &qya_ph[0],
			                      pars.qya.size() * sizeof(pars.qya[0]), cudaMemcpyHostToDevice, streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(alphaInd_d[g], &alphaInd_ph[0],
			                      pars.alphaInd.size() * sizeof(pars.alphaInd[0]), cudaMemcpyHostToDevice, streams[stream_id]));
		}

		// make sure transfers are complete
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g){
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaDeviceSynchronize());
		}

		size_t psi_size = pars.psiProbeInit.size();
		int stream_count = 0;
//		setWorkStartStop(0, pars.xp.size() * pars.yp.size());
	    const size_t PRISM_PRINT_FREQUENCY_PROBES = pars.xp.size() * pars.yp.size() / 10; // for printing status
		WorkDispatcher dispatcher(0, pars.xp.size() * pars.yp.size());
//		setWorkStartStop(0, 1);
		cout << " pars.xp.size()  = " << pars.xp.size()  << endl;
		cout << " pars.yp.size()  = " << pars.yp.size()  << endl;

		for (auto t = 0; t < total_num_streams; ++t){
			int GPU_num = stream_count % pars.meta.NUM_GPUS; // determine which GPU handles this job
			cudaStream_t& current_stream = streams[stream_count];
			cout << "Launching GPU worker on stream #" << stream_count << " on GPU #" << GPU_num << '\n';

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
			workers_GPU.push_back(thread([&pars, current_trans_d, current_PsiProbeInit_d, current_alphaInd_d, &dispatcher,
					                                current_psi_ds, current_psi_intensity_ds, current_integratedOutput_ds,
					                                GPU_num, current_qya_d, current_qxa_d, current_output_ph, &current_cufft_plan,
					                                current_prop_d, &current_stream, &psi_size, stream_count, &PRISM_PRINT_FREQUENCY_PROBES]() {

				// set the GPU context
				cudaErrchk(cudaSetDevice(GPU_num)); // set current GPU

#ifndef NDEBUG
				{
//					 check memory usage on the GPU
					std::lock_guard<mutex> lock(PRISM::mem_lock);
					size_t free_mem, total_mem;
					cudaErrchk(cudaMemGetInfo(&free_mem, &total_mem));
					pars.max_mem = std::max(total_mem - free_mem, pars.max_mem);
//					cout << "max_mem = " << pars.max_mem << '\n';
				}
#endif // NDEBUG

				size_t Nstart, Nstop;
				Nstart=Nstop=0;
//				while (getWorkID(pars, Nstart, Nstop)){ // synchronously get work assignment
				while (dispatcher.getWork(Nstart, Nstop, pars.meta.batch_size_GPU)){ // synchronously get work assignment
					while (Nstart < Nstop){
						if (Nstart % PRISM_PRINT_FREQUENCY_PROBES < pars.meta.batch_size_GPU | Nstart == 100){
							cout << "Computing Probe Position #" << Nstart << "/" << pars.xp.size() * pars.yp.size() << '\n';
						}
//						ay = Nstart / pars.xp.size();
//						ax = Nstart % pars.xp.size();
//						cout << "outside ax = " << ax << endl;
//						cout << "outside ay = " << ay <<vi  endl;
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
#ifdef PRISM_BUILDING_GUI
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
			PRISM_FFTW_INIT_THREADS();
			PRISM_FFTW_PLAN_WITH_NTHREADS(pars.meta.NUM_THREADS);vector<thread> workers_CPU;
			workers_CPU.reserve(pars.meta.NUM_THREADS); // prevents multiple reallocations
			// If the batch size is too big, the work won't be spread over the threads, which will usually hurt more than the benefit
			// of batch FFT
			pars.meta.batch_size_CPU = min(pars.meta.batch_size_target_CPU, pars.xp.size() * pars.yp.size() / pars.meta.NUM_THREADS);
			cout << "multislice pars.meta.batch_size_CPU = " << pars.meta.batch_size_CPU << endl;
			for (auto t = 0; t < pars.meta.NUM_THREADS; ++t) {
				cout << "Launching CPU worker #" << t << endl;
				// push_back is better whenever constructing a new object
				workers_CPU.push_back(thread([&pars, &dispatcher, t, &PRISM_PRINT_FREQUENCY_PROBES]() {
				size_t Nstart, Nstop, early_CPU_stop;
				Nstart=Nstop=0;
				// stop the CPU workers earlier than the GPU ones to prevent slower workers taking the last jobs and having to
				// wait longer for everything to complete
                    if (pars.meta.NUM_GPUS > 0){
                      // if there are no GPUs, make sure to do all work on CPU
                        early_CPU_stop = (size_t)std::max((PRISM_FLOAT_PRECISION)0.0, pars.xp.size() * pars.yp.size() - pars.meta.gpu_cpu_ratio);
                    } else {
                        early_CPU_stop = pars.xp.size() * pars.yp.size();
                    }
					if (dispatcher.getWork(Nstart, Nstop, pars.meta.batch_size_CPU, early_CPU_stop)) { // synchronously get work assignment
                        Array1D<std::complex<PRISM_FLOAT_PRECISION> > psi_stack = zeros_ND<1, complex<PRISM_FLOAT_PRECISION> >({{pars.psiProbeInit.size() * pars.meta.batch_size_CPU}});

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

                        gatekeeper.unlock();
						do {
							//	cout << "Nstop = " << Nstop << endl;
							while (Nstart < Nstop) {
                                if (Nstart % PRISM_PRINT_FREQUENCY_PROBES  < pars.meta.batch_size_CPU | Nstart == 100){
                                    cout << "Computing Probe Position #" << Nstart << "/" << pars.xp.size() * pars.yp.size() << endl;
                                }
//							ay = Nstart / pars.xp.size();
//							ax = Nstart % pars.xp.size();
//                            if (ay==7){
//                                cout << "ax = " << ax << endl;
//								cout << "ay = " << ay << endl;
//                            }
//							getMultisliceProbe_CPU(pars, ay, ax, plan_forward, plan_inverse, psi);
                                getMultisliceProbe_CPU_batch(pars, Nstart, Nstop, plan_forward, plan_inverse, psi_stack);
#ifdef PRISM_BUILDING_GUI
                                pars.progressbar->signalOutputUpdate(Nstart, pars.xp.size() * pars.yp.size());
#endif
                                //++Nstart;
                                Nstart=Nstop;
							}
							if (Nstop >= early_CPU_stop) break;
						} while(dispatcher.getWork(Nstart, Nstop, pars.meta.batch_size_CPU, early_CPU_stop));
						gatekeeper.lock();
						PRISM_FFTW_DESTROY_PLAN(plan_forward);
						PRISM_FFTW_DESTROY_PLAN(plan_inverse);
						gatekeeper.unlock();
					}
//					cout << "CPU worker #" << t << " finished\n";
			
					}));
				
			}
			cout << "Waiting on CPU threads..." << endl;
			for (auto& t:workers_CPU)t.join();
			PRISM_FFTW_CLEANUP_THREADS();
		}
		// synchronize threads
		cout << "Waiting on GPU threads..." << endl;
		for (auto& t:workers_GPU)t.join();



		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g){
			cudaSetDevice(g);
			cudaDeviceSynchronize();
		}

		// copy the results of the GPU, which are in pinned memory, back to the actual stack. The CPU work populates the
		// beginning, so make sure to copy from the offset of where the GPU started. Launch this copy on a background thread
		// while we cleanup the GPU
//		const size_t GPU_start_offset = (size_t)CPU_stop*pars.output.get_dimk()*pars.output.get_dimj()*pars.output.get_dimi();
//		std::thread copy_t([&GPU_start_offset, &pars, &stack_ph](){
//			memcpy(&pars.output[GPU_start_offset],
//			       &stack_ph[GPU_start_offset],
//			       (pars.output.size()-GPU_start_offset) * sizeof(PRISM_FLOAT_PRECISION));
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
		delete[] streams;
		delete[] cufft_plan;
		delete[] PsiProbeInit_d;
		delete[] trans_d;
		delete[] prop_d;
		delete[] qxa_d;
		delete[] qya_d;
		delete[] alphaInd_d;
		delete[] psi_ds;
		delete[] psi_intensity_ds;
		delete[] integratedOutput_ds;
		delete[] output_ph;
	}




	__host__ void buildMultisliceOutput_GPU_streaming(Parameters <PRISM_FLOAT_PRECISION> &pars){
#ifdef PRISM_BUILDING_GUI
		pars.progressbar->signalDescriptionMessage("Computing final output");
#endif

		// populate the Multislice output stack dividing the work between GPUs and CPU cores.
		// this version assumes the full trans array fits into DRAM on each GPU

		using namespace std;

		// determine batch size
		pars.meta.batch_size_GPU = min(pars.meta.batch_size_target_GPU, pars.xp.size() * pars.yp.size() / (pars.meta.NUM_STREAMS_PER_GPU*pars.meta.NUM_GPUS));
		cout << "multislice pars.meta.batch_size_CPU = " << pars.meta.batch_size_CPU << endl;
		// create CUDA streams
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;
        cudaStream_t *streams   = new cudaStream_t[total_num_streams];
        cufftHandle *cufft_plan = new cufftHandle[total_num_streams];
		cout <<"total_num_streams = " << total_num_streams<< endl;
		for (auto j = 0; j < total_num_streams; ++j){
			cudaSetDevice(j % pars.meta.NUM_GPUS);
			cudaErrchk(cudaStreamCreate(&streams[j]));
			cufftErrchk(cufftPlan2d(&cufft_plan[j], pars.psiProbeInit.get_dimj(), pars.psiProbeInit.get_dimi(), PRISM_CUFFT_PLAN_TYPE));
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
//		PRISM_FLOAT_PRECISION                *output_ph[total_num_streams];
		PRISM_FLOAT_PRECISION                **output_ph = new PRISM_FLOAT_PRECISION*[total_num_streams];
		// allocate pinned memory
		cudaErrchk(cudaMallocHost((void **)&PsiProbeInit_ph, pars.psiProbeInit.size()*sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
		cudaErrchk(cudaMallocHost((void **)&trans_ph,        pars.transmission.size()*sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
		cudaErrchk(cudaMallocHost((void **)&prop_ph,         pars.prop.size()*sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
		cudaErrchk(cudaMallocHost((void **)&qxa_ph,          pars.qxa.size()*sizeof(PRISM_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **)&qya_ph,          pars.qya.size()*sizeof(PRISM_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **)&alphaInd_ph,     pars.alphaInd.size()*sizeof(PRISM_FLOAT_PRECISION)));
		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaMallocHost((void **) &output_ph[s], pars.output.get_dimi() * sizeof(PRISM_FLOAT_PRECISION)));
		}
		// copy host memory to pinned
		memcpy(PsiProbeInit_ph, &pars.psiProbeInit[0], pars.psiProbeInit.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>));
		memcpy(trans_ph,        &pars.transmission[0],        pars.transmission.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>));
		memcpy(prop_ph,         &pars.prop[0],    pars.prop.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>));
		memcpy(qxa_ph,          &pars.qxa[0],     pars.qxa.size() * sizeof(PRISM_FLOAT_PRECISION));
		memcpy(qya_ph,          &pars.qya[0],     pars.qya.size() * sizeof(PRISM_FLOAT_PRECISION));
		memcpy(alphaInd_ph,     &pars.alphaInd[0],     pars.alphaInd.size() * sizeof(PRISM_FLOAT_PRECISION));


		// pointers to read-only GPU memory (one copy per GPU)
		PRISM_CUDA_COMPLEX_FLOAT **PsiProbeInit_d = new PRISM_CUDA_COMPLEX_FLOAT*[pars.meta.NUM_GPUS];
		PRISM_CUDA_COMPLEX_FLOAT **prop_d 	   	  = new PRISM_CUDA_COMPLEX_FLOAT*[pars.meta.NUM_GPUS];
		PRISM_FLOAT_PRECISION    **qxa_d 		  = new PRISM_FLOAT_PRECISION*[pars.meta.NUM_GPUS];
		PRISM_FLOAT_PRECISION    **qya_d 		  = new PRISM_FLOAT_PRECISION*[pars.meta.NUM_GPUS];
		PRISM_FLOAT_PRECISION    **alphaInd_d 	  = new PRISM_FLOAT_PRECISION*[pars.meta.NUM_GPUS];

		// pointers to read/write GPU memory (one per stream)
		PRISM_CUDA_COMPLEX_FLOAT **trans_ds 		   = new PRISM_CUDA_COMPLEX_FLOAT*[total_num_streams];
		PRISM_CUDA_COMPLEX_FLOAT **psi_ds  		       = new PRISM_CUDA_COMPLEX_FLOAT*[total_num_streams];
		PRISM_FLOAT_PRECISION    **psi_intensity_ds    = new PRISM_FLOAT_PRECISION*[total_num_streams];
		PRISM_FLOAT_PRECISION    **integratedOutput_ds = new PRISM_FLOAT_PRECISION*[total_num_streams];
//		// pointers to read-only GPU memory (one copy per GPU)
//		PRISM_CUDA_COMPLEX_FLOAT *PsiProbeInit_d[pars.meta.NUM_GPUS];
//		PRISM_CUDA_COMPLEX_FLOAT *prop_d[pars.meta.NUM_GPUS];
//		PRISM_FLOAT_PRECISION    *qxa_d[pars.meta.NUM_GPUS];
//		PRISM_FLOAT_PRECISION    *qya_d[pars.meta.NUM_GPUS];
//		PRISM_FLOAT_PRECISION    *alphaInd_d[pars.meta.NUM_GPUS];
//
//		// pointers to read/write GPU memory (one per stream)
//		PRISM_CUDA_COMPLEX_FLOAT *trans_ds[total_num_streams];
//		PRISM_CUDA_COMPLEX_FLOAT *psi_ds[total_num_streams];
//		PRISM_FLOAT_PRECISION    *psi_intensity_ds[total_num_streams];
//		PRISM_FLOAT_PRECISION    *integratedOutput_ds[total_num_streams];

		// allocate memory on each GPU
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaMalloc((void **) &PsiProbeInit_d[g],     pars.psiProbeInit.size()        * sizeof(pars.psiProbeInit[0])));
			cudaErrchk(cudaMalloc((void **) &prop_d[g],             pars.prop.size()           * sizeof(pars.prop[0])));
			cudaErrchk(cudaMalloc((void **) &qxa_d[g],              pars.qxa.size()            * sizeof(pars.qxa[0])));
			cudaErrchk(cudaMalloc((void **) &qya_d[g],              pars.qya.size()            * sizeof(pars.qya[0])));
			cudaErrchk(cudaMalloc((void **) &alphaInd_d[g],         pars.alphaInd.size()            * sizeof(pars.alphaInd[0])));
		}

		cout << "pars.psiProbeInit.size()  = " << pars.psiProbeInit.size()<< endl;
		cout << "pars.meta.batch_size_GPU  = " << pars.meta.batch_size_GPU<< endl;
		cout << "pars.meta.batch_size_GPU*pars.psiProbeInit.size() * sizeof(pars.psiProbeInit[0]) = " << pars.meta.batch_size_GPU*pars.psiProbeInit.size()        * sizeof(pars.psiProbeInit[0]) << endl;
		cout << "pars.meta.batch_size_GPU*pars.psiProbeInit.size() *  sizeof(PRISM_FLOAT_PRECISION) = " << pars.meta.batch_size_GPU*pars.psiProbeInit.size()        * sizeof(PRISM_FLOAT_PRECISION)<< endl;
		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaSetDevice(s % pars.meta.NUM_GPUS));
			cudaErrchk(cudaMalloc((void **) &trans_ds[s],            pars.transmission.get_dimj() * pars.transmission.get_dimi() * sizeof(pars.transmission[0])));
			cudaErrchk(cudaMalloc((void **) &psi_ds[s],              pars.meta.batch_size_GPU*pars.psiProbeInit.size()        * sizeof(pars.psiProbeInit[0])));
			cudaErrchk(cudaMalloc((void **) &psi_intensity_ds[s],    pars.meta.batch_size_GPU*pars.psiProbeInit.size()        * sizeof(PRISM_FLOAT_PRECISION)));
			cudaErrchk(cudaMalloc((void **) &integratedOutput_ds[s], pars.detectorAngles.size() * sizeof(PRISM_FLOAT_PRECISION)));
			cudaErrchk(cudaMemset(psi_ds[s], 0, pars.meta.batch_size_GPU*pars.psiProbeInit.size()        * sizeof(pars.psiProbeInit[0])));
			cudaErrchk(cudaMemset(psi_intensity_ds[s], 0, pars.meta.batch_size_GPU*pars.psiProbeInit.size()        * sizeof(PRISM_FLOAT_PRECISION)));
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
			                           pars.psiProbeInit.size() * sizeof(pars.psiProbeInit[0]), cudaMemcpyHostToDevice, streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
//			cudaErrchk(cudaMemcpyAsync(trans_d[g], &trans_ph[0],
//			                           trans.size() * sizeof(trans[0]), cudaMemcpyHostToDevice, streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(prop_d[g], &prop_ph[0],
			                           pars.prop.size() * sizeof(pars.prop[0]), cudaMemcpyHostToDevice, streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(qxa_d[g], &qxa_ph[0],
			                           pars.qxa.size() * sizeof(pars.qxa[0]), cudaMemcpyHostToDevice, streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(qya_d[g], &qya_ph[0],
			                           pars.qya.size() * sizeof(pars.qya[0]), cudaMemcpyHostToDevice, streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cudaErrchk(cudaMemcpyAsync(alphaInd_d[g], &alphaInd_ph[0],
			                           pars.alphaInd.size() * sizeof(pars.alphaInd[0]), cudaMemcpyHostToDevice, streams[stream_id]));
		}

		// make sure transfers are complete
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g){
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaDeviceSynchronize());
		}

		size_t psi_size = pars.psiProbeInit.size();
		int stream_count = 0;
//		setWorkStartStop(0, pars.xp.size() * pars.yp.size());
		const size_t PRISM_PRINT_FREQUENCY_PROBES = pars.xp.size() * pars.yp.size() / 10; // for printing status
		WorkDispatcher dispatcher(0, pars.xp.size() * pars.yp.size());
		// If the batch size is too big, the work won't be spread over the threads, which will usually hurt more than the benefit
		// of batch FFT

		for (auto t = 0; t < total_num_streams; ++t){
			int GPU_num = stream_count % pars.meta.NUM_GPUS; // determine which GPU handles this job
			cudaStream_t& current_stream = streams[stream_count];
			cout << "Launching GPU worker on stream #" << stream_count << " of GPU #" << GPU_num << endl;

			// get pointers to the pre-copied arrays, making sure to get those on the current GPU
			PRISM_CUDA_COMPLEX_FLOAT *current_PsiProbeInit_d = PsiProbeInit_d[GPU_num];
			PRISM_CUDA_COMPLEX_FLOAT *current_prop_d   = prop_d[GPU_num];
			PRISM_FLOAT_PRECISION *current_qxa_d       = qxa_d[GPU_num];
			PRISM_FLOAT_PRECISION *current_qya_d       = qya_d[GPU_num];
			PRISM_FLOAT_PRECISION *current_alphaInd_d  = alphaInd_d[GPU_num];

			// get pointers to per-stream arrays
			PRISM_CUDA_COMPLEX_FLOAT *current_trans_ds         = trans_ds[stream_count];
			PRISM_CUDA_COMPLEX_FLOAT *current_psi_ds           = psi_ds[stream_count];
			PRISM_FLOAT_PRECISION *current_psi_intensity_ds    = psi_intensity_ds[stream_count];
			PRISM_FLOAT_PRECISION *current_integratedOutput_ds = integratedOutput_ds[stream_count];
			PRISM_FLOAT_PRECISION *current_output_ph           = output_ph[stream_count];
			cufftHandle & current_cufft_plan                   = cufft_plan[stream_count];
			// launch a new thread
			// push_back is better whenever constructing a new object
			workers_GPU.push_back(thread([&pars, current_trans_ds, trans_ph, current_PsiProbeInit_d, current_alphaInd_d, &dispatcher,
					                                current_psi_ds, current_psi_intensity_ds, current_integratedOutput_ds,
					                                GPU_num, current_qya_d, current_qxa_d, current_output_ph, current_cufft_plan,
					                                current_prop_d, &current_stream, &psi_size, stream_count, &PRISM_PRINT_FREQUENCY_PROBES]()  {

				// set the GPU context
				cudaErrchk(cudaSetDevice(GPU_num)); // set current GPU


#ifndef NDEBUG
				{
//					 check memory usage on the GPU
					std::lock_guard<mutex> lock(PRISM::mem_lock);
					size_t free_mem, total_mem;
					cudaErrchk(cudaMemGetInfo(&free_mem, &total_mem));
					pars.max_mem = std::max(total_mem - free_mem, pars.max_mem);
//					cout << "max_mem = " << pars.max_mem << endl;
				}
#endif // NDEBUG

				size_t Nstart, Nstop;
				Nstart=Nstop=0;
//				while (getWorkID(pars, Nstart, Nstop)){ // synchronously get work assignment
				while (dispatcher.getWork(Nstart, Nstop, pars.meta.batch_size_GPU)){ // synchronously get work assignment
					while (Nstart < Nstop){
						if (Nstart % PRISM_PRINT_FREQUENCY_PROBES < pars.meta.batch_size_GPU | Nstart == 100){
							cout << "Computing Probe Position #" << Nstart << "/" << pars.xp.size() * pars.yp.size() << endl;
						}
//						ay = Nstart / pars.xp.size();
//						ax = Nstart % pars.xp.size();
//
//						getMultisliceProbe_GPU_streaming(pars, current_trans_ds, trans_ph, current_PsiProbeInit_d, current_psi_ds,
//						                                 current_output_ph, current_psi_intensity_ds,
//						                                 current_integratedOutput_ds, current_qya_d, current_qxa_d,
//						                                 current_prop_d, ay, ax, pars.psiProbeInit.get_dimj(), pars.psiProbeInit.get_dimi(),
//						                                 current_alphaInd_d, current_cufft_plan, current_stream);
						getMultisliceProbe_GPU_streaming_batch(pars, current_trans_ds, trans_ph, current_PsiProbeInit_d, current_psi_ds,
						                                       current_output_ph, current_psi_intensity_ds,
						                                       current_integratedOutput_ds, current_qya_d, current_qxa_d,
						                                       current_prop_d, Nstart, Nstop, pars.psiProbeInit.get_dimj(), pars.psiProbeInit.get_dimi(),
						                                       current_alphaInd_d, current_cufft_plan, current_stream);
#ifdef PRISM_BUILDING_GUI
						pars.progressbar->signalOutputUpdate(Nstart, pars.xp.size() * pars.yp.size());
#endif
//						++Nstart;
						Nstart = Nstop;
					}
				}
				cout << "GPU worker on stream #" << stream_count << " of GPU #" << GPU_num << "finished\n";
			}));

			++stream_count;
		}


		// now launch CPU work

		if (pars.meta.also_do_CPU_work){
			PRISM_FFTW_INIT_THREADS();
			PRISM_FFTW_PLAN_WITH_NTHREADS(pars.meta.NUM_THREADS);vector<thread> workers_CPU;
			workers_CPU.reserve(pars.meta.NUM_THREADS); // prevents multiple reallocations
			for (auto t = 0; t < pars.meta.NUM_THREADS; ++t) {
				cout << "Launching CPU worker #" << t << endl;
				// push_back is better whenever constructing a new object
				workers_CPU.push_back(thread([&pars, &dispatcher, t, &PRISM_PRINT_FREQUENCY_PROBES]() {
				size_t Nstart, Nstop, early_CPU_stop;
				Nstart=Nstop=0;
				// stop the CPU workers earlier than the GPU ones to prevent slower workers taking the last jobs and having to
				// wait longer for everything to complete
                if (pars.meta.NUM_GPUS > 0){
                        // if there are no GPUs, make sure to do all work on CPU
                            early_CPU_stop = (size_t)std::max((PRISM_FLOAT_PRECISION)0.0, pars.xp.size() * pars.yp.size() - pars.meta.gpu_cpu_ratio);
                } else {
                            early_CPU_stop = pars.xp.size() * pars.yp.size();
                }
				if (dispatcher.getWork(Nstart, Nstop, pars.meta.batch_size_CPU, early_CPU_stop)) { // synchronously get work assignment
					Array1D<std::complex<PRISM_FLOAT_PRECISION> > psi_stack = zeros_ND<1, complex<PRISM_FLOAT_PRECISION> >({{pars.psiProbeInit.size() * pars.meta.batch_size_CPU}});

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

					gatekeeper.unlock();
					do {
						//	cout << "Nstop = " << Nstop << endl;
						while (Nstart < Nstop) {
							if (Nstart % PRISM_PRINT_FREQUENCY_PROBES  < pars.meta.batch_size_CPU | Nstart == 100){
								cout << "Computing Probe Position #" << Nstart << "/" << pars.xp.size() * pars.yp.size() << endl;
							}
//							ay = Nstart / pars.xp.size();
//							ax = Nstart % pars.xp.size();
//                            if (ay==7){
//                                cout << "ax = " << ax << endl;
//								cout << "ay = " << ay << endl;
//                            }
//							getMultisliceProbe_CPU(pars, ay, ax, plan_forward, plan_inverse, psi);
							getMultisliceProbe_CPU_batch(pars, Nstart, Nstop, plan_forward, plan_inverse, psi_stack);
#ifdef PRISM_BUILDING_GUI
							pars.progressbar->signalOutputUpdate(Nstart, pars.xp.size() * pars.yp.size());
#endif
							//++Nstart;
							Nstart=Nstop;
						}
						if (Nstop >= early_CPU_stop) break;
					} while(dispatcher.getWork(Nstart, Nstop, pars.meta.batch_size_CPU, early_CPU_stop));
					gatekeeper.lock();
					PRISM_FFTW_DESTROY_PLAN(plan_forward);
					PRISM_FFTW_DESTROY_PLAN(plan_inverse);
					gatekeeper.unlock();
				}
				cout << "CPU worker #" << t << " finished\n";

			}));

			}
			cout << "Waiting on GPU threads..." << endl;
			for (auto& t:workers_CPU)t.join();
			PRISM_FFTW_CLEANUP_THREADS();
		}
		// synchronize threads
		cout << "Waiting on GPU threads..." << endl;
		for (auto& t:workers_GPU)t.join();



		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g){
			cudaSetDevice(g);
			cudaDeviceSynchronize();
		}

		// copy the results of the GPU, which are in pinned memory, back to the actual stack. The CPU work populates the
		// beginning, so make sure to copy from the offset of where the GPU started. Launch this copy on a background thread
		// while we cleanup the GPU
//		const size_t GPU_start_offset = (size_t)CPU_stop*pars.output.get_dimk()*pars.output.get_dimj()*pars.output.get_dimi();
//		std::thread copy_t([&GPU_start_offset, &pars, &stack_ph](){
//			memcpy(&pars.output[GPU_start_offset],
//			       &stack_ph[GPU_start_offset],
//			       (pars.output.size()-GPU_start_offset) * sizeof(PRISM_FLOAT_PRECISION));
//		});

		// synchronize GPUs and cleanup data
		for (auto j = 0; j < pars.meta.NUM_GPUS; ++j){
			cudaErrchk(cudaSetDevice(j));
//			cudaErrchk(cudaDeviceSynchronize());
			cudaErrchk(cudaFree(PsiProbeInit_d[j]));
			cudaErrchk(cudaFree(trans_ds[j]));
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
		delete[] streams;
		delete[] cufft_plan;
		delete[] PsiProbeInit_d;
		delete[] trans_ds;
		delete[] prop_d;
		delete[] qxa_d;
		delete[] qya_d;
		delete[] alphaInd_d;
		delete[] psi_ds;
		delete[] psi_intensity_ds;
		delete[] integratedOutput_ds;
		delete[] output_ph;
	}

}
