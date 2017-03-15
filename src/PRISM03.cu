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
namespace PRISM {
	using namespace std;

	void buildPRISMOutput_GPU(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                          const PRISM_FLOAT_PRECISION xTiltShift,
	                          const PRISM_FLOAT_PRECISION yTiltShift,
	                          const Array2D<PRISM_FLOAT_PRECISION> &alphaInd,
	                          const Array2D<std::complex<PRISM_FLOAT_PRECISION> > &PsiProbeInit) {

		// create CUDA streams
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;
		cudaStream_t streams[total_num_streams];
		cufftHandle cufft_plan[total_num_streams];

		cout << "pars.imageSizeReduce[0] = " << pars.imageSizeReduce[0] << endl;
		cout << "pars.imageSizeReduce[1] = " << pars.imageSizeReduce[1] << endl;
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
		PRISM_FLOAT_PRECISION               *xVec_ph;
		PRISM_FLOAT_PRECISION               *yVec_ph;
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
		cudaErrchk(cudaMallocHost((void **) &xVec_ph, pars.xVec.size() * sizeof(PRISM_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **) &yVec_ph, pars.xVec.size() * sizeof(PRISM_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **) &qxaReduce_ph, pars.qxaReduce.size() * sizeof(PRISM_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **) &qyaReduce_ph, pars.qyaReduce.size() * sizeof(PRISM_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **) &alphaInd_ph,  alphaInd.size()       * sizeof(PRISM_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **) &xBeams_ph, pars.xyBeams.get_dimj()  * sizeof(size_t)));
		cudaErrchk(cudaMallocHost((void **) &yBeams_ph, pars.xyBeams.get_dimj()  * sizeof(size_t)));



		// copy host memory to pinned
		for (auto s = 0; s < total_num_streams; ++s) {
			memset(output_ph[s], 0, pars.stack.get_dimj() * pars.stack.get_dimi() *
			                        sizeof(PRISM_FLOAT_PRECISION));
		}

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
		memcpy(xVec_ph, &(*pars.xVec.begin()), pars.xVec.size() * sizeof(PRISM_FLOAT_PRECISION));
		memcpy(yVec_ph, &(*pars.yVec.begin()), pars.yVec.size() * sizeof(PRISM_FLOAT_PRECISION));
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


		for (auto i = 0; i < 10; ++i){
			cout << "PsiProbeInit[" << i << "] = " << PsiProbeInit[i] << endl;
			cout << "PsiProbeInit_ph[" << i << "] = " << PsiProbeInit_ph[i] << endl;
		}
		for (auto i = 0; i < 10; ++i){
			cout << "permuted_Scompact_ph[" << i << "] = " << permuted_Scompact_ph[i] << endl;
			cout << "pars.Scompact.at("  << i << ",0,0) = " << pars.Scompact.at(i,0,0) << endl;
		}

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
		size_t                   *y_ds[total_num_streams];
		size_t                   *x_ds[total_num_streams];
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
			                      pars.xyBeams.get_dimj() * sizeof(size_t)));
			cudaErrchk(cudaMalloc((void **) &x_ds[s],
			                      pars.xyBeams.get_dimj() * sizeof(size_t)));

			cudaErrchk(cudaMemset(psi_ds[s], 0,
			                      pars.imageSizeReduce[0] * pars.imageSizeReduce[1] * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMemset(phaseCoeffs_ds[s], 0,
			                      pars.numberBeams * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMemset(psi_intensity_ds[s], 0,
			                      pars.imageSizeReduce[0] * pars.imageSizeReduce[1] * sizeof(PRISM_FLOAT_PRECISION)));
			cudaErrchk(cudaMemset(y_ds[s], 0,
			                      pars.xyBeams.get_dimj() * sizeof(size_t)));
			cudaErrchk(cudaMemset(x_ds[s], 0,
			                      pars.xyBeams.get_dimj() * sizeof(size_t)));
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
		for (auto t = 0; t < pars.meta.NUM_THREADS; ++t) {
			cout << "Launching thread #" << t << " to result\n";

			int GPU_num = stream_count % pars.meta.NUM_GPUS; // determine which GPU handles this job
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
			size_t *current_y_ds                                  = y_ds[stream_count];
			size_t *current_x_ds                                  = x_ds[stream_count];
			cufftHandle &current_cufft_plan                       = cufft_plan[stream_count];

			// get pointer to output pinned memory
			PRISM_FLOAT_PRECISION *current_output_ph              = output_ph[total_num_streams];


			// emplace_back is better whenever constructing a new object
			workers_GPU.emplace_back(thread([&pars, GPU_num, stream_count, &yTiltShift, &xTiltShift, current_permuted_Scompact_d,
					                                current_alphaInd_d, current_PsiProbeInit_d, current_qxaReduce_d, current_qyaReduce_d,
					                                current_yBeams_d, current_xBeams_d, current_psi_ds, current_phaseCoeffs_ds,
					                                current_psi_intensity_ds, current_y_ds, current_x_ds,
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
						                current_x_ds, current_output_ph, current_cufft_plan, current_stream );



//						buildSignal_CPU(pars, ay, ax, yTiltShift, xTiltShift, alphaInd, PsiProbeInit);
						++Nstart;
					}
				}
			}));
			++stream_count;
		}
		// synchronize
		cout << "Waiting for threads...\n";
		for (auto &t:workers_GPU)t.join();

		// free pinned memory
		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaFreeHost(output_ph[s]));
		}
		cudaErrchk(cudaFreeHost(permuted_Scompact_ph));
		cudaErrchk(cudaFreeHost(PsiProbeInit_ph));
		cudaErrchk(cudaFreeHost(xVec_ph));
		cudaErrchk(cudaFreeHost(yVec_ph));
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
		}

		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaSetDevice(s % pars.meta.NUM_GPUS));
			cudaErrchk(cudaFree(psi_ds[s]));
			cudaErrchk(cudaFree(phaseCoeffs_ds[s]));
			cudaErrchk(cudaFree(psi_intensity_ds[s]));
			cudaErrchk(cudaFree(alphaInd_d[s]));
			cudaErrchk(cudaFree(y_ds[s]));
			cudaErrchk(cudaFree(x_ds[s]));
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
	                     size_t *y_ds,
	                     size_t *x_ds,
	                     PRISM_FLOAT_PRECISION *output_ph,
	                     const cufftHandle &cufft_plan,
	                     const cudaStream_t& stream){

		
		cout <<"test\n";
	}
}