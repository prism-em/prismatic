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
		cudaErrchk(cudaMallocHost((void **) &xBeams_ph, pars.xyBeams.get_dimj() * sizeof(size_t)));
		cudaErrchk(cudaMallocHost((void **) &yBeams_ph, pars.xyBeams.get_dimj() * sizeof(size_t)));



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
//
//		// pointers to read-only GPU memory (one copy per GPU)
//		PRISM_CUDA_COMPLEX_FLOAT *trans_d[pars.meta.NUM_GPUS];
//		PRISM_CUDA_COMPLEX_FLOAT *prop_d[pars.meta.NUM_GPUS];
//		size_t *qxInd_d[pars.meta.NUM_GPUS];
//		size_t *qyInd_d[pars.meta.NUM_GPUS];
//		size_t *beamsIndex_d[pars.meta.NUM_GPUS];
//
//		// pointers to read/write GPU memory (one per stream)
//		PRISM_CUDA_COMPLEX_FLOAT *psi_ds[total_num_streams];
//		PRISM_CUDA_COMPLEX_FLOAT *psi_small_ds[total_num_streams];
//
//		// allocate memory on each GPU
//		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
//			cudaErrchk(cudaSetDevice(g));
//			cudaErrchk(cudaMalloc((void **) &trans_d[g], trans.size() * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
//			cudaErrchk(cudaMalloc((void **) &prop_d[g], pars.prop.size() * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
//			cudaErrchk(cudaMalloc((void **) &qxInd_d[g], pars.qxInd.size() * sizeof(size_t)));
//			cudaErrchk(cudaMalloc((void **) &qyInd_d[g], pars.qyInd.size() * sizeof(size_t)));
//			cudaErrchk(cudaMalloc((void **) &beamsIndex_d[g], pars.beamsIndex.size() * sizeof(size_t)));
//		}
//
//		// allocate memory per stream and 0 it
//		for (auto s = 0; s < total_num_streams; ++s) {
//			cudaErrchk(cudaSetDevice(s % pars.meta.NUM_GPUS));
//			cudaErrchk(cudaMalloc((void **) &psi_ds[s],
//			                      pars.imageSize[0] * pars.imageSize[1] * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
//			cudaErrchk(cudaMalloc((void **) &psi_small_ds[s],
//			                      pars.qxInd.size() * pars.qyInd.size() * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
//			cudaErrchk(
//					cudaMemset(psi_ds[s], 0, pars.imageSize[0] * pars.imageSize[1] * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
//			cudaErrchk(cudaMemset(psi_small_ds[s], 0,
//			                      pars.qxInd.size() * pars.qyInd.size() * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
//		}


		// launch threads to compute results for batches of xp, yp
		// I do this by dividing the xp points among threads, and each computes
		// all of the relevant yp for each of its xp. This seems an okay strategy
		// as long as the number of xp and yp are similar. If that is not the case
		// this may need to be adapted
		vector<thread> workers;
		workers.reserve(pars.meta.NUM_THREADS); // prevents multiple reallocations
		setWorkStartStop(0, pars.xp.size() * pars.yp.size());
		for (auto t = 0; t < pars.meta.NUM_THREADS; ++t) {
			cout << "Launching thread #" << t << " to result\n";
			// emplace_back is better whenever constructing a new object
			workers.emplace_back(thread([&pars, &xTiltShift, &yTiltShift,
					                            &alphaInd, &PsiProbeInit]() {
				size_t Nstart, Nstop, ay, ax;
				while (getWorkID(pars, Nstart, Nstop)) { // synchronously get work assignment
					while (Nstart != Nstop) {
						ay = Nstart / pars.xp.size();
						ax = Nstart % pars.xp.size();
						buildSignal_GPU(pars, ay, ax, yTiltShift, xTiltShift, alphaInd, PsiProbeInit);
						++Nstart;
					}
				}
			}));
		}
		// synchronize
		cout << "Waiting for threads...\n";
		for (auto &t:workers)t.join();





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



	}

	void buildSignal_GPU(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                     const size_t &ay,
	                     const size_t &ax,
	                     const PRISM_FLOAT_PRECISION &yTiltShift,
	                     const PRISM_FLOAT_PRECISION &xTiltShift,
	                     const Array2D<PRISM_FLOAT_PRECISION> &alphaInd,
	                     const Array2D<std::complex<PRISM_FLOAT_PRECISION> > &PsiProbeInit) {
//		cout << "DUMMY CODE" << endl;
	}
}