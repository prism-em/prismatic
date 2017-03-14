#include "PRISM02.cuh"
#include "configure.h"
#include <mutex>
#include <thread>
#include <complex>
#include <vector>
#include "getWorkID.h"
#include "fftw3.h"
#include "defines.h"
#include "cufft.h"
namespace PRISM {
	using namespace std;
	void propagatePlaneWave_GPU(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                            PRISM_CUDA_COMPLEX_FLOAT* trans_d,
	                            PRISM_FLOAT_PRECISION* Scompact_slice_ph,
	                            PRISM_CUDA_COMPLEX_FLOAT* psi_d,
	                            const PRISM_FLOAT_PRECISION* qyInd_d,
	                            const PRISM_FLOAT_PRECISION* qxInd_d,
	                            const PRISM_CUDA_COMPLEX_FLOAT* prop_d,
	                            const size_t dimj,
	                            const size_t dimi,
	                            const size_t& beamNumber,
	                            const cufftHandle& plan,
	                            const cufftHandle& plan_small,
	                            cudaStream_t& stream){

//
//		psi[pars.beamsIndex[a0]] = 1;
//		const PRISM_FLOAT_PRECISION N = (PRISM_FLOAT_PRECISION)psi.size();
//
//
//		PRISM_FFTW_EXECUTE(plan_inverse);
//		for (auto &i : psi)i /= N; // fftw scales by N, need to correct
//		const complex<PRISM_FLOAT_PRECISION>* trans_t = &trans[0];
//		for (auto a2 = 0; a2 < pars.numPlanes; ++a2){
//
//			for (auto& p:psi)p*=(*trans_t++); // transmit
//			PRISM_FFTW_EXECUTE(plan_forward); // FFT
//			for (auto i = psi.begin(), j = pars.prop.begin(); i != psi.end();++i,++j)*i*=(*j); // propagate
//			PRISM_FFTW_EXECUTE(plan_inverse); // IFFT
//			for (auto &i : psi)i /= N; // fftw scales by N, need to correct
//		}
//		PRISM_FFTW_EXECUTE(plan_forward);
//
//		Array2D< complex<PRISM_FLOAT_PRECISION> > psi_small = zeros_ND<2, complex<PRISM_FLOAT_PRECISION> >({{pars.qyInd.size(), pars.qxInd.size()}});
//
//
//		unique_lock<mutex> gatekeeper(fftw_plan_lock);
//		PRISM_FFTW_PLAN plan_final = PRISM_FFTW_PLAN_DFT_2D(psi_small.get_dimj(), psi_small.get_dimi(),
//		                                                    reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi_small[0]),
//		                                                    reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi_small[0]),
//		                                                    FFTW_BACKWARD, FFTW_ESTIMATE);
//		gatekeeper.unlock();
//		for (auto y = 0; y < pars.qyInd.size(); ++y){
//			for (auto x = 0; x < pars.qxInd.size(); ++x){
//				psi_small.at(y,x) = psi.at(pars.qyInd[y], pars.qxInd[x]);
//			}
//		}
//		PRISM_FFTW_EXECUTE(plan_final);
//		gatekeeper.lock();
//		PRISM_FFTW_DESTROY_PLAN(plan_final);
//		gatekeeper.unlock();
//
//
//		complex<PRISM_FLOAT_PRECISION>* S_t = &pars.Scompact[a0 * pars.Scompact.get_dimj() * pars.Scompact.get_dimi()];
//		const PRISM_FLOAT_PRECISION N_small = (PRISM_FLOAT_PRECISION)psi_small.size();
//		for (auto& i:psi_small) {
//			*S_t++ = i/N_small;
//		}
//	};

	}

	void fill_Scompact_GPU(Parameters <PRISM_FLOAT_PRECISION> &pars){

		//initialize data

		// psi per stream
		// trans per gpu
		// prop per gpu
		// Scompact pinned memory chunk per stream
		// 2 cufft plans, one for big and one for small psi array
		// qxind/qyind per gpu

//		mutex fftw_plan_lock;
		cout << "entering GPU" << endl;
		const PRISM_FLOAT_PRECISION pi = acos(-1);
		const std::complex<PRISM_FLOAT_PRECISION> i(0, 1);
		pars.Scompact = zeros_ND<3, complex<PRISM_FLOAT_PRECISION> > ({{pars.numberBeams,pars.imageSize[0]/2, pars.imageSize[1]/2}});
		Array3D<complex<PRISM_FLOAT_PRECISION> > trans = zeros_ND<3, complex<PRISM_FLOAT_PRECISION> >(
				{{pars.pot.get_dimk(), pars.pot.get_dimj(), pars.pot.get_dimi()}});
		{
			auto p = pars.pot.begin();
			for (auto &j:trans)j = exp(i * pars.sigma * (*p++));
		}




		// create CUDA streams
		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;
		cudaStream_t streams[total_num_streams];
		cufftHandle cufft_plan[total_num_streams];
		cufftHandle cufft_plan_small[total_num_streams];
		cout <<"total_num_streams = " << total_num_streams<< endl;
		for (auto j = 0; j < total_num_streams; ++j){
			cudaSetDevice(j % pars.meta.NUM_GPUS);
			cudaErrchk(cudaStreamCreate(&streams[j]));
			cufftErrchk(cufftPlan2d(&cufft_plan[j], pars.imageSize[1], pars.imageSize[0], CUFFT_C2C));
			cufftErrchk(cufftPlan2d(&cufft_plan_small[j], pars.imageSize[1], pars.imageSize[0], CUFFT_C2C));
			cufftErrchk(cufftSetStream(cufft_plan[j], streams[j]));
			cufftErrchk(cufftSetStream(cufft_plan_small[j], streams[j]));
		}

		// pointers to pinned host memory for async transfers
		std::complex<PRISM_FLOAT_PRECISION>  *trans_ph;
		std::complex<PRISM_FLOAT_PRECISION>  *prop_ph;
		std::complex<PRISM_FLOAT_PRECISION>  *Scompact_slice_ph[total_num_streams];
		PRISM_FLOAT_PRECISION                *qxInd_ph;
		PRISM_FLOAT_PRECISION                *qyInd_ph;
		size_t					             *beamsIndex_ph;

		// allocate pinned memory
		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaMallocHost((void **) &Scompact_slice_ph[s],
			                          pars.Scompact.get_dimj() * pars.Scompact.get_dimi() *
			                          sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
		}
		cudaErrchk(cudaMallocHost((void **)&trans_ph, trans.size()*sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
		cudaErrchk(cudaMallocHost((void **)&prop_ph,  pars.prop.size()*sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
		cudaErrchk(cudaMallocHost((void **)&qxInd_ph, pars.qxInd.size()*sizeof(PRISM_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **)&qyInd_ph, pars.qyInd.size()*sizeof(PRISM_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocHost((void **)&beamsIndex_ph, pars.beamsIndex.size()*sizeof(size_t)));

		// copy host memory to pinned
		for (auto s = 0; s < total_num_streams; ++s) {
			memset(&Scompact_slice_ph[s], 0, pars.Scompact.get_dimj() * pars.Scompact.get_dimi() *
			                                 sizeof(std::complex<PRISM_FLOAT_PRECISION>));
		}
		memcpy(trans_ph,   &trans[0],       trans.size()      * sizeof(std::complex<PRISM_FLOAT_PRECISION>));
		memcpy(prop_ph,    &pars.prop[0],   pars.prop.size()  * sizeof(std::complex<PRISM_FLOAT_PRECISION>));
		memcpy(qxInd_ph,   &pars.qxInd[0],  pars.qxInd.size() * sizeof(PRISM_FLOAT_PRECISION));
		memcpy(qyInd_ph,   &pars.qyInd[0],  pars.qyInd.size() * sizeof(PRISM_FLOAT_PRECISION));
		memcpy(beamsIndex_ph,   &pars.beamsIndex[0],  pars.beamsIndex.size() * sizeof(size_t));

		// pointers to read-only GPU memory (one copy per GPU)
		PRISM_CUDA_COMPLEX_FLOAT *trans_d[pars.meta.NUM_GPUS];
		PRISM_CUDA_COMPLEX_FLOAT *prop_d[pars.meta.NUM_GPUS];
		PRISM_FLOAT_PRECISION    *qxInd_d[pars.meta.NUM_GPUS];
		PRISM_FLOAT_PRECISION    *qyInd_d[pars.meta.NUM_GPUS];
		size_t                   *beamsIndex[pars.meta.NUM_GPUS];

		// pointers to read/write GPU memory (one per stream)
		PRISM_CUDA_COMPLEX_FLOAT *psi_ds[total_num_streams];
		PRISM_CUDA_COMPLEX_FLOAT *psi_small_ds[total_num_streams];

		// allocate memory on each GPU
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaMalloc((void **) &trans_d[g],    trans.size()      * sizeof(trans[0])));
			cudaErrchk(cudaMalloc((void **) &prop_d[g],     pars.prop.size()  * sizeof(pars.prop[0])));
			cudaErrchk(cudaMalloc((void **) &qxInd_d[g],    pars.qxInd.size() * sizeof(pars.qxInd[0])));
			cudaErrchk(cudaMalloc((void **) &qyInd_d[g],    pars.qyInd.size() * sizeof(pars.qyInd[0])));
			cudaErrchk(cudaMalloc((void **) &beamsIndex[g], pars.beamsIndex.size() * sizeof(pars.beamsIndex[0])));
		}

		// allocate memory per stream and 0 it
		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaSetDevice(s % pars.meta.NUM_GPUS));
			cudaErrchk(cudaMalloc((void **) &psi_ds[s],       pars.imageSize[0] * pars.imageSize[1] * sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
			cudaErrchk(cudaMalloc((void **) &psi_small_ds[s], pars.qxInd.size() * pars.qyInd.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
			cudaErrchk(cudaMemset(psi_ds[s], 0, pars.imageSize[0] * pars.imageSize[1] * sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
			cudaErrchk(cudaMemset(psi_small_ds[s], 0, pars.qxInd.size() * pars.qyInd.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
		}


		// Copy memory to each GPU asynchronously from the pinned host memory spaces.
		// The streams are laid out so that consecutive streams represent different GPUs. If we
		// have more than one stream per GPU, then we want to interleave as much as possible
		int stream_id = 0;
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			stream_id = g;
			cudaErrchk(cudaSetDevice(g));
			cout << "stream_id = " << stream_id << endl;
			cudaErrchk(cudaMemcpyAsync(trans_d[g], &trans_ph[0],
			                           trans.size() * sizeof(trans[0]), cudaMemcpyHostToDevice, streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cout << "stream_id = " << stream_id << endl;
			cudaErrchk(cudaMemcpyAsync(prop_d[g], &prop_ph[0],
		                           pars.prop.size() * sizeof(pars.prop[0]), cudaMemcpyHostToDevice, streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cout << "stream_id = " << stream_id << endl;
			cudaErrchk(cudaMemcpyAsync(qxInd_d[g], &qxInd_ph[0],
			                           pars.qxInd.size() * sizeof(pars.qxInd[0]), cudaMemcpyHostToDevice, streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cout << "stream_id = " << stream_id << endl;
			cudaErrchk(cudaMemcpyAsync(qyInd_d[g], &qyInd_ph[0],
			                           pars.qyInd.size() * sizeof(pars.qyInd[0]), cudaMemcpyHostToDevice, streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cout << "stream_id = " << stream_id << endl;
			cudaErrchk(cudaMemcpyAsync(beamsIndex[g], &beamsIndex_ph[0],
			                           pars.qyInd.size() * sizeof(pars.qyInd[0]), cudaMemcpyHostToDevice, streams[stream_id]));
		}









		// launch GPU work
		vector<thread> workers_GPU;
		workers_GPU.reserve(total_num_streams); // prevents multiple reallocations

		int stream_count = 0;
		setWorkStartStop(0, pars.numberBeams);
		for (auto t = 0; t < total_num_streams; ++t){

			int GPU_num = stream_count % pars.meta.NUM_GPUS; // determine which GPU handles this job
			cudaStream_t& current_stream = streams[stream_count];
			cout << "Launching GPU worker on stream #" << stream_count << " of GPU #" << GPU_num << '\n';

			// get pointers to the pre-copied arrays, making sure to get those on the current GPU
			PRISM_CUDA_COMPLEX_FLOAT *current_trans_d = trans_d[GPU_num];
			PRISM_CUDA_COMPLEX_FLOAT *current_prop_d  = prop_d[GPU_num];
			PRISM_FLOAT_PRECISION *current_qxa_d      = qxInd_d[GPU_num];
			PRISM_FLOAT_PRECISION *current_qya_d      = qyInd_d[GPU_num];

			// get pointers to per-stream arrays
			PRISM_CUDA_COMPLEX_FLOAT *current_psi_ds        = psi_ds[stream_count];
			PRISM_CUDA_COMPLEX_FLOAT *current_psi_small_ds  = psi_small_ds[stream_count];
			cufftHandle& current_cufft_plan                 = cufft_plan[stream_count];
			cufftHandle& current_cufft_plan_small           = cufft_plan_small[stream_count];

			workers_GPU.emplace_back(thread([&pars, current_trans_d, current_prop_d, current_qxa_d, current_qya_d,
					                                current_psi_ds, current_psi_small_ds, current_cufft_plan, current_cufft_plan_small](){
				size_t currentBeam, stop;
				while (getWorkID(pars, currentBeam, stop)){
					while(currentBeam != stop){
//						propagatePlaneWave_GPU(pars, trans, currentBeam, psi);
						++currentBeam;
					}
				}

			}));
		}

		// launch CPU work
		vector<thread> workers_CPU;
		workers_CPU.reserve(pars.meta.NUM_THREADS); // prevents multiple reallocations



		for (auto &t:workers_GPU)t.join();
		for (auto &t:workers_CPU)t.join();

		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaSetDevice(s % pars.meta.NUM_GPUS));
//			cudaErrchk(cudaFree(psi_ds[s]));
//			cudaErrchk(cudaFree(psi_intensity_ds[s]));
//			cudaErrchk(cudaFree(integratedOutput_ds[s]));
			cufftErrchk(cufftDestroy(cufft_plan[s]));
			cufftErrchk(cufftDestroy(cufft_plan_small[s]));
		}


		// free pinned memory
//		cudaErrchk(cudaFreeHost(PsiProbeInit_ph));
//		cudaErrchk(cudaFreeHost(trans_ph));
//		cudaErrchk(cudaFreeHost(prop_ph));
//		cudaErrchk(cudaFreeHost(qxInd_ph));
//		cudaErrchk(cudaFreeHost(qyInd_ph));
//		cudaErrchk(cudaFreeHost(alphaInd_ph));
//		cudaErrchk(cudaFreeHost(stack_ph));


		// destroy CUDA streams
		for (auto j = 0; j < total_num_streams; ++j){
			cudaSetDevice(j % pars.meta.NUM_GPUS);
			cudaErrchk(cudaStreamDestroy(streams[j]));
		}
	}
}