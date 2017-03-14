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
#define PI 3.14159265359
#define BLOCK_SIZE1D 1024

namespace PRISM {
	using namespace std;

	// define some constants
	__device__ __constant__ PRISM_FLOAT_PRECISION pi       = PI;
	__device__ __constant__ PRISM_CUDA_COMPLEX_FLOAT i     = {0, 1};
	__device__ __constant__ PRISM_CUDA_COMPLEX_FLOAT pi_cx = {PI, 0};
	__device__ __constant__ PRISM_CUDA_COMPLEX_FLOAT minus_2pii = {0, -2*PI};

	// computes exp(real(a) + i * imag(a))
	__device__ __forceinline__ cuDoubleComplex exp_cx(const cuDoubleComplex a){
		double e = exp(a.x);
		double s,c;
		sincos(a.y, &s, &c);
		return make_cuDoubleComplex(e*c, e*s);
	}
	__device__ __forceinline__ cuFloatComplex exp_cx(const cuFloatComplex a){
		float e = expf(a.x);
		float s,c;
		sincosf(a.y, &s, &c);
		return make_cuFloatComplex(e*c, e*s);
	}

	// creates initial probe using existing GPU memory rather than streaming each probe
	__global__ void initializePsi(cuFloatComplex *psi_d, const size_t N, const size_t beamLoc){
		int idx = threadIdx.x + blockDim.x*blockIdx.x;
		if (idx < N) {
			psi_d[idx] = (idx == beamLoc) ? make_cuFloatComplex(1,0):make_cuFloatComplex(0,0);
		}
	}

	__global__ void initializePsi(cuDoubleComplex *psi_d, const size_t N, const size_t beamLoc){
		int idx = threadIdx.x + blockDim.x*blockIdx.x;
		if (idx < N) {
			psi_d[idx] = (idx == beamLoc) ? make_cuDoubleComplex(1,0):make_cuDoubleComplex(0,0);
		}
	}

	// multiply two complex arrays
	__global__ void multiply_inplace(PRISM_CUDA_COMPLEX_FLOAT* arr,
	                                 const PRISM_CUDA_COMPLEX_FLOAT* other,
	                                 const size_t N){
		int idx = threadIdx.x + blockDim.x*blockIdx.x;
		if (idx < N) {
			PRISM_CUDA_COMPLEX_FLOAT a = arr[idx];
			PRISM_CUDA_COMPLEX_FLOAT o = other[idx];
			arr[idx].x = a.x * o.x - a.y * o.y;
			arr[idx].y = a.x * o.y + a.y * o.x;
		}
	}

	// divide two complex arrays
	__global__ void divide_inplace(PRISM_CUDA_COMPLEX_FLOAT* arr,
	                               const PRISM_FLOAT_PRECISION val,
	                               const size_t N){
		int idx = threadIdx.x + blockDim.x*blockIdx.x;
		if (idx < N) {
			arr[idx].x /= val;
			arr[idx].y /= val;
		}
	}

	__global__ void array_subset(const PRISM_CUDA_COMPLEX_FLOAT* psi_d,
	                             PRISM_CUDA_COMPLEX_FLOAT* psi_small_d,
	                             const PRISM_FLOAT_PRECISION* qyInd_d,
	                             const PRISM_FLOAT_PRECISION* qxInd_d,
	                             const size_t dimj,
	                             const size_t dimi,
	                             const size_t N){
		int idx = threadIdx.x + blockDim.x*blockIdx.x;
		if (idx < N) {
//			arr[idx].x /= val;
//			arr[idx].y /= val;
		}
	}

	void propagatePlaneWave_GPU(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                            PRISM_CUDA_COMPLEX_FLOAT* trans_d,
	                            PRISM_CUDA_COMPLEX_FLOAT* psi_d,
	                            PRISM_CUDA_COMPLEX_FLOAT* psi_small_d,
	                            complex<PRISM_FLOAT_PRECISION>* Scompact_slice_ph,
	                            const PRISM_FLOAT_PRECISION* qyInd_d,
	                            const PRISM_FLOAT_PRECISION* qxInd_d,
	                            const PRISM_CUDA_COMPLEX_FLOAT* prop_d,
	                            const size_t* beamsIndex,
	                            const size_t& beamNumber,
	                            const cufftHandle& plan,
	                            const cufftHandle& plan_small,
	                            cudaStream_t& stream){

		const size_t psi_size = pars.imageSize[0] * pars.imageSize[1];
		initializePsi<<< (psi_size - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D >>>(psi_d, psi_size, pars.beamsIndex[beamNumber]);
//		complex<float> ans;
//		cudaMemcpy(&ans, psi_d, sizeof(ans), cudaMemcpyDeviceToHost);
//		cout << "ans[0] = " << ans << endl;
//		cudaMemcpy(&ans, (psi_d + 1), sizeof(ans), cudaMemcpyDeviceToHost);
//		cout << "ans[1] = " << ans << endl;
//		cudaMemcpy(&ans, (psi_d + 15), sizeof(ans), cudaMemcpyDeviceToHost);
//		cout << "ans[15] = " << ans << endl;

		for (auto planeNum = 0; planeNum < pars.numPlanes; ++planeNum) {
			cufftErrchk(cufftExecC2C(plan, &psi_d[0], &psi_d[0], CUFFT_INVERSE));
			multiply_inplace<<<(psi_size-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_d, &trans_d[planeNum*psi_size], psi_size);
			cufftErrchk(cufftExecC2C(plan, &psi_d[0], &psi_d[0], CUFFT_FORWARD));
			multiply_inplace<<<(psi_size-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_d, prop_d, psi_size);
			divide_inplace<<<(psi_size-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_d, psi_size, psi_size);
		}

//		formatOutput_GPU(pars, psi_intensity_ds, alphaInd_d, stack_ph, integratedOutput_ds, ay, ax, dimj, dimi,stream);

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
		size_t                               *qxInd_ph;
		size_t                               *qyInd_ph;
		size_t					             *beamsIndex_ph;

		// allocate pinned memory
		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaMallocHost((void **) &Scompact_slice_ph[s],
			                          pars.Scompact.get_dimj() * pars.Scompact.get_dimi() *
			                          sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
		}
		cudaErrchk(cudaMallocHost((void **)&trans_ph, trans.size()*sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
		cudaErrchk(cudaMallocHost((void **)&prop_ph,  pars.prop.size()*sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
		cudaErrchk(cudaMallocHost((void **)&qxInd_ph, pars.qxInd.size()*sizeof(size_t)));
		cudaErrchk(cudaMallocHost((void **)&qyInd_ph, pars.qyInd.size()*sizeof(size_t)));
		cudaErrchk(cudaMallocHost((void **)&beamsIndex_ph, pars.beamsIndex.size()*sizeof(size_t)));

		cout << "copy to pinned" << endl;
		// copy host memory to pinned
		for (auto s = 0; s < total_num_streams; ++s) {
			memset(Scompact_slice_ph[s], 0, pars.Scompact.get_dimj() * pars.Scompact.get_dimi() *
			                                 sizeof(std::complex<PRISM_FLOAT_PRECISION>));
		}
		memcpy(trans_ph,   &trans[0],       trans.size()      * sizeof(std::complex<PRISM_FLOAT_PRECISION>));
		memcpy(prop_ph,    &pars.prop[0],   pars.prop.size()  * sizeof(std::complex<PRISM_FLOAT_PRECISION>));
		memcpy(qxInd_ph,   &pars.qxInd[0],  pars.qxInd.size() * sizeof(size_t));
		memcpy(qyInd_ph,   &pars.qyInd[0],  pars.qyInd.size() * sizeof(size_t));
		memcpy(beamsIndex_ph,   &pars.beamsIndex[0],  pars.beamsIndex.size() * sizeof(size_t));
		cout << "pars.qxInd.size() = " << pars.qxInd.size() << endl;
		cout << "memcpy done" << endl;

		for (auto i = 0; i < 10; ++i){
			cout << "qxInd_ph[i] = " << qxInd_ph[i] << endl;
			cout << "pars.qxInd[[i] = " << pars.qxInd[i] << endl;
		}
		// pointers to read-only GPU memory (one copy per GPU)
		PRISM_CUDA_COMPLEX_FLOAT *trans_d[pars.meta.NUM_GPUS];
		PRISM_CUDA_COMPLEX_FLOAT *prop_d[pars.meta.NUM_GPUS];
		PRISM_FLOAT_PRECISION    *qxInd_d[pars.meta.NUM_GPUS];
		PRISM_FLOAT_PRECISION    *qyInd_d[pars.meta.NUM_GPUS];
		size_t                   *beamsIndex_d[pars.meta.NUM_GPUS];

		// pointers to read/write GPU memory (one per stream)
		PRISM_CUDA_COMPLEX_FLOAT *psi_ds[total_num_streams];
		PRISM_CUDA_COMPLEX_FLOAT *psi_small_ds[total_num_streams];

		cout << "allocate gpu memory" << endl;
		// allocate memory on each GPU
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaMalloc((void **) &trans_d[g],      trans.size()      * sizeof(trans[0])));
			cudaErrchk(cudaMalloc((void **) &prop_d[g],       pars.prop.size()  * sizeof(pars.prop[0])));
			cudaErrchk(cudaMalloc((void **) &qxInd_d[g],      pars.qxInd.size() * sizeof(pars.qxInd[0])));
			cudaErrchk(cudaMalloc((void **) &qyInd_d[g],      pars.qyInd.size() * sizeof(pars.qyInd[0])));
			cudaErrchk(cudaMalloc((void **) &beamsIndex_d[g], pars.beamsIndex.size() * sizeof(pars.beamsIndex[0])));
		}

		// allocate memory per stream and 0 it
		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaSetDevice(s % pars.meta.NUM_GPUS));
			cudaErrchk(cudaMalloc((void **) &psi_ds[s],       pars.imageSize[0] * pars.imageSize[1] * sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
			cudaErrchk(cudaMalloc((void **) &psi_small_ds[s], pars.qxInd.size() * pars.qyInd.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
			cudaErrchk(cudaMemset(psi_ds[s], 0, pars.imageSize[0] * pars.imageSize[1] * sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
			cudaErrchk(cudaMemset(psi_small_ds[s], 0, pars.qxInd.size() * pars.qyInd.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>)));
		}


		cout << "asynch copies" << endl;
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
			cout << "qxInd_ph[1] = " << qxInd_ph[1]<< endl;
			cout << "pars.qxInd[1] = " << pars.qxInd[1]<< endl;
			cudaErrchk(cudaMemcpyAsync(qxInd_d[g], &qxInd_ph[0],
			                           pars.qxInd.size() * sizeof(pars.qxInd[0]), cudaMemcpyHostToDevice, streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cout << "stream_id = " << stream_id << endl;
			cudaErrchk(cudaMemcpyAsync(qyInd_d[g], &qyInd_ph[0],
			                           pars.qyInd.size() * sizeof(pars.qyInd[0]), cudaMemcpyHostToDevice, streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cout << "stream_id = " << stream_id << endl;
			cudaErrchk(cudaMemcpyAsync(beamsIndex_d[g], &beamsIndex_ph[0],
			                           pars.beamsIndex.size() * sizeof(size_t), cudaMemcpyHostToDevice, streams[stream_id]));
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
			PRISM_FLOAT_PRECISION *current_qxInd_d    = qxInd_d[GPU_num];
			PRISM_FLOAT_PRECISION *current_qyInd_d    = qyInd_d[GPU_num];
			size_t                *current_beamsIndex = beamsIndex_d[GPU_num];
			// get pointers to per-stream arrays
			PRISM_CUDA_COMPLEX_FLOAT *current_psi_ds             = psi_ds[stream_count];
			PRISM_CUDA_COMPLEX_FLOAT *current_psi_small_ds       = psi_small_ds[stream_count];
			cufftHandle& current_cufft_plan                      = cufft_plan[stream_count];
			cufftHandle& current_cufft_plan_small                = cufft_plan_small[stream_count];
			complex<PRISM_FLOAT_PRECISION > *current_S_slice_ph  = Scompact_slice_ph[stream_count];

			workers_GPU.emplace_back(thread([&pars, current_trans_d, current_prop_d, current_qxInd_d, current_qyInd_d,
					                                current_psi_ds, current_psi_small_ds, &current_cufft_plan, &current_cufft_plan_small,
					                                current_S_slice_ph, current_beamsIndex, GPU_num, stream_count, &current_stream](){
				size_t currentBeam, stop;
				while (getWorkID(pars, currentBeam, stop)){
					while(currentBeam != stop){
						propagatePlaneWave_GPU(pars,
						                       current_trans_d,
						                       current_psi_ds,
						                       current_psi_small_ds,
						                       current_S_slice_ph,
						                       current_qyInd_d,
						                       current_qxInd_d,
						                       current_prop_d,
						                       current_beamsIndex,
						                       currentBeam,
						                       current_cufft_plan,
						                       current_cufft_plan_small,
						                       current_stream);
//
//						(Parameters<PRISM_FLOAT_PRECISION>& pars,
//								PRISM_CUDA_COMPLEX_FLOAT* trans_d,
//								PRISM_CUDA_COMPLEX_FLOAT* PsiProbeInit_d,
//								PRISM_CUDA_COMPLEX_FLOAT* psi_ds,
//								PRISM_FLOAT_PRECISION* stack_ph,
//								PRISM_FLOAT_PRECISION* psi_intensity_ds,
//								PRISM_FLOAT_PRECISION* integratedOutput_ds,
//						const PRISM_FLOAT_PRECISION* qyInd_d,
//						const PRISM_FLOAT_PRECISION* qxInd_d,
//						const PRISM_CUDA_COMPLEX_FLOAT* prop_d,
//						const size_t& ay,
//						const size_t& ax,
//						const size_t dimj,
//						const size_t dimi,
//						const PRISM_FLOAT_PRECISION* alphaInd_d,
//						const cufftHandle& plan,
//						cudaStream_t& stream){
						++currentBeam;
					}
				}
				cout << "GPU worker on stream #" << stream_count << " of GPU #" << GPU_num << "finished\n";
			}));
			++stream_count;
		}

		// launch CPU work
		vector<thread> workers_CPU;
		workers_CPU.reserve(pars.meta.NUM_THREADS); // prevents multiple reallocations



		for (auto &t:workers_GPU)t.join();
		for (auto &t:workers_CPU)t.join();


		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaFree(trans_d[g]));
			cudaErrchk(cudaFree(prop_d[g]));
			cudaErrchk(cudaFree(qxInd_d[g]));
			cudaErrchk(cudaFree(qyInd_d[g]));
			cudaErrchk(cudaFree(beamsIndex_d[g]));
		}

		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaSetDevice(s % pars.meta.NUM_GPUS));
			cudaErrchk(cudaFree(psi_ds[s]));
			cudaErrchk(cudaFree(psi_small_ds[s]));
			cufftErrchk(cufftDestroy(cufft_plan[s]));
			cufftErrchk(cufftDestroy(cufft_plan_small[s]));
		}



		// allocate memory on each GPU
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaMalloc((void **) &trans_d[g],      trans.size()      * sizeof(trans[0])));
			cudaErrchk(cudaMalloc((void **) &prop_d[g],       pars.prop.size()  * sizeof(pars.prop[0])));
			cudaErrchk(cudaMalloc((void **) &qxInd_d[g],      pars.qxInd.size() * sizeof(pars.qxInd[0])));
			cudaErrchk(cudaMalloc((void **) &qyInd_d[g],      pars.qyInd.size() * sizeof(pars.qyInd[0])));
			cudaErrchk(cudaMalloc((void **) &beamsIndex_d[g], pars.beamsIndex.size() * sizeof(pars.beamsIndex[0])));
		}


		// free pinned memory
		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaFreeHost(Scompact_slice_ph[s]));
		}
		cudaErrchk(cudaFreeHost(trans_ph));
		cudaErrchk(cudaFreeHost(prop_ph));
		cudaErrchk(cudaFreeHost(qxInd_ph));
		cudaErrchk(cudaFreeHost(qyInd_ph));
		cudaErrchk(cudaFreeHost(beamsIndex_ph));


		// destroy CUDA streams
		for (auto j = 0; j < total_num_streams; ++j){
			cudaSetDevice(j % pars.meta.NUM_GPUS);
			cudaErrchk(cudaStreamDestroy(streams[j]));
		}
	}
}