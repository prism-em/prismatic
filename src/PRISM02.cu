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
#include "utility.cuh"
#define PI 3.14159265359
#define BLOCK_SIZE1D 1024

namespace PRISM {
	using namespace std;



	void propagatePlaneWave_GPU(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                            PRISM_CUDA_COMPLEX_FLOAT* trans_d,
	                            PRISM_CUDA_COMPLEX_FLOAT* psi_d,
	                            PRISM_CUDA_COMPLEX_FLOAT* psi_small_d,
	                            complex<PRISM_FLOAT_PRECISION>* Scompact_slice_ph,
	                            const size_t* qyInd_d,
	                            const size_t* qxInd_d,
	                            const PRISM_CUDA_COMPLEX_FLOAT* prop_d,
	                            const size_t* beamsIndex,
	                            const size_t& beamNumber,
	                            const cufftHandle& plan,
	                            const cufftHandle& plan_small,
	                            cudaStream_t& stream){


		if (beamNumber==0) {
			cout << "TEST" << endl;
			size_t ans_sizet;
			cudaMemcpy(&ans_sizet, qyInd_d, sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qyInd_d[0] = " << ans_sizet << endl;
			cudaMemcpy(&ans_sizet, (qyInd_d + 1), sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qyInd_d[1] = " << ans_sizet << endl;
			cudaMemcpy(&ans_sizet, (qyInd_d + 15), sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qyInd_d[15] = " << ans_sizet << endl;
			cudaMemcpy(&ans_sizet, qxInd_d, sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qxInd_d[0] = " << ans_sizet << endl;
			cudaMemcpy(&ans_sizet, (qxInd_d + 1), sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qxInd_d[1] = " << ans_sizet << endl;
			cudaMemcpy(&ans_sizet, (qxInd_d + 15), sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qxInd_d[15] = " << ans_sizet << endl;


		}
		if (beamNumber==0) {
			cout << "block 0 " << endl;
			size_t ans_sizet;
			cudaMemcpy(&ans_sizet, qyInd_d, sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qyInd_d[0] = " << ans_sizet << endl;
			cudaMemcpy(&ans_sizet, (qyInd_d + 1), sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qyInd_d[1] = " << ans_sizet << endl;
			cudaMemcpy(&ans_sizet, (qyInd_d + 15), sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qyInd_d[15] = " << ans_sizet << endl;
			cudaMemcpy(&ans_sizet, qxInd_d, sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qxInd_d[0] = " << ans_sizet << endl;
			cudaMemcpy(&ans_sizet, (qxInd_d + 1), sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qxInd_d[1] = " << ans_sizet << endl;
			cudaMemcpy(&ans_sizet, (qxInd_d + 15), sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qxInd_d[15] = " << ans_sizet << endl;


		}
		const size_t psi_size = pars.imageSize[0] * pars.imageSize[1];
		initializePsi_oneNonzero<<< (psi_size - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream>>>(psi_d, psi_size, pars.beamsIndex[beamNumber]);

		if (beamNumber==0) {
			cout << "pars.beamsIndex[beamNumber]) = " << pars.beamsIndex[beamNumber] << endl;
			complex<double> ans;
			cudaMemcpy(&ans, psi_d, sizeof(ans), cudaMemcpyDeviceToHost);
			cout << "ans[0] = " << ans << endl;
			cudaMemcpy(&ans, (psi_d + 1), sizeof(ans), cudaMemcpyDeviceToHost);
			cout << "ans[1] = " << ans << endl;
			cudaMemcpy(&ans, (psi_d + 15), sizeof(ans), cudaMemcpyDeviceToHost);
			cout << "ans[15] = " << ans << endl;
		}
		if (beamNumber==0) {
			cout << "block 1 " << endl;
			size_t ans_sizet;
			cudaMemcpy(&ans_sizet, qyInd_d, sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qyInd_d[0] = " << ans_sizet << endl;
			cudaMemcpy(&ans_sizet, (qyInd_d + 1), sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qyInd_d[1] = " << ans_sizet << endl;
			cudaMemcpy(&ans_sizet, (qyInd_d + 15), sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qyInd_d[15] = " << ans_sizet << endl;
			cudaMemcpy(&ans_sizet, qxInd_d, sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qxInd_d[0] = " << ans_sizet << endl;
			cudaMemcpy(&ans_sizet, (qxInd_d + 1), sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qxInd_d[1] = " << ans_sizet << endl;
			cudaMemcpy(&ans_sizet, (qxInd_d + 15), sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qxInd_d[15] = " << ans_sizet << endl;


		}
		for (auto planeNum = 0; planeNum < 1; ++planeNum) {
			cufftErrchk(cufftExecC2C(plan, &psi_d[0], &psi_d[0], CUFFT_INVERSE));
			//multiply_inplace<<<(psi_size-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_d, &trans_d[planeNum*psi_size], psi_size);
			multiply_cx<<<(psi_size-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_d, &trans_d[planeNum*psi_size], psi_size);
//			multiply_inplace<<<(psi_size-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_d, &trans_d[planeNum*psi_size], psi_size);
			divide_inplace<<<(psi_size-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_d, PRISM_MAKE_CU_COMPLEX(psi_size, 0), psi_size);
			cufftErrchk(cufftExecC2C(plan, &psi_d[0], &psi_d[0], CUFFT_FORWARD));
			//multiply_inplace<<<(psi_size-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_d, prop_d, psi_size);
			multiply_cx<<<(psi_size-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_d, prop_d, psi_size);

		}

		{
			if (beamNumber == 0) {
				vector<complex<PRISM_FLOAT_PRECISION> > debug_cx(pars.imageSize[0] * pars.imageSize[1], 0);
				cout << "debug_cx.size() = " << debug_cx.size() << endl;
				cudaMemcpy(&debug_cx[0], (psi_d), psi_size * sizeof(complex<PRISM_FLOAT_PRECISION>),
				           cudaMemcpyDeviceToHost);
				Array2D<PRISM_FLOAT_PRECISION> debug = zeros_ND<2, PRISM_FLOAT_PRECISION>(
						{{pars.imageSize[0], pars.imageSize[1]}});
				for (auto i = 0; i < psi_size; ++i) {
					debug[i] = abs(debug_cx[i]);
				}
				debug.toMRC_f("debugPSI.mrc");
			}
		}

//		for (auto planeNum = 0; planeNum < pars.numPlanes; ++planeNum) {
//			cufftErrchk(cufftExecC2C(plan, &psi_d[0], &psi_d[0], CUFFT_INVERSE));
//			multiply_inplace<<<(psi_size-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_d, &trans_d[planeNum*psi_size], psi_size);
//			cufftErrchk(cufftExecC2C(plan, &psi_d[0], &psi_d[0], CUFFT_FORWARD));
//			multiply_inplace<<<(psi_size-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_d, prop_d, psi_size);
//			divide_inplace<<<(psi_size-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>>(psi_d, psi_size, psi_size);
//		}
		if (beamNumber==0) {
			cout << "pars.beamsIndex[beamNumber]) = " << pars.beamsIndex[beamNumber] << endl;
			complex<PRISM_FLOAT_PRECISION> ans;
			cudaMemcpy(&ans, psi_d, sizeof(ans), cudaMemcpyDeviceToHost);
			cout << "after ans[0] = " << ans << endl;
			cudaMemcpy(&ans, (psi_d + 1), sizeof(ans), cudaMemcpyDeviceToHost);
			cout << "ans[1] = " << ans << endl;
			cudaMemcpy(&ans, (psi_d + 15), sizeof(ans), cudaMemcpyDeviceToHost);
			cout << "ans[15] = " << ans << endl;

			complex<PRISM_FLOAT_PRECISION> ans_sum = {0.0,0.0};
			for (auto i = 0; i < pars.imageSize[1] * pars.imageSize[0]; ++i){
				cudaMemcpy(&ans, (psi_d + i), sizeof(ans), cudaMemcpyDeviceToHost);
//				cout << "psi_d[" << i << "] = " << ans << endl;
				ans_sum += (complex<PRISM_FLOAT_PRECISION>)ans;
			}
			cout << "psi sum = " << ans_sum << endl;

//			complex<PRISM_FLOAT_PRECISION> ans_sum = {0.0f,0.0f};
//			for (auto i = 0; i < pars.imageSize[1] * pars.imageSize[0]; ++i){
//				cudaMemcpy(&ans, (psi_d + i), sizeof(ans), cudaMemcpyDeviceToHost);
//				cout << "psi_d[" << i << "] = " << ans << endl;
//				ans_sum += ans;
//			}
//			cout << "psi sum = " << ans_sum << endl;
		}

		if (beamNumber==0) {
			cout << "block 2 " << endl;
			size_t ans_sizet;
			cudaMemcpy(&ans_sizet, qyInd_d, sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qyInd_d[0] = " << ans_sizet << endl;
			cudaMemcpy(&ans_sizet, (qyInd_d + 1), sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qyInd_d[1] = " << ans_sizet << endl;
			cudaMemcpy(&ans_sizet, (qyInd_d + 15), sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qyInd_d[15] = " << ans_sizet << endl;
			cudaMemcpy(&ans_sizet, qxInd_d, sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qxInd_d[0] = " << ans_sizet << endl;
			cudaMemcpy(&ans_sizet, (qxInd_d + 1), sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qxInd_d[1] = " << ans_sizet << endl;
			cudaMemcpy(&ans_sizet, (qxInd_d + 15), sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qxInd_d[15] = " << ans_sizet << endl;


		}
		array_subset<<<(pars.qyInd.size()*pars.qxInd.size()-1) / BLOCK_SIZE1D + 1,BLOCK_SIZE1D, 0, stream>>> ( psi_d,
																					psi_small_d,
																					qyInd_d,
																					qxInd_d,
																					pars.imageSize[1],
																					pars.qyInd.size(),
																					pars.qxInd.size());
		if (beamNumber==0) {
			complex<PRISM_FLOAT_PRECISION> ans;
			cudaMemcpy(&ans, trans_d, sizeof(ans), cudaMemcpyDeviceToHost);
			cout << "trans_d[0] = " << ans << endl;
			cudaMemcpy(&ans, (trans_d + 1), sizeof(ans), cudaMemcpyDeviceToHost);
			cout << "trans_d[1] = " << ans << endl;
			cudaMemcpy(&ans, (trans_d + 15), sizeof(ans), cudaMemcpyDeviceToHost);
			cout << "trans_d[15] = " << ans << endl;
			cudaMemcpy(&ans, (trans_d + 3049200-1), sizeof(ans), cudaMemcpyDeviceToHost);
			cout << "trans_d[3049199] = " << ans << endl;

//			complex<PRISM_FLOAT_PRECISION> ans_sum = {0,0};
//			for (auto i = 0; i < pars.pot.size(); ++i){
//				cudaMemcpy(&ans, (trans_d + i), sizeof(ans), cudaMemcpyDeviceToHost);
//				ans_sum += ans;
//			}
//			for (auto i = 0; i < 100; ++i){
//				cudaMemcpy(&ans, (trans_d + i), sizeof(ans), cudaMemcpyDeviceToHost);
//				cout << "trans_d[" << i << "] = " << ans << endl;
//			}
//
//			cout << "trans_d sum = " << ans_sum << endl;
		}
		if (beamNumber==0) {
			complex<PRISM_FLOAT_PRECISION> ans;
			cudaMemcpy(&ans, prop_d, sizeof(ans), cudaMemcpyDeviceToHost);
			cout << "prop_d[0] = " << ans << endl;
			cudaMemcpy(&ans, (prop_d + 1), sizeof(ans), cudaMemcpyDeviceToHost);
			cout << "prop_d[1] = " << ans << endl;
			cudaMemcpy(&ans, (prop_d + 15), sizeof(ans), cudaMemcpyDeviceToHost);
			cout << "prop_d[15] = " << ans << endl;
//			complex<PRISM_FLOAT_PRECISION> ans_sum = {0,0};
//			for (auto i = 0; i < pars.prop.size(); ++i){
//				cudaMemcpy(&ans, (prop_d + i), sizeof(ans), cudaMemcpyDeviceToHost);
//				ans_sum += ans;
//			}
//			cout << "prop sum = " << ans_sum << endl;
		}

		if (beamNumber==0) {
			cout << "block 3 " << endl;
			size_t ans_sizet;
			cudaMemcpy(&ans_sizet, qyInd_d, sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qyInd_d[0] = " << ans_sizet << endl;
			cudaMemcpy(&ans_sizet, (qyInd_d + 1), sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qyInd_d[1] = " << ans_sizet << endl;
			cudaMemcpy(&ans_sizet, (qyInd_d + 15), sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qyInd_d[15] = " << ans_sizet << endl;
			cudaMemcpy(&ans_sizet, qxInd_d, sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qxInd_d[0] = " << ans_sizet << endl;
			cudaMemcpy(&ans_sizet, (qxInd_d + 1), sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qxInd_d[1] = " << ans_sizet << endl;
			cudaMemcpy(&ans_sizet, (qxInd_d + 15), sizeof(ans_sizet), cudaMemcpyDeviceToHost);
			cout << "qxInd_d[15] = " << ans_sizet << endl;


		}

		if (beamNumber==0) {
			complex<PRISM_FLOAT_PRECISION> ans;
//			for (auto i = 0; i < pars.imageSize[1] * pars.imageSize[0]; ++i) {
//				cudaMemcpy(&ans, psi_small_d + i, sizeof(ans), cudaMemcpyDeviceToHost);
//				cout << "ans[" <<i << "] = " << ans << endl;
//			}
			cudaMemcpy(&ans, (psi_small_d + 1), sizeof(ans), cudaMemcpyDeviceToHost);
			cout << "ans[1] = " << ans << endl;
			cudaMemcpy(&ans, (psi_small_d + 15), sizeof(ans), cudaMemcpyDeviceToHost);
			cout << "ans[15] = " << ans << endl;
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
		complex<PRISM_FLOAT_PRECISION> trans_sum = {0,0};
		for (auto& j : trans)trans_sum+=j;
		cout << "trans sum = " << trans_sum << endl;



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
			cufftErrchk(cufftPlan2d(&cufft_plan_small[j], pars.qyInd.size(), pars.qxInd.size(), CUFFT_C2C));
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


		complex<PRISM_FLOAT_PRECISION> trans_ph_sum = {0,0};
		for (auto jj = 0; jj < trans.size(); ++jj)trans_ph_sum += trans_ph[jj];
		cout << "trans_ph sum = " << trans_ph_sum << endl;
		for (auto i = 0; i < 10; ++i){
			cout << "qxInd_ph[i] = " << qxInd_ph[i] << endl;
			cout << "pars.qxInd[[i] = " << pars.qxInd[i] << endl;
		}
		// pointers to read-only GPU memory (one copy per GPU)
		PRISM_CUDA_COMPLEX_FLOAT *trans_d[pars.meta.NUM_GPUS];
		PRISM_CUDA_COMPLEX_FLOAT *prop_d[pars.meta.NUM_GPUS];
		size_t                   *qxInd_d[pars.meta.NUM_GPUS];
		size_t                   *qyInd_d[pars.meta.NUM_GPUS];
		size_t                   *beamsIndex_d[pars.meta.NUM_GPUS];

		// pointers to read/write GPU memory (one per stream)
		PRISM_CUDA_COMPLEX_FLOAT *psi_ds[total_num_streams];
		PRISM_CUDA_COMPLEX_FLOAT *psi_small_ds[total_num_streams];

		cout << "allocate gpu memory" << endl;
		// allocate memory on each GPU
		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g) {
			cudaErrchk(cudaSetDevice(g));
			cudaErrchk(cudaMalloc((void **) &trans_d[g],      trans.size()      * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMalloc((void **) &prop_d[g],       pars.prop.size()  * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMalloc((void **) &qxInd_d[g],      pars.qxInd.size() * sizeof(size_t)));
			cudaErrchk(cudaMalloc((void **) &qyInd_d[g],      pars.qyInd.size() * sizeof(size_t)));
			cudaErrchk(cudaMalloc((void **) &beamsIndex_d[g], pars.beamsIndex.size() * sizeof(size_t)));
		}

		// allocate memory per stream and 0 it
		for (auto s = 0; s < total_num_streams; ++s) {
			cudaErrchk(cudaSetDevice(s % pars.meta.NUM_GPUS));
			cudaErrchk(cudaMalloc((void **) &psi_ds[s],       pars.imageSize[0] * pars.imageSize[1] * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMalloc((void **) &psi_small_ds[s], pars.qxInd.size() * pars.qyInd.size() * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMemset(psi_ds[s], 0, pars.imageSize[0] * pars.imageSize[1] * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
			cudaErrchk(cudaMemset(psi_small_ds[s], 0, pars.qxInd.size() * pars.qyInd.size() * sizeof(PRISM_CUDA_COMPLEX_FLOAT)));
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
			trans_sum = {0,0};
			for (auto& j : trans)trans_sum+=j;
			cout << "trans sum again = " << trans_sum << endl;
			cudaErrchk(cudaMemcpyAsync(trans_d[g], &trans_ph[0],
			                           trans.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>), cudaMemcpyHostToDevice, streams[stream_id]));

			if (g==0) {
				complex<PRISM_FLOAT_PRECISION> ans;
				cudaMemcpy(&ans, trans_d[g], sizeof(ans), cudaMemcpyDeviceToHost);
				cout << "trans_d[0] = " << ans << endl;
				cudaMemcpy(&ans, (trans_d[g] + 1), sizeof(ans), cudaMemcpyDeviceToHost);
				cout << "trans_d[1] = " << ans << endl;
				cudaMemcpy(&ans, (trans_d[g] + 15), sizeof(ans), cudaMemcpyDeviceToHost);
				cout << "trans_d[15] = " << ans << endl;
				cudaMemcpy(&ans, (trans_d[g] + 3049200-1), sizeof(ans), cudaMemcpyDeviceToHost);
				cout << "trans_d[3049199] = " << ans << endl;

//				complex<PRISM_FLOAT_PRECISION> ans_sum = {0,0};
//				for (auto i = 0; i < trans.size(); ++i){
//					cudaMemcpy(&ans, (trans_d[g] + i), sizeof(ans), cudaMemcpyDeviceToHost);
//					ans_sum += ans;
//				}
//				for (auto i = 0; i < 100; ++i){
//					cudaMemcpy(&ans, (trans_d[g] + i), sizeof(ans), cudaMemcpyDeviceToHost);
//					cout << "trans_d[" << i << "] = " << ans << endl;
//				}
//
//				cout << "trans_d sum = " << ans_sum << endl;
			}





			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cout << "stream_id = " << stream_id << endl;
			cudaErrchk(cudaMemcpyAsync(prop_d[g], &prop_ph[0],
		                           pars.prop.size() * sizeof(std::complex<PRISM_FLOAT_PRECISION>), cudaMemcpyHostToDevice, streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cout << "stream_id = " << stream_id << endl;
			cout << "qxInd_ph[1] = " << qxInd_ph[1]<< endl;
			cout << "pars.qxInd[1] = " << pars.qxInd[1]<< endl;
			cudaErrchk(cudaMemcpyAsync(qxInd_d[g], &qxInd_ph[0],
			                           pars.qxInd.size() * sizeof(size_t), cudaMemcpyHostToDevice, streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cout << "stream_id = " << stream_id << endl;
			cudaErrchk(cudaMemcpyAsync(qyInd_d[g], &qyInd_ph[0],
			                           pars.qyInd.size() * sizeof(size_t), cudaMemcpyHostToDevice, streams[stream_id]));
			stream_id = (stream_id + pars.meta.NUM_GPUS) % total_num_streams;
			cout << "stream_id = " << stream_id << endl;
//			cudaErrchk(cudaMemcpy(beamsIndex_d[g], &beamsIndex_ph[0],
//			                           pars.beamsIndex.size() * sizeof(size_t), cudaMemcpyHostToDevice));
			cudaErrchk(cudaMemcpyAsync(beamsIndex_d[g], &beamsIndex_ph[0],
			                           pars.beamsIndex.size() * sizeof(size_t), cudaMemcpyHostToDevice, streams[stream_id]));
		}

		for (auto g = 0; g < pars.meta.NUM_GPUS; ++g){
			cudaSetDevice(g);
			cudaDeviceSynchronize();
		}







		// launch GPU work
		vector<thread> workers_GPU;
		workers_GPU.reserve(total_num_streams); // prevents multiple reallocations
		cout << "pars.numberBeams = " << pars.numberBeams<< endl;
		int stream_count = 0;
		setWorkStartStop(0, pars.numberBeams);
		for (auto t = 0; t < total_num_streams; ++t){

			int GPU_num = stream_count % pars.meta.NUM_GPUS; // determine which GPU handles this job
			cudaStream_t& current_stream = streams[stream_count];
			cout << "Launching GPU worker on stream #" << stream_count << " of GPU #" << GPU_num << '\n';

			// get pointers to the pre-copied arrays, making sure to get those on the current GPU
			PRISM_CUDA_COMPLEX_FLOAT *current_trans_d = trans_d[GPU_num];
			PRISM_CUDA_COMPLEX_FLOAT *current_prop_d  = prop_d[GPU_num];
			size_t                   *current_qxInd_d    = qxInd_d[GPU_num];
			size_t                   *current_qyInd_d    = qyInd_d[GPU_num];
			size_t                   *current_beamsIndex = beamsIndex_d[GPU_num];
			// get pointers to per-stream arrays
			PRISM_CUDA_COMPLEX_FLOAT *current_psi_ds             = psi_ds[stream_count];
			PRISM_CUDA_COMPLEX_FLOAT *current_psi_small_ds       = psi_small_ds[stream_count];
			cufftHandle& current_cufft_plan                      = cufft_plan[stream_count];
			cufftHandle& current_cufft_plan_small                = cufft_plan_small[stream_count];
			complex<PRISM_FLOAT_PRECISION > *current_S_slice_ph  = Scompact_slice_ph[stream_count];



			if (t==0) {
				size_t ans_sizet;
				cudaMemcpy(&ans_sizet, current_qyInd_d, sizeof(ans_sizet), cudaMemcpyDeviceToHost);
				cout << "current_qyInd_d[0] = " << ans_sizet << endl;
				cudaMemcpy(&ans_sizet, (current_qyInd_d + 1), sizeof(ans_sizet), cudaMemcpyDeviceToHost);
				cout << "current_qyInd_d[1] = " << ans_sizet << endl;
				cudaMemcpy(&ans_sizet, (current_qyInd_d + 15), sizeof(ans_sizet), cudaMemcpyDeviceToHost);
				cout << "current_qyInd_d[15] = " << ans_sizet << endl;
				cudaMemcpy(&ans_sizet, current_qyInd_d, sizeof(ans_sizet), cudaMemcpyDeviceToHost);
				cout << "current_qyInd_d[0] = " << ans_sizet << endl;
				cudaMemcpy(&ans_sizet, (current_qyInd_d + 1), sizeof(ans_sizet), cudaMemcpyDeviceToHost);
				cout << "current_qyInd_d[1] = " << ans_sizet << endl;
				cudaMemcpy(&ans_sizet, (current_qyInd_d + 15), sizeof(ans_sizet), cudaMemcpyDeviceToHost);
				cout << "current_qyInd_d[15] = " << ans_sizet << endl;


			}


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
//						cudaDeviceSynchronize();
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