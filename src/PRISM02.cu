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
	                            PRISM_CUDA_COMPLEX_FLOAT* psi_ds,
	                            const PRISM_FLOAT_PRECISION* qya_d,
	                            const PRISM_FLOAT_PRECISION* qxa_d,
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

		const int total_num_streams = pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU;
		vector<thread> workers_GPU;
		workers_GPU.reserve(total_num_streams); // prevents multiple reallocations
		setWorkStartStop(0, pars.numberBeams);
		for (auto t = 0; t < total_num_streams; ++t){

			workers_GPU.emplace_back([&pars, &trans](){
				Array2D< complex<PRISM_FLOAT_PRECISION> > psi = zeros_ND<2, complex<PRISM_FLOAT_PRECISION> >({{pars.imageSize[0], pars.imageSize[1]}});
				size_t currentBeam, stop;
				while (getWorkID(pars, currentBeam, stop)){
					while(currentBeam != stop){
						propagatePlaneWave_GPU(pars, trans, currentBeam, psi);
						++currentBeam;
					}
				}

			});
		}

		for (auto &t:workers_GPU)t.join();
	}
}