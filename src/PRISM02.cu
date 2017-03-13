#include "PRISM02.cuh"
#include "configure.h"
#include <mutex>
#include <thread>
#include <complex>
#include <vector>
#include "getWorkID.h"
#include "fftw3.h"
#include "defines.h"
namespace PRISM {
	using namespace std;
	void propagatePlaneWave_GPU(Parameters<PRISM_FLOAT_PRECISION> &pars,
						        Array3D<complex<PRISM_FLOAT_PRECISION> >& trans,
								size_t a0,
								Array2D<complex<PRISM_FLOAT_PRECISION> > &psi){

	}

	void fill_Scompact_GPU(Parameters <PRISM_FLOAT_PRECISION> &pars){

//		mutex fftw_plan_lock;

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
//		for (auto t = 0; t < total_num_streams; ++t){
//
//			workers_GPU.emplace_back([&pars](){
//				Array2D< complex<PRISM_FLOAT_PRECISION> > psi = zeros_ND<2, complex<PRISM_FLOAT_PRECISION> >({{pars.imageSize[0], pars.imageSize[1]}});
//				size_t Nstart, Nstop, ay, ax;
//				while (getWorkID(pars, Nstart, Nstop)){
//					while(Nstart != Nstop){
//						ay = Nstart / pars.xp.size();
//						ax = Nstart % pars.xp.size();
//						++Nstart;
//						propagatePlaneWave_GPU(pars, trans, Nstart, psi);
//					}
//				}
//
//			});
//		}

		for (auto &t:workers_GPU)t.join();
	}
}