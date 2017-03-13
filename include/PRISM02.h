//
// Created by AJ Pryor on 2/24/17.
//

#ifndef PRISM_PRISM02_H
#define PRISM_PRISM02_H
#include "params.h"
#include <iostream>
#include <vector>
#include <thread>
#include "fftw3.h"
#include <mutex>
#include "ArrayND.h"
#include <complex>
#include "utility.h"
#include "configure.h"
#include "getWorkID.h"
namespace PRISM {

	using namespace std;

//	template<class T>
//	using Array3D = PRISM::ArrayND<3, std::vector<T> >;
//	template<class T>
//	using Array2D = PRISM::ArrayND<2, std::vector<T> >;
//	template<class T>
//	using Array1D = PRISM::ArrayND<1, std::vector<T> >;
/*
	inline Array1D<PRISM_FLOAT_PRECISION> makeFourierCoords(const size_t &N, const PRISM_FLOAT_PRECISION &pixel_size) {
		Array1D<PRISM_FLOAT_PRECISION> result = zeros_ND<1, PRISM_FLOAT_PRECISION>({{N}});
		long long nc = (size_t) floor((PRISM_FLOAT_PRECISION) N / 2);

		PRISM_FLOAT_PRECISION dp = 1 / (N * pixel_size);
		for (auto i = 0; i < N; ++i) {
			result[(nc + (size_t) i) % N] = (i - nc) * dp;
		}
		return result;
	};
*/
	inline void propagatePlaneWave_CPU(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                        Array3D<complex<PRISM_FLOAT_PRECISION> >& trans,
	                        size_t a0,
	                        Array2D<complex<PRISM_FLOAT_PRECISION> > &psi,
	                        const PRISM_FFTW_PLAN &plan_forward,
	                        const PRISM_FFTW_PLAN &plan_inverse,
	                        mutex& fftw_plan_lock){
		psi[pars.beamsIndex[a0]] = 1;
		const PRISM_FLOAT_PRECISION N = (PRISM_FLOAT_PRECISION)psi.size();


		PRISM_FFTW_EXECUTE(plan_inverse);
		for (auto &i : psi)i /= N; // fftw scales by N, need to correct
		const complex<PRISM_FLOAT_PRECISION>* trans_t = &trans[0];
		for (auto a2 = 0; a2 < pars.numPlanes; ++a2){

			for (auto& p:psi)p*=(*trans_t++); // transmit
			PRISM_FFTW_EXECUTE(plan_forward); // FFT
			for (auto i = psi.begin(), j = pars.prop.begin(); i != psi.end();++i,++j)*i*=(*j); // propagate
			PRISM_FFTW_EXECUTE(plan_inverse); // IFFT
			for (auto &i : psi)i /= N; // fftw scales by N, need to correct
		}
		PRISM_FFTW_EXECUTE(plan_forward);

		Array2D< complex<PRISM_FLOAT_PRECISION> > psi_small = zeros_ND<2, complex<PRISM_FLOAT_PRECISION> >({{pars.qyInd.size(), pars.qxInd.size()}});


		unique_lock<mutex> gatekeeper(fftw_plan_lock);
		PRISM_FFTW_PLAN plan_final = PRISM_FFTW_PLAN_DFT_2D(psi_small.get_dimj(), psi_small.get_dimi(),
		                                                    reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi_small[0]),
		                                                    reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi_small[0]),
		                                                    FFTW_BACKWARD, FFTW_ESTIMATE);
		gatekeeper.unlock();
		for (auto y = 0; y < pars.qyInd.size(); ++y){
			for (auto x = 0; x < pars.qxInd.size(); ++x){
				psi_small.at(y,x) = psi.at(pars.qyInd[y], pars.qxInd[x]);
			}
		}
		PRISM_FFTW_EXECUTE(plan_final);
		gatekeeper.lock();
		PRISM_FFTW_DESTROY_PLAN(plan_final);
		gatekeeper.unlock();


		complex<PRISM_FLOAT_PRECISION>* S_t = &pars.Scompact[a0 * pars.Scompact.get_dimj() * pars.Scompact.get_dimi()];
		const PRISM_FLOAT_PRECISION N_small = (PRISM_FLOAT_PRECISION)psi_small.size();
		for (auto& i:psi_small) {
			*S_t++ = i/N_small;
		}
	};

	 void fill_Scompact_CPUOnly(Parameters<PRISM_FLOAT_PRECISION> &pars) {
		mutex fftw_plan_lock;

		const PRISM_FLOAT_PRECISION pi = acos(-1);
		const std::complex<PRISM_FLOAT_PRECISION> i(0, 1);
		pars.Scompact = zeros_ND<3, complex<PRISM_FLOAT_PRECISION> > ({{pars.numberBeams,pars.imageSize[0]/2, pars.imageSize[1]/2}});
		Array3D<complex<PRISM_FLOAT_PRECISION> > trans = zeros_ND<3, complex<PRISM_FLOAT_PRECISION> >(
				{{pars.pot.get_dimk(), pars.pot.get_dimj(), pars.pot.get_dimi()}});
		{
			auto p = pars.pot.begin();
			for (auto &j:trans)j = exp(i * pars.sigma * (*p++));
		}

//		for (auto& i:trans){i.real(1);i.imag(2);};


		vector<thread> workers;
		workers.reserve(pars.meta.NUM_THREADS); // prevents multiple reallocations
//		auto WORK_CHUNK_SIZE = ( (pars.numberBeams-1) / pars.meta.NUM_THREADS) + 1;
//		auto start = 0;
//		auto stop = start + WORK_CHUNK_SIZE;
		 setWorkStartStop(0, pars.numberBeams);
		 for (auto t = 0; t < pars.meta.NUM_THREADS; ++t){
			cout << "Launching thread #" << t << " to compute beams\n";
			workers.emplace_back([&pars, &fftw_plan_lock, &trans](){
				// allocate array for psi just once per thread
				Array2D< complex<PRISM_FLOAT_PRECISION> > psi = zeros_ND<2, complex<PRISM_FLOAT_PRECISION> >({{pars.imageSize[0], pars.imageSize[1]}});

				unique_lock<mutex> gatekeeper(fftw_plan_lock);
				PRISM_FFTW_PLAN plan_forward = PRISM_FFTW_PLAN_DFT_2D(psi.get_dimj(), psi.get_dimi(),
				                                          reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi[0]),
				                                          reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi[0]),
				                                          FFTW_FORWARD, FFTW_ESTIMATE);
				PRISM_FFTW_PLAN plan_inverse = PRISM_FFTW_PLAN_DFT_2D(psi.get_dimj(), psi.get_dimi(),
				                                          reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi[0]),
				                                          reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi[0]),
				                                          FFTW_BACKWARD, FFTW_ESTIMATE);
				gatekeeper.unlock(); // unlock it so we only block as long as necessary to deal with plans
				size_t currentBeam, stop;
				while (getWorkID(pars, currentBeam, stop)) { // synchronously get work assignment
					while (currentBeam != stop) {
						// re-zero psi each iteration
						memset((void *) &psi[0], 0, psi.size() * sizeof(complex<PRISM_FLOAT_PRECISION>));
						propagatePlaneWave_CPU(pars, trans, currentBeam, psi, plan_forward, plan_inverse, fftw_plan_lock);
						++currentBeam;
					}
				}
				// clean up
				gatekeeper.lock();
				PRISM_FFTW_DESTROY_PLAN(plan_forward);
				PRISM_FFTW_DESTROY_PLAN(plan_inverse);
				gatekeeper.unlock();
			});
		}
		cout << "Waiting for threads...\n";
		for (auto& t:workers)t.join();
	}

	inline void PRISM02(Parameters<PRISM_FLOAT_PRECISION>& pars){
		// propagate plane waves to construct compact S-matrix

		cout << "Entering PRISM02" << endl;
		// constants
		constexpr PRISM_FLOAT_PRECISION m = 9.109383e-31;
		constexpr PRISM_FLOAT_PRECISION e = 1.602177e-19;
		constexpr PRISM_FLOAT_PRECISION c = 299792458;
		constexpr PRISM_FLOAT_PRECISION h = 6.62607e-34;
		const PRISM_FLOAT_PRECISION pi = acos(-1);
		const std::complex<PRISM_FLOAT_PRECISION> i(0, 1);

		// setup some coordinates
		pars.imageSize[0] = pars.pot.get_dimj();
		pars.imageSize[1] = pars.pot.get_dimi();
		Array1D<PRISM_FLOAT_PRECISION> qx = makeFourierCoords(pars.imageSize[1], pars.pixelSize[1]);
		Array1D<PRISM_FLOAT_PRECISION> qy = makeFourierCoords(pars.imageSize[0], pars.pixelSize[0]);

		pair< Array2D<PRISM_FLOAT_PRECISION>, Array2D<PRISM_FLOAT_PRECISION> > mesh = meshgrid(qx,qy);
		pars.qxa = mesh.first;
		pars.qya = mesh.second;
		Array2D<PRISM_FLOAT_PRECISION> q2(pars.qya);
		transform(pars.qxa.begin(), pars.qxa.end(),
		          pars.qya.begin(), q2.begin(), [](const PRISM_FLOAT_PRECISION& a, const PRISM_FLOAT_PRECISION& b){
					return a*a + b*b;
				});

		// get qMax
		pars.qMax = 0;
		{
			PRISM_FLOAT_PRECISION qx_max;
			PRISM_FLOAT_PRECISION qy_max;
			for (auto i = 0; i < qx.size(); ++i) {
				qx_max = ( abs(qx[i]) > qx_max) ? abs(qx[i]) : qx_max;
				qy_max = ( abs(qy[i]) > qy_max) ? abs(qy[i]) : qy_max;
			}
			pars.qMax = min(qx_max, qy_max) / 2;
		}

		pars.qMask = zeros_ND<2, unsigned int>({{pars.imageSize[1], pars.imageSize[0]}});
		{
			long offset_x = pars.qMask.get_dimi()/4;
			long offset_y = pars.qMask.get_dimj()/4;
			long ndimy = (long)pars.qMask.get_dimj();
			long ndimx = (long)pars.qMask.get_dimi();
			for (long y = 0; y < pars.qMask.get_dimj() / 2; ++y) {
				for (long x = 0; x < pars.qMask.get_dimi() / 2; ++x) {
					pars.qMask.at( ((y-offset_y) % ndimy + ndimy) % ndimy,
					               ((x-offset_x) % ndimx + ndimx) % ndimx) = 1;
				}
			}
		}

		// build propagators
		pars.prop     = zeros_ND<2, std::complex<PRISM_FLOAT_PRECISION> >({{pars.imageSize[1], pars.imageSize[0]}});
		pars.propBack = zeros_ND<2, std::complex<PRISM_FLOAT_PRECISION> >({{pars.imageSize[1], pars.imageSize[0]}});
		for (auto y = 0; y < pars.qMask.get_dimj(); ++y) {
			for (auto x = 0; x < pars.qMask.get_dimi(); ++x) {
				if (pars.qMask.at(y,x)==1)
				{
					pars.prop.at(y,x)     = exp(-i * pi * complex<PRISM_FLOAT_PRECISION>(pars.lambda, 0) *
					                            complex<PRISM_FLOAT_PRECISION >(pars.meta.sliceThickness, 0) *
					                            complex<PRISM_FLOAT_PRECISION>(q2.at(y, x), 0));
					pars.propBack.at(y,x) = exp(i * pi * complex<PRISM_FLOAT_PRECISION>(pars.lambda, 0) *
					                            complex<PRISM_FLOAT_PRECISION>(pars.meta.cellDim[0], 0) *
					                            complex<PRISM_FLOAT_PRECISION>(q2.at(y, x), 0));
				}
			}
		}

		Array1D<PRISM_FLOAT_PRECISION> xv = makeFourierCoords(pars.imageSize[1], (PRISM_FLOAT_PRECISION)1/pars.imageSize[1]);
		Array1D<PRISM_FLOAT_PRECISION> yv = makeFourierCoords(pars.imageSize[0], (PRISM_FLOAT_PRECISION)1/pars.imageSize[0]);
		pair< Array2D<PRISM_FLOAT_PRECISION>, Array2D<PRISM_FLOAT_PRECISION> > mesh_a = meshgrid(xv, yv);

		// create beam mask and count beams
		PRISM::Array2D<unsigned int> mask;
		mask = zeros_ND<2, unsigned int>({{pars.imageSize[1], pars.imageSize[0]}});
		pars.numberBeams = 0;
		long interp_f = (long)pars.meta.interpolationFactor;
		for (auto y = 0; y < pars.qMask.get_dimj(); ++y) {
			for (auto x = 0; x < pars.qMask.get_dimi(); ++x) {
				if (q2.at(y,x) < pow(pars.meta.alphaBeamMax / pars.lambda,2) &&
				    pars.qMask.at(y,x)==1 &&
				    (long)round(mesh_a.first.at(y,x))  % interp_f == 0 &&
				    (long)round(mesh_a.second.at(y,x)) % interp_f == 0){
					mask.at(y,x)=1;
					++pars.numberBeams;
				}
			}
		}

		pars.beams     = zeros_ND<2, PRISM_FLOAT_PRECISION>({{pars.imageSize[0], pars.imageSize[1]}});
		{
			int beam_count = 1;
			for (auto y = 0; y < pars.qMask.get_dimj(); ++y) {
				for (auto x = 0; x < pars.qMask.get_dimi(); ++x) {
					if (mask.at(y,x)==1){
						pars.beamsIndex.push_back((size_t)y*pars.qMask.get_dimj() + (size_t)x);
						pars.beams.at(y,x) = beam_count++;
					}
				}
			}
		}

		// TODO: ensure this block is correct for arbitrary dimension
		// get the indices for the compact S-matrix
		pars.qxInd = zeros_ND<1, size_t>({{pars.imageSize[1]/2}});
		pars.qyInd = zeros_ND<1, size_t>({{pars.imageSize[0]/2}});
		{
			long n_0        = pars.imageSize[0];
			long n_1        = pars.imageSize[1];
			long n_half0    = pars.imageSize[0]/2;
			long n_half1    = pars.imageSize[1]/2;
			long n_quarter0 = pars.imageSize[0]/4;
			long n_quarter1 = pars.imageSize[1]/4;
			for (auto i = 0; i < n_quarter0; ++i) {
				pars.qyInd[i] = i;
				pars.qyInd[i + n_quarter0] = (i-n_quarter0) + n_0;
			}
			for (auto i = 0; i < n_quarter1; ++i) {
				pars.qxInd[i] = i;
				pars.qxInd[i + n_quarter1] = (i-n_quarter1) + n_1;
			}
		}

		cout << "Computing compact S matrix" << endl;

		// populate compact-S matrix
		fill_Scompact(pars);

		// downsample Fourier components by x2 to match output
		pars.imageSizeOutput = pars.imageSize;
		pars.imageSizeOutput[0]/=2;
		pars.imageSizeOutput[1]/=2;
		pars.pixelSizeOutput = pars.pixelSize;
		pars.pixelSizeOutput[0]*=2;
		pars.pixelSizeOutput[1]*=2;

		pars.qxaOutput   = zeros_ND<2, PRISM_FLOAT_PRECISION>({{pars.qyInd.size(), pars.qxInd.size()}});
		pars.qyaOutput   = zeros_ND<2, PRISM_FLOAT_PRECISION>({{pars.qyInd.size(), pars.qxInd.size()}});
		pars.beamsOutput = zeros_ND<2, PRISM_FLOAT_PRECISION>({{pars.qyInd.size(), pars.qxInd.size()}});
		for (auto y = 0; y < pars.qyInd.size(); ++y){
			for (auto x = 0; x < pars.qxInd.size(); ++x){
				pars.qxaOutput.at(y,x)   = pars.qxa.at(pars.qyInd[y],pars.qxInd[x]);
				pars.qyaOutput.at(y,x)   = pars.qya.at(pars.qyInd[y],pars.qxInd[x]);
				pars.beamsOutput.at(y,x) = pars.beams.at(pars.qyInd[y],pars.qxInd[x]);
			}
		}
	}


}
#endif //PRISM_PRISM02_H
