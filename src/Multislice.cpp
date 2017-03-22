//
// Created by AJ Pryor on 3/13/17.
//

//
// Created by AJ Pryor on 3/6/17.
//

#include <iostream>
#include <thread>
#include <vector>
#include <mutex>
#include <numeric>
#include "configure.h"
#include <numeric>
#include "meta.h"
#include "ArrayND.h"
#include "params.h"
#include "utility.h"
#include "fftw3.h"
#include "getWorkID.h"
#include "Multislice.h"
namespace PRISM{
	using namespace std;
	static const PRISM_FLOAT_PRECISION pi = acos(-1);
	static const std::complex<PRISM_FLOAT_PRECISION> i(0, 1);

	void setupCoordinates_multislice(Parameters<PRISM_FLOAT_PRECISION>& pars){

		// setup coordinates and build propagators
		Array1D<PRISM_FLOAT_PRECISION> xR = zeros_ND<1, PRISM_FLOAT_PRECISION>({{2}});
		xR[0] = 0.1 * pars.meta.cellDim[2];
		xR[1] = 0.9 * pars.meta.cellDim[2];
		Array1D<PRISM_FLOAT_PRECISION> yR = zeros_ND<1, PRISM_FLOAT_PRECISION>({{2}});
		yR[0] = 0.1 * pars.meta.cellDim[1];
		yR[1] = 0.9 * pars.meta.cellDim[1];
		vector<PRISM_FLOAT_PRECISION> xp_d = vecFromRange(xR[0] + pars.meta.dxy / 2, pars.meta.dxy, xR[1] - pars.meta.dxy / 2);
		vector<PRISM_FLOAT_PRECISION> yp_d = vecFromRange(yR[0] + pars.meta.dxy / 2, pars.meta.dxy, yR[1] - pars.meta.dxy / 2);

		Array1D<PRISM_FLOAT_PRECISION> xp(xp_d, {{xp_d.size()}});
		Array1D<PRISM_FLOAT_PRECISION> yp(yp_d, {{yp_d.size()}});
		pars.xp = xp;
		pars.yp = yp;
		pars.imageSize[0] = pars.pot.get_dimj();
		pars.imageSize[1] = pars.pot.get_dimi();
		Array1D<PRISM_FLOAT_PRECISION> qx = makeFourierCoords(pars.imageSize[1], pars.pixelSize[1]);
		Array1D<PRISM_FLOAT_PRECISION> qy = makeFourierCoords(pars.imageSize[0], pars.pixelSize[0]);

		pair< Array2D<PRISM_FLOAT_PRECISION>, Array2D<PRISM_FLOAT_PRECISION> > mesh = meshgrid(qy,qx);
		pars.qya = mesh.first;
		pars.qxa = mesh.second;
		Array2D<PRISM_FLOAT_PRECISION> q2(pars.qya);
		transform(pars.qxa.begin(), pars.qxa.end(),
		          pars.qya.begin(), q2.begin(), [](const PRISM_FLOAT_PRECISION& a, const PRISM_FLOAT_PRECISION& b){
					return a*a + b*b;
				});
		Array2D<PRISM_FLOAT_PRECISION> q1(q2);
		pars.q2 = q2;
		pars.q1 = q1;
		for (auto& q : pars.q1)q=sqrt(q);

		// get qMax
		pars.qMax = 0;
		{
			PRISM_FLOAT_PRECISION qx_max;
			PRISM_FLOAT_PRECISION qy_max;
			for (auto ii = 0; ii < qx.size(); ++ii) {
				qx_max = ( abs(qx[ii]) > qx_max) ? abs(qx[ii]) : qx_max;
				qy_max = ( abs(qy[ii]) > qy_max) ? abs(qy[ii]) : qy_max;
			}
			pars.qMax = min(qx_max, qy_max) / 2;
		}

		pars.qMask = zeros_ND<2, unsigned int>({{pars.imageSize[0], pars.imageSize[1]}});
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
		pars.prop     = zeros_ND<2, std::complex<PRISM_FLOAT_PRECISION> >({{pars.imageSize[0], pars.imageSize[1]}});
		pars.propBack = zeros_ND<2, std::complex<PRISM_FLOAT_PRECISION> >({{pars.imageSize[0], pars.imageSize[1]}});
		for (auto y = 0; y < pars.qMask.get_dimj(); ++y) {
			for (auto x = 0; x < pars.qMask.get_dimi(); ++x) {
				if (pars.qMask.at(y,x)==1)
				{
					pars.prop.at(y,x)     = exp(-i * pi * complex<PRISM_FLOAT_PRECISION>(pars.lambda, 0) *
					                            complex<PRISM_FLOAT_PRECISION>(pars.meta.sliceThickness, 0) *
					                            complex<PRISM_FLOAT_PRECISION>(pars.q2.at(y, x), 0));
					pars.propBack.at(y,x) = exp(i * pi * complex<PRISM_FLOAT_PRECISION>(pars.lambda, 0) *
					                            complex<PRISM_FLOAT_PRECISION>(pars.meta.cellDim[0], 0) *
					                            complex<PRISM_FLOAT_PRECISION>(pars.q2.at(y, x), 0));
				}
			}
		}

	}

	void setupDetector_multislice(Parameters<PRISM_FLOAT_PRECISION>& pars){
		pars.dr = (PRISM_FLOAT_PRECISION)2.5 / 1000;
		pars.alphaMax = pars.qMax * pars.lambda;
		vector<PRISM_FLOAT_PRECISION> detectorAngles_d = vecFromRange(pars.dr / 2, pars.dr, pars.alphaMax - pars.dr / 2);
		Array1D<PRISM_FLOAT_PRECISION> detectorAngles(detectorAngles_d, {{detectorAngles_d.size()}});
		pars.detectorAngles = detectorAngles;
		pars.Ndet = pars.detectorAngles.size();
		Array2D<PRISM_FLOAT_PRECISION> alpha = pars.q1 * pars.lambda;
		pars.alphaInd = (alpha + pars.dr/2) / pars.dr;
		for (auto& q : pars.alphaInd) q = std::round(q);
		pars.dq = (pars.qxa.at(0, 1) + pars.qya.at(1, 0)) / 2;
	}

	void setupProbes_multislice(Parameters<PRISM_FLOAT_PRECISION>& pars){

		PRISM_FLOAT_PRECISION qProbeMax = pars.probeSemiangle/ pars.lambda; // currently a single semiangle
		pars.psiProbeInit = zeros_ND<2, complex<PRISM_FLOAT_PRECISION> >({{pars.q1.get_dimj(), pars.q1.get_dimi()}});
		Array2D<complex<PRISM_FLOAT_PRECISION> > psi;
		psi = zeros_ND<2, complex<PRISM_FLOAT_PRECISION> >({{pars.q1.get_dimj(), pars.q1.get_dimi()}});
		transform(pars.psiProbeInit.begin(), pars.psiProbeInit.end(),
		          pars.q1.begin(), pars.psiProbeInit.begin(),
		          [&pars, &qProbeMax](std::complex<PRISM_FLOAT_PRECISION> &a, PRISM_FLOAT_PRECISION &q1_t) {
			          a.real(erf((qProbeMax - q1_t) / (0.5 * pars.dq)) * 0.5 + 0.5);
			          a.imag(0);
			          return a;
		          });

		transform(pars.psiProbeInit.begin(), pars.psiProbeInit.end(),
		          pars.q2.begin(), pars.psiProbeInit.begin(),
		          [&pars](std::complex<PRISM_FLOAT_PRECISION> &a, PRISM_FLOAT_PRECISION &q2_t) {
			          a = a * exp(-i * pi * pars.lambda * pars.probeDefocus * q2_t); // TODO: fix hardcoded length-1 defocus
			          return a;
		          });
		PRISM_FLOAT_PRECISION norm_constant = sqrt(accumulate(pars.psiProbeInit.begin(), pars.psiProbeInit.end(),
		                                                      (PRISM_FLOAT_PRECISION)0.0, [](PRISM_FLOAT_PRECISION accum, std::complex<PRISM_FLOAT_PRECISION> &a) {
					return accum + abs(a) * abs(a);
				})); // make sure to initialize with 0.0 and NOT 0 or it won't be a float and answer will be wrong
		PRISM_FLOAT_PRECISION a = 0;
		for (auto &i : pars.psiProbeInit) { a += i.real(); };
		transform(pars.psiProbeInit.begin(), pars.psiProbeInit.end(),
		          pars.psiProbeInit.begin(), [&norm_constant](std::complex<PRISM_FLOAT_PRECISION> &a) {
					return a / norm_constant;
				});
	}

	void createTransmission(Parameters<PRISM_FLOAT_PRECISION>& pars){
		pars.transmission = zeros_ND<3, complex<PRISM_FLOAT_PRECISION> >(
				{{pars.pot.get_dimk(), pars.pot.get_dimj(), pars.pot.get_dimi()}});
		{
			auto p = pars.pot.begin();
			for (auto &j:pars.transmission)j = exp(i * pars.sigma * (*p++));
		}
	}

	void createStack(Parameters<PRISM_FLOAT_PRECISION>& pars){
		pars.stack = zeros_ND<4, PRISM_FLOAT_PRECISION>({{pars.yp.size(), pars.xp.size(), pars.Ndet, 1}});
	}

	void formatOutput_CPU_integrate(Parameters<PRISM_FLOAT_PRECISION>& pars,
	                                       Array2D< complex<PRISM_FLOAT_PRECISION> >& psi,
	                                       const Array2D<PRISM_FLOAT_PRECISION> &alphaInd,
	                                       const size_t& ay,
	                                       const size_t& ax){

		Array2D<PRISM_FLOAT_PRECISION> intOutput = zeros_ND<2, PRISM_FLOAT_PRECISION>({{psi.get_dimj(), psi.get_dimi()}});
		auto psi_ptr = psi.begin();
		for (auto& j:intOutput) j = pow(abs(*psi_ptr++),2);
		//update stack -- ax,ay are unique per thread so this write is thread-safe without a lock
		auto idx = alphaInd.begin();
		for (auto counts = intOutput.begin(); counts != intOutput.end(); ++counts){
			if (*idx <= pars.Ndet){
				pars.stack.at(ay,ax,(*idx)-1, 0) += *counts;
			}
			++idx;
		};
	}
	void getMultisliceProbe_CPU(Parameters<PRISM_FLOAT_PRECISION>& pars,
	                                   Array3D<complex<PRISM_FLOAT_PRECISION> >& trans,
	                                   const Array2D<complex<PRISM_FLOAT_PRECISION> >& psiProbeInit,
	                                   const size_t& ay,
	                                   const size_t& ax,
	                                   const Array2D<PRISM_FLOAT_PRECISION> &alphaInd){
		static mutex fftw_plan_lock; // for synchronizing access to shared FFTW resources

		// populates the output stack for Multislice simulation using the CPU. The number of
		// threads used is determined by pars.meta.NUM_THREADS

		Array2D<complex<PRISM_FLOAT_PRECISION> > psi(pars.psiProbeInit);


		// fftw_execute is the only thread-safe function in the library, so we need to synchronize access
		// to the plan creation methods
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
		{
			auto qxa_ptr = pars.qxa.begin();
			auto qya_ptr = pars.qya.begin();
			for (auto& p:psi)p*=exp(-2 * pi * i * ( (*qxa_ptr++)*pars.xp[ax] +
			                                        (*qya_ptr++)*pars.yp[ay]));
		}


		for (auto a2 = 0; a2 < pars.numPlanes; ++a2){
			PRISM_FFTW_EXECUTE(plan_inverse);
			complex<PRISM_FLOAT_PRECISION>* t_ptr = &pars.transmission[a2 * pars.transmission.get_dimj() * pars.transmission.get_dimi()];
			for (auto& p:psi)p *= (*t_ptr++); // transmit
			PRISM_FFTW_EXECUTE(plan_forward);
			auto p_ptr = pars.prop.begin();
			for (auto& p:psi)p *= (*p_ptr++); // propagate
			for (auto& p:psi)p /= psi.size(); // scale FFT

		}
		gatekeeper.lock();
		PRISM_FFTW_DESTROY_PLAN(plan_forward);
		PRISM_FFTW_DESTROY_PLAN(plan_inverse);
		gatekeeper.unlock();
		formatOutput_CPU(pars, psi, pars.alphaInd, ay, ax);
	}

	void buildMultisliceOutput_CPUOnly(Parameters<PRISM_FLOAT_PRECISION>& pars){
		cout << "CPU version" << endl;
		vector<thread> workers;
		workers.reserve(pars.meta.NUM_THREADS); // prevents multiple reallocations
		PRISM_FFTW_INIT_THREADS();
		PRISM_FFTW_PLAN_WITH_NTHREADS(pars.meta.NUM_THREADS);
		setWorkStartStop(0, pars.xp.size() * pars.yp.size());
		for (auto t = 0; t < pars.meta.NUM_THREADS; ++t){
			cout << "Launching CPU worker #" << t << '\n';

			// emplace_back is better whenever constructing a new object
			workers.emplace_back(thread([&pars, t]() {
				size_t Nstart, Nstop, ay, ax;
				while (getWorkID(pars, Nstart, Nstop)) { // synchronously get work assignment
					while (Nstart != Nstop) {
						ay = Nstart / pars.xp.size();
						ax = Nstart % pars.xp.size();
						getMultisliceProbe_CPU(pars, pars.transmission, pars.psiProbeInit, ay, ax, pars.alphaInd);
						++Nstart;
					}
				}
				cout << "CPU worker #" << t << " finished\n";
			}));
		}
		for (auto& t:workers)t.join();
		PRISM_FFTW_CLEANUP_THREADS();
	};


	void Multislice(Parameters<PRISM_FLOAT_PRECISION>& pars){
		using namespace std;

		// setup coordinates and build propagators
		setupCoordinates_multislice(pars);

		// setup detector coordinates and angles
		setupDetector_multislice(pars);

		// create initial probes
		setupProbes_multislice(pars);

		// create transmission array
		createTransmission(pars);

		// initialize output stack
		createStack(pars);

		// create the output
		buildMultisliceOutput(pars);
	}
}
