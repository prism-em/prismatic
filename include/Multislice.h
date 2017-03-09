//
// Created by AJ Pryor on 3/6/17.
//

#ifndef PRISM_MULTISLICE_H
#define PRISM_MULTISLICE_H
#include <iostream>
#include <thread>
#include <vector>
#include <mutex>
#include "configure.h"
#include "meta.h"
#include "ArrayND.h"
#include "params.h"
#include "utility.h"
#include "fftw3.h"
namespace PRISM{
       using namespace std;

	inline void formatOutput_cpu_integrate(Parameters<PRISM_FLOAT_PRECISION>& pars,
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
	inline void getMultisliceProbe_cpu(Parameters<PRISM_FLOAT_PRECISION>& pars,
                                Array3D<complex<PRISM_FLOAT_PRECISION> >& trans,
                                const Array2D<complex<PRISM_FLOAT_PRECISION> >& PsiProbeInit,
                                const size_t& ay,
                                const size_t& ax,
                                const Array2D<PRISM_FLOAT_PRECISION> &alphaInd){
        static mutex fftw_plan_lock; // for synchronizing access to shared FFTW resources
		static const PRISM_FLOAT_PRECISION pi = acos(-1);
		static const std::complex<PRISM_FLOAT_PRECISION> i(0, 1);
		// populates the output stack for Multislice simulation using the CPU. The number of
		// threads used is determined by pars.meta.NUM_THREADS

		Array2D<complex<PRISM_FLOAT_PRECISION> > psi(PsiProbeInit);

		{
			auto qxa_ptr = pars.qxa.begin();
			auto qya_ptr = pars.qya.begin();
            for (auto& p:psi)p*=exp(-2 * pi * i * ( (*qxa_ptr++)*pars.xp[ax] +
                                                    (*qya_ptr++)*pars.yp[ay]));
		}

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

        for (auto a2 = 0; a2 < pars.numPlanes; ++a2){
	        PRISM_FFTW_EXECUTE(plan_inverse);
            complex<PRISM_FLOAT_PRECISION>* t_ptr = &trans[a2 * trans.get_dimj() * trans.get_dimi()];
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

        formatOutput_cpu(pars, psi, alphaInd, ay, ax);

//		Array2D<PRISM_FLOAT_PRECISION> intOutput = zeros_ND<2, PRISM_FLOAT_PRECISION>({{psi.get_dimj(), psi.get_dimi()}});
//        auto psi_ptr = psi.begin();
//        for (auto& j:intOutput) j = pow(abs(*psi_ptr++),2);
//        //update stack -- ax,ay are unique per thread so this write is thread-safe without a lock
//        auto idx = alphaInd.begin();
//        for (auto counts = intOutput.begin(); counts != intOutput.end(); ++counts){
//            if (*idx <= pars.Ndet){
//                pars.stack.at(ay,ax,(*idx)-1, 0) += *counts;
//            }
//            ++idx;
//        };


//
//		% Propgate through all potential planes
//		for a2 = 1:emdSTEM.numPlanes
//		psi = fft2(ifft2(psi).*trans(:,:,a2)).*emdSTEM.prop;
//		end
//
//		% Record output
//		emdSTEM.MULTIstack(a0,a1,:) = ...
//		accumarray(alphaIndsSub,abs(psi(alphaMask)).^2,[Ndet 1]);
//		using namespace std;
//		cout << "test cpu" << endl;
	}

	inline void buildMultisliceOutput_cpuOnly(Parameters<PRISM_FLOAT_PRECISION>& pars,
                                       Array3D<complex<PRISM_FLOAT_PRECISION> >& trans,
                                       Array2D<complex<PRISM_FLOAT_PRECISION> >& PsiProbeInit,
                                       Array2D<PRISM_FLOAT_PRECISION> &alphaInd){
		cout << "cpu version" << endl;
		vector<thread> workers;
		workers.reserve(pars.meta.NUM_THREADS); // prevents multiple reallocations
		auto WORK_CHUNK_SIZE = ((pars.yp.size() - 1) / pars.meta.NUM_THREADS) + 1; //TODO: divide work more generally than just splitting up by yp. If input isn't square this might not do a good job
		auto start = 0;
		auto stop = start + WORK_CHUNK_SIZE;
		while (start < pars.yp.size()) {
			cout << "Launching thread to compute all x-probe positions for y-probes "
				 << start << "/" << min(stop,pars.yp.size()) << '\n';
			// emplace_back is better whenever constructing a new object
			workers.emplace_back(thread([&pars, &trans,
												&alphaInd, &PsiProbeInit,
												start, stop]() {
				for (auto ay = start; ay < min((size_t) stop, pars.yp.size()); ++ay) {
					for (auto ax = 0; ax < pars.xp.size(); ++ax) {
						getMultisliceProbe_cpu(pars, trans, PsiProbeInit, ay, ax, alphaInd);
					}
				}
			}));
			start += WORK_CHUNK_SIZE;
			if (start >= pars.yp.size())break;
			stop += WORK_CHUNK_SIZE;
		}
		for (auto& t:workers)t.join();

	};


	inline void Multislice(Parameters<PRISM_FLOAT_PRECISION>& pars){
		using namespace std;
		static const PRISM_FLOAT_PRECISION pi = acos(-1);
		static const std::complex<PRISM_FLOAT_PRECISION> i(0, 1);

		const PRISM_FLOAT_PRECISION dxy = (PRISM_FLOAT_PRECISION)0.25 * 2; // TODO: move this

		// should move these elsewhere and in PRISM03
		pars.probeDefocusArray = zeros_ND<1, PRISM_FLOAT_PRECISION>({{1}});
		pars.probeSemiangleArray = zeros_ND<1, PRISM_FLOAT_PRECISION>({{1}});
		pars.probeXtiltArray = zeros_ND<1, PRISM_FLOAT_PRECISION>({{1}});
		pars.probeYtiltArray = zeros_ND<1, PRISM_FLOAT_PRECISION>({{1}});
		pars.probeDefocusArray[0] = (PRISM_FLOAT_PRECISION)0.0;
		pars.probeSemiangleArray[0] = (PRISM_FLOAT_PRECISION)20.0 / 1000;
		pars.probeXtiltArray[0] = (PRISM_FLOAT_PRECISION)0.0 / 1000;
		pars.probeYtiltArray[0] = (PRISM_FLOAT_PRECISION)0.0 / 1000;

		Array1D<PRISM_FLOAT_PRECISION> xR = zeros_ND<1, PRISM_FLOAT_PRECISION>({{2}});
		xR[0] = 0.1 * pars.meta.cellDim[2];
		xR[1] = 0.9 * pars.meta.cellDim[2];
		Array1D<PRISM_FLOAT_PRECISION> yR = zeros_ND<1, PRISM_FLOAT_PRECISION>({{2}});
		yR[0] = 0.1 * pars.meta.cellDim[1];
		yR[1] = 0.9 * pars.meta.cellDim[1];
		vector<PRISM_FLOAT_PRECISION> xp_d = vecFromRange(xR[0] + dxy / 2, dxy, xR[1] - dxy / 2);
		vector<PRISM_FLOAT_PRECISION> yp_d = vecFromRange(yR[0] + dxy / 2, dxy, yR[1] - dxy / 2);

		Array1D<PRISM_FLOAT_PRECISION> xp(xp_d, {{xp_d.size()}});
		Array1D<PRISM_FLOAT_PRECISION> yp(yp_d, {{yp_d.size()}});
		pars.xp = xp;
		pars.yp = yp;

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
		Array2D<PRISM_FLOAT_PRECISION> q1(q2);
		for (auto& q : q1)q=sqrt(q);

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
					                            complex<PRISM_FLOAT_PRECISION>(pars.meta.sliceThickness, 0) *
					                            complex<PRISM_FLOAT_PRECISION>(q2.at(y, x), 0));
					pars.propBack.at(y,x) = exp(i * pi * complex<PRISM_FLOAT_PRECISION>(pars.lambda, 0) *
					                            complex<PRISM_FLOAT_PRECISION>(pars.meta.cellDim[0], 0) *
					                            complex<PRISM_FLOAT_PRECISION>(q2.at(y, x), 0));
				}
			}
		}

		pars.dr = (PRISM_FLOAT_PRECISION)2.5 / 1000;
		pars.alphaMax = pars.qMax * pars.lambda;
		vector<PRISM_FLOAT_PRECISION> detectorAngles_d = vecFromRange(pars.dr / 2, pars.dr, pars.alphaMax - pars.dr / 2);
		Array1D<PRISM_FLOAT_PRECISION> detectorAngles(detectorAngles_d, {{detectorAngles_d.size()}});
		pars.detectorAngles = detectorAngles;
		Array2D<PRISM_FLOAT_PRECISION> alpha = q1 * pars.lambda;
//		Array2D<PRISM_FLOAT_PRECISION> alphaInd(q1); // copy constructor more efficient than assignment
//		transform(alphaInd.begin(), alphaInd.end(),
//				  alphaInd.begin(),
//				  [&pars](const PRISM_FLOAT_PRECISION &a) {
//					  return 1 + round((a * pars.lambda - pars.detectorAngles[0]) / pars.dr);
//				  });
//		transform(alphaInd.begin(), alphaInd.end(),
//				  alphaInd.begin(),
//				  [](const PRISM_FLOAT_PRECISION &a) { return a < 1 ? 1 : a; });
		Array2D<PRISM_FLOAT_PRECISION> alphaInd = (alpha + pars.dr/2) / pars.dr;
		for (auto& q : alphaInd) q = std::round(q);

		pars.Ndet = pars.detectorAngles.size();
		if (pars.probeSemiangleArray.size() > 1)throw std::domain_error("Currently only scalar probeSemiangleArray supported. Multiple inputs received.\n");
		PRISM_FLOAT_PRECISION qProbeMax = pars.probeSemiangleArray[0] / pars.lambda; // currently a single semiangle
		Array2D<complex<PRISM_FLOAT_PRECISION> > PsiProbeInit, psi;
		PsiProbeInit = zeros_ND<2, complex<PRISM_FLOAT_PRECISION> >({{q1.get_dimj(), q1.get_dimi()}});
		psi = zeros_ND<2, complex<PRISM_FLOAT_PRECISION> >({{q1.get_dimj(), q1.get_dimi()}});
		pars.dq = (pars.qxa.at(0, 1) + pars.qya.at(1, 0)) / 2;

		transform(PsiProbeInit.begin(), PsiProbeInit.end(),
		          q1.begin(), PsiProbeInit.begin(),
		          [&pars, &qProbeMax](std::complex<PRISM_FLOAT_PRECISION> &a, PRISM_FLOAT_PRECISION &q1_t) {
			          a.real(erf((qProbeMax - q1_t) / (0.5 * pars.dq)) * 0.5 + 0.5);
			          a.imag(0);
			          return a;
		          });


		transform(PsiProbeInit.begin(), PsiProbeInit.end(),
		          q2.begin(), PsiProbeInit.begin(),
		          [&pars](std::complex<PRISM_FLOAT_PRECISION> &a, PRISM_FLOAT_PRECISION &q2_t) {
			          a = a * exp(-i * pi * pars.lambda * pars.probeDefocusArray[0] * q2_t); // TODO: fix hardcoded length-1 defocus
			          return a;
		          });
		PRISM_FLOAT_PRECISION norm_constant = sqrt(accumulate(PsiProbeInit.begin(), PsiProbeInit.end(),
		                                                      (PRISM_FLOAT_PRECISION)0.0, [](PRISM_FLOAT_PRECISION accum, std::complex<PRISM_FLOAT_PRECISION> &a) {
					return accum + abs(a) * abs(a);
				})); // make sure to initialize with 0.0 and NOT 0 or it won't be a float and answer will be wrong
		PRISM_FLOAT_PRECISION a = 0;
		for (auto &i : PsiProbeInit) { a += i.real(); };
		transform(PsiProbeInit.begin(), PsiProbeInit.end(),
		          PsiProbeInit.begin(), [&norm_constant](std::complex<PRISM_FLOAT_PRECISION> &a) {
					return a / norm_constant;
				});

		Array3D<complex<PRISM_FLOAT_PRECISION> > trans = zeros_ND<3, complex<PRISM_FLOAT_PRECISION> >(
				{{pars.pot.get_dimk(), pars.pot.get_dimj(), pars.pot.get_dimi()}});
		{
			auto p = pars.pot.begin();
			for (auto &j:trans)j = exp(i * pars.sigma * (*p++));
		}

		pars.stack = zeros_ND<4, PRISM_FLOAT_PRECISION>({{pars.yp.size(), pars.xp.size(), pars.Ndet, 1}}); // TODO: encapsulate stack creation for 3D/4D output

		buildMultisliceOutput(pars, trans, PsiProbeInit, alphaInd);
		//int debug=0;
	}

}
#endif //PRISM_MULTISLICE_H
