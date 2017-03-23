//
// Created by AJ Pryor on 2/13/17.
//

#include "PRISM03.h"
#include "params.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <thread>
#include <mutex>
#include <numeric>
#include <vector>
#include "fftw3.h"
#include "utility.h"
#include "getWorkID.h"

namespace PRISM {
	using namespace std;
	const static std::complex<PRISM_FLOAT_PRECISION> i(0, 1);
	// this might seem a strange way to get pi, but it's slightly more future proof
	const static PRISM_FLOAT_PRECISION pi = std::acos(-1);

	Array2D<PRISM_FLOAT_PRECISION> array2D_subset(const Array2D<PRISM_FLOAT_PRECISION> &arr,
	                                              const size_t &starty, const size_t &stepy, const size_t &stopy,
	                                              const size_t &startx, const size_t &stepx, const size_t &stopx) {
		vector<PRISM_FLOAT_PRECISION> _d;
		size_t dimx = 0;
		size_t dimy = 0;
		for (auto y = starty; y < stopy; y += stepy) {
			for (auto x = startx; x < stopx; x += stepx) {
				_d.push_back(arr.at(y, x));
			}
			++dimy;
		}
		for (auto x = startx; x < stopx; x += stepx)++dimx;
		Array2D<PRISM_FLOAT_PRECISION> result(_d, {{dimy, dimx}});
		return result;
	}

	inline void setupCoordinates(Parameters<PRISM_FLOAT_PRECISION> &pars) {
		Array1D<PRISM_FLOAT_PRECISION> xR = zeros_ND<1, PRISM_FLOAT_PRECISION>({{2}});
		xR[0] = 0.1 * pars.meta.cellDim[2];
		xR[1] = 0.9 * pars.meta.cellDim[2];
		Array1D<PRISM_FLOAT_PRECISION> yR = zeros_ND<1, PRISM_FLOAT_PRECISION>({{2}});
		yR[0] = 0.1 * pars.meta.cellDim[1];
		yR[1] = 0.9 * pars.meta.cellDim[1];

		vector<PRISM_FLOAT_PRECISION> xp_d = vecFromRange(xR[0] + pars.meta.dxy / 2, pars.meta.dxy,
		                                                  xR[1] - pars.meta.dxy / 2);
		vector<PRISM_FLOAT_PRECISION> yp_d = vecFromRange(yR[0] + pars.meta.dxy / 2, pars.meta.dxy,
		                                                  yR[1] - pars.meta.dxy / 2);

		Array1D<PRISM_FLOAT_PRECISION> xp(xp_d, {{xp_d.size()}});
		Array1D<PRISM_FLOAT_PRECISION> yp(yp_d, {{yp_d.size()}});
		pars.xp = xp;
		pars.yp = yp;


	}

	inline void setupDetector(Parameters<PRISM_FLOAT_PRECISION> &pars) {
		// setup coordinates related to detector size, angles, and image output

		pars.dr = 2.5 / 1000;
		pars.alphaMax = pars.qMax * pars.lambda;

		vector<PRISM_FLOAT_PRECISION> detectorAngles_d = vecFromRange(pars.dr / 2, pars.dr,
		                                                              pars.alphaMax - pars.dr / 2);
		Array1D<PRISM_FLOAT_PRECISION> detectorAngles(detectorAngles_d, {{detectorAngles_d.size()}});
		pars.detectorAngles = detectorAngles;
		PRISM_FLOAT_PRECISION r_0 = pars.imageSizeOutput[0] / pars.meta.interpolationFactor / 2;
		PRISM_FLOAT_PRECISION r_1 = pars.imageSizeOutput[1] / pars.meta.interpolationFactor / 2;
		vector<PRISM_FLOAT_PRECISION> yVec_d = vecFromRange(-r_0, (PRISM_FLOAT_PRECISION) 1.0, r_0 - 1);
		vector<PRISM_FLOAT_PRECISION> xVec_d = vecFromRange(-r_1, (PRISM_FLOAT_PRECISION) 1.0, r_1 - 1);
		Array1D<PRISM_FLOAT_PRECISION> yVec(yVec_d, {{yVec_d.size()}});
		Array1D<PRISM_FLOAT_PRECISION> xVec(xVec_d, {{xVec_d.size()}});
		pars.yVec = yVec;
		pars.xVec = xVec;
		pars.Ndet = pars.detectorAngles.size();
	}

	inline void setupBeams(Parameters<PRISM_FLOAT_PRECISION> &pars) {
		// setup some coordinates for the beams

		Array2D<PRISM_FLOAT_PRECISION> beamsReduce = array2D_subset(pars.beamsOutput,
		                                                            0, pars.meta.interpolationFactor,
		                                                            pars.beamsOutput.get_dimj(),
		                                                            0, pars.meta.interpolationFactor,
		                                                            pars.beamsOutput.get_dimi());

		vector<size_t> imageSizeReduce{beamsReduce.get_dimj(), beamsReduce.get_dimi()};
		pars.imageSizeReduce = imageSizeReduce;
		pars.xyBeams = zeros_ND<2, long>({{pars.beamsIndex.size(), 2}});

		for (auto a0 = 1; a0 <= pars.beamsIndex.size(); ++a0) {
			for (auto y = 0; y < beamsReduce.get_dimj(); ++y) {
				for (auto x = 0; x < beamsReduce.get_dimi(); ++x) {
					if (beamsReduce.at(y, x) == a0) {
						pars.xyBeams.at(a0 - 1, 0) = y;
						pars.xyBeams.at(a0 - 1, 1) = x;
					}
				}
			}
		}
	}

	inline void createStack_integrate(Parameters<PRISM_FLOAT_PRECISION> &pars) {
		// create output stack of a size corresponding to 3D mode (integration)

		pars.stack = zeros_ND<4, PRISM_FLOAT_PRECISION>({{pars.yp.size(), pars.xp.size(), pars.Ndet, 1}});
	}

	inline void setupFourierCoordinates(Parameters<PRISM_FLOAT_PRECISION> &pars) {
		// create Fourier space coordinates

		pars.qxaReduce = array2D_subset(pars.qxaOutput,
		                                0, pars.meta.interpolationFactor, pars.qxaOutput.get_dimj(),
		                                0, pars.meta.interpolationFactor, pars.qxaOutput.get_dimi());
		pars.qyaReduce = array2D_subset(pars.qyaOutput,
		                                0, pars.meta.interpolationFactor, pars.qyaOutput.get_dimj(),
		                                0, pars.meta.interpolationFactor, pars.qyaOutput.get_dimi());
		pars.q1 = zeros_ND<2, PRISM_FLOAT_PRECISION>({{pars.imageSizeReduce[0], pars.imageSizeReduce[1]}});
		pars.q2 = zeros_ND<2, PRISM_FLOAT_PRECISION>({{pars.imageSizeReduce[0], pars.imageSizeReduce[1]}});
	}


	void buildPRISMOutput_CPUOnly(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                              const PRISM_FLOAT_PRECISION xTiltShift,
	                              const PRISM_FLOAT_PRECISION yTiltShift) {

		// launch threads to compute results for batches of xp, yp
		// I do this by dividing the xp points among threads, and each computes
		// all of the relevant yp for each of its xp. This seems an okay strategy
		// as long as the number of xp and yp are similar. If that is not the case
		// this may need to be adapted
		vector<thread> workers;
		workers.reserve(pars.meta.NUM_THREADS); // prevents multiple reallocations
		setWorkStartStop(0, pars.xp.size() * pars.yp.size(), 1);
		for (auto t = 0; t < pars.meta.NUM_THREADS; ++t) {
			cout << "Launching CPU worker thread #" << t << " to compute partial PRISM result\n";
			// emplace_back is better whenever constructing a new object
			workers.emplace_back(thread([&pars, &xTiltShift, &yTiltShift]() {
				size_t Nstart, Nstop, ay, ax;
				while (getWorkID(pars, Nstart, Nstop)) { // synchronously get work assignment
					while (Nstart != Nstop) {
						ay = Nstart / pars.xp.size();
						ax = Nstart % pars.xp.size();
						buildSignal_CPU(pars, ay, ax, yTiltShift, xTiltShift);
						++Nstart;
					}
				}
			}));
		}
		// synchronize
		cout << "Waiting for threads...\n";
		for (auto &t:workers)t.join();
	}


	void buildSignal_CPU(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                     const size_t &ay,
	                     const size_t &ax,
	                     const PRISM_FLOAT_PRECISION &yTiltShift,
	                     const PRISM_FLOAT_PRECISION &xTiltShift) {
		static mutex fftw_plan_lock; // for synchronizing access to shared FFTW resources
		const static std::complex<PRISM_FLOAT_PRECISION> i(0, 1);
		const static PRISM_FLOAT_PRECISION pi = std::acos(-1);

		PRISM_FLOAT_PRECISION x0 = pars.xp[ax] / pars.pixelSizeOutput[1];
		PRISM_FLOAT_PRECISION y0 = pars.yp[ay] / pars.pixelSizeOutput[0];
		Array1D<PRISM_FLOAT_PRECISION> x = pars.xVec + round(x0);
		transform(x.begin(), x.end(), x.begin(), [&pars](PRISM_FLOAT_PRECISION &a) {
			return fmod(a, (PRISM_FLOAT_PRECISION) pars.imageSizeOutput[1]);
		});
		Array1D<PRISM_FLOAT_PRECISION> y = pars.yVec + round(y0);
		transform(y.begin(), y.end(), y.begin(), [&pars](PRISM_FLOAT_PRECISION &a) {
			return fmod(a, (PRISM_FLOAT_PRECISION) pars.imageSizeOutput[0]);
		});
		Array2D<PRISM_FLOAT_PRECISION> intOutput = PRISM::zeros_ND<2, PRISM_FLOAT_PRECISION>(
				{{pars.imageSizeReduce[0], pars.imageSizeReduce[1]}});
		for (auto a5 = 0; a5 < pars.meta.numFP; ++a5) {
			Array2D<std::complex<PRISM_FLOAT_PRECISION> > psi = PRISM::zeros_ND<2, std::complex<PRISM_FLOAT_PRECISION> >(
					{{pars.imageSizeReduce[0], pars.imageSizeReduce[1]}});
			for (auto a4 = 0; a4 < pars.beamsIndex.size(); ++a4) {
				PRISM_FLOAT_PRECISION yB = pars.xyBeams.at(a4, 0);
				PRISM_FLOAT_PRECISION xB = pars.xyBeams.at(a4, 1);
				if (abs(pars.psiProbeInit.at(yB, xB)) > 0) {
					PRISM_FLOAT_PRECISION q0_0 = pars.qxaReduce.at(yB, xB);
					PRISM_FLOAT_PRECISION q0_1 = pars.qyaReduce.at(yB, xB);
					std::complex<PRISM_FLOAT_PRECISION> phaseShift = exp(
							-2 * pi * i * (q0_0 * (pars.xp[ax] + xTiltShift) +
							               q0_1 * (pars.yp[ay] + yTiltShift)));
//                    std::complex<PRISM_FLOAT_PRECISION> phaseShift = exp(-2 * pi * i);
					// caching this constant made a 5x performance improvement even with
					// full compiler optimization turned on. Optimizing compilers aren't perfect...
					const std::complex<PRISM_FLOAT_PRECISION> tmp_const = pars.psiProbeInit.at(yB, xB) * phaseShift;
					auto psi_ptr = psi.begin();
					for (auto j = 0; j < y.size(); ++j) {
						for (auto i = 0; i < x.size(); ++i) {
							// access contiguously for performance
							*psi_ptr++ += (tmp_const * pars.Scompact.at(a4, y[j], x[i]));
						}
					}
				}
			}

			// fftw_execute is the only thread-safe function in the library, so we need to synchronize access
			// to the plan creation methods
			unique_lock<mutex> gatekeeper(fftw_plan_lock);
			PRISM_FFTW_PLAN plan = PRISM_FFTW_PLAN_DFT_2D(psi.get_dimj(), psi.get_dimi(),
			                                              reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi[0]),
			                                              reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi[0]),
			                                              FFTW_FORWARD, FFTW_ESTIMATE);

			gatekeeper.unlock(); // unlock it so we only block as long as necessary to deal with plans
			PRISM_FFTW_EXECUTE(plan);
			gatekeeper.lock();
			PRISM_FFTW_DESTROY_PLAN(plan);
			gatekeeper.unlock();

			for (auto jj = 0; jj < intOutput.get_dimj(); ++jj) {
				for (auto ii = 0; ii < intOutput.get_dimi(); ++ii) {
					intOutput.at(jj, ii) += pow(abs(psi.at(jj, ii)), 2);
//	                intOutput.at(ii,jj) += pow(abs(psi.at(jj,ii)),2);
				}
			}
		}

//         update stack -- ax,ay are unique per thread so this write is thread-safe without a lock
		auto idx = pars.alphaInd.begin();
		for (auto counts = intOutput.begin(); counts != intOutput.end(); ++counts) {
			if (*idx <= pars.Ndet) {
				pars.stack.at(ay, ax, (*idx) - 1, 0) += *counts * pars.scale;
			}
			++idx;
		};
	}

	void PRISM03(Parameters<PRISM_FLOAT_PRECISION> &pars) {
		// compute final image

		cout << "Entering PRISM03" << endl;



		// setup necessary coordinates
		setupCoordinates(pars);

		// setup angles of detector and image sizes
		setupDetector(pars);

		// setup coordinates and indices for the beams
		setupBeams(pars);

		// setup Fourier coordinates for the S-matrix
		setupFourierCoordinates(pars);

		// initialize the output stack to the correct size for the output mode
		createStack_integrate(pars);

		Array2D<PRISM_FLOAT_PRECISION> intOutput = zeros_ND<2, PRISM_FLOAT_PRECISION>(
				{{pars.imageSizeReduce[0], pars.imageSizeReduce[1]}});
		Array2D<complex<PRISM_FLOAT_PRECISION> > psi;
		psi = zeros_ND<2, complex<PRISM_FLOAT_PRECISION> >({{pars.imageSizeReduce[0], pars.imageSizeReduce[1]}});

		pars.psiProbeInit = zeros_ND<2, complex<PRISM_FLOAT_PRECISION> >(
				{{pars.imageSizeReduce[0], pars.imageSizeReduce[1]}});
		pars.dq = (pars.qxaReduce.at(0, 1) + pars.qyaReduce.at(1, 0)) / 2;
		PRISM_FLOAT_PRECISION scale = pow(pars.meta.interpolationFactor, 4);
		pars.scale = scale;

		// The operators +, -, /, * return PRISM arrays by value, so to avoid unnecessary memory
		// allocations/copies for chained operations I try to do things like create variables
		// initially with at most one operation, and then perform in-place transforms if more is needed
		PRISM_FLOAT_PRECISION qProbeMax = pars.probeSemiangle / pars.lambda;
		Array2D<PRISM_FLOAT_PRECISION> qxaShift = pars.qxaReduce - (pars.probeXtilt / pars.lambda);
		Array2D<PRISM_FLOAT_PRECISION> qyaShift = pars.qyaReduce - (pars.probeYtilt / pars.lambda);
		transform(qxaShift.begin(), qxaShift.end(),
		          qyaShift.begin(), pars.q2.begin(),
		          [](const PRISM_FLOAT_PRECISION &a, const PRISM_FLOAT_PRECISION &b) { return a * a + b * b; });
		transform(pars.q2.begin(), pars.q2.end(),
		          pars.q1.begin(),
		          [](const PRISM_FLOAT_PRECISION &a) { return sqrt(a); });
//						Array2D<PRISM_FLOAT_PRECISION> alphaInd(pars.q1); // copy constructor more efficient than assignment
		pars.alphaInd = pars.q1;
		transform(pars.alphaInd.begin(), pars.alphaInd.end(),
		          pars.alphaInd.begin(),
		          [&pars](const PRISM_FLOAT_PRECISION &a) {
			          return 1 + round((a * pars.lambda - pars.detectorAngles[0]) / pars.dr);
		          });
		transform(pars.alphaInd.begin(), pars.alphaInd.end(),
		          pars.alphaInd.begin(),
		          [](const PRISM_FLOAT_PRECISION &a) { return a < 1 ? 1 : a; });
		PRISM::ArrayND<2, std::vector<unsigned short> > alphaMask(
				std::vector<unsigned short>(pars.alphaInd.size(), 0),
				{{pars.alphaInd.get_dimj(), pars.alphaInd.get_dimi()}});
		transform(pars.alphaInd.begin(), pars.alphaInd.end(),
		          alphaMask.begin(),
		          [&pars](const PRISM_FLOAT_PRECISION &a) { return (a < pars.Ndet) ? 1 : 0; });

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
			          a = a * exp(-i * pi * pars.lambda * pars.probeDefocus * q2_t);
			          return a;
		          });
		PRISM_FLOAT_PRECISION norm_constant = sqrt(accumulate(pars.psiProbeInit.begin(), pars.psiProbeInit.end(),
		                                                      (PRISM_FLOAT_PRECISION) 0.0,
		                                                      [](PRISM_FLOAT_PRECISION accum,
		                                                         std::complex<PRISM_FLOAT_PRECISION> &a) {
			                                                      return accum + abs(a) * abs(a);
		                                                      })); // make sure to initialize with 0.0 and NOT 0 or it won't be a float and answer will be wrong
		PRISM_FLOAT_PRECISION a = 0;
		for (auto &i : pars.psiProbeInit) { a += i.real(); };
		transform(pars.psiProbeInit.begin(), pars.psiProbeInit.end(),
		          pars.psiProbeInit.begin(), [&norm_constant](std::complex<PRISM_FLOAT_PRECISION> &a) {
					return a / norm_constant;
				});


//		PRISM_FLOAT_PRECISION zTotal = pars.meta.cellDim[0];
//		PRISM_FLOAT_PRECISION xTiltShift = -zTotal * tan(pars.probeXtilt);
//		PRISM_FLOAT_PRECISION yTiltShift = -zTotal * tan(pars.probeYtilt);

		buildPRISMOutput(pars, xTiltShift, yTiltShift);
	}
}
