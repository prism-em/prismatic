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

	}

	void buildSignal_GPU(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                     const size_t &ay,
	                     const size_t &ax,
	                     const PRISM_FLOAT_PRECISION &yTiltShift,
	                     const PRISM_FLOAT_PRECISION &xTiltShift,
	                     const Array2D<PRISM_FLOAT_PRECISION> &alphaInd,
	                     const Array2D<std::complex<PRISM_FLOAT_PRECISION> > &PsiProbeInit) {
		cout << "DUMMY CODE" << endl;
	}
}