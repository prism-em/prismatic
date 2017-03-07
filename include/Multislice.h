//
// Created by AJ Pryor on 3/6/17.
//

#ifndef PRISM_MULTISLICE_H
#define PRISM_MULTISLICE_H
#include <iostream>
#include <vector>
#include "configure.h"
#include "meta.h"
#include "ArrayND.h"
#include "params.h"
#include "utility.h"
namespace PRISM{
	template <class T>
	void buildMultisliceOutput_cpuOnly(Parameters<T>& pars){
		using namespace std;
		cout << " test " << endl;

	};
	template <class T>
	void Multislice(Parameters<T>& pars){
		using namespace std;
		const double pi = acos(-1);
		const std::complex<double> i(0, 1);
		const T dxy = 0.25 * 2;

		// should move these elsewhere and in PRISM03
		pars.probeDefocusArray = zeros_ND<1, T>({{1}});
		pars.probeSemiangleArray = zeros_ND<1, T>({{1}});
		pars.probeXtiltArray = zeros_ND<1, T>({{1}});
		pars.probeYtiltArray = zeros_ND<1, T>({{1}});
		pars.probeDefocusArray[0] = 0.0;
		pars.probeSemiangleArray[0] = 20.0 / 1000;
		pars.probeXtiltArray[0] = 0.0 / 1000;
		pars.probeYtiltArray[0] = 0.0 / 1000;

		Array1D<T> xR = zeros_ND<1, T>({{2}});
		xR[0] = 0.1 * pars.meta.cellDim[2];
		xR[1] = 0.9 * pars.meta.cellDim[2];
		Array1D<T> yR = zeros_ND<1, T>({{2}});
		yR[0] = 0.1 * pars.meta.cellDim[1];
		yR[1] = 0.9 * pars.meta.cellDim[1];
		vector<T> xp_d = vecFromRange(xR[0] + dxy / 2, dxy, xR[1] - dxy / 2);
		vector<T> yp_d = vecFromRange(yR[0] + dxy / 2, dxy, yR[1] - dxy / 2);

		Array1D<T> xp(xp_d, {{xp_d.size()}});
		Array1D<T> yp(yp_d, {{yp_d.size()}});
		pars.xp = xp;
		pars.yp = yp;

		// setup some coordinates
		pars.imageSize[0] = pars.pot.get_dimj();
		pars.imageSize[1] = pars.pot.get_dimi();
		Array1D<T> qx = makeFourierCoords(pars.imageSize[1], pars.pixelSize[1]);
		Array1D<T> qy = makeFourierCoords(pars.imageSize[0], pars.pixelSize[0]);

		pair< Array2D<T>, Array2D<T> > mesh = meshgrid(qx,qy);
		pars.qxa = mesh.first;
		pars.qya = mesh.second;
		Array2D<T> q2(pars.qya);
		transform(pars.qxa.begin(), pars.qxa.end(),
		          pars.qya.begin(), q2.begin(), [](const T& a, const T& b){
					return a*a + b*b;
				});
		Array2D<T> q1(q2);
		for (auto& q : q1)q=sqrt(q);

		// get qMax
		pars.qMax = 0;
		{
			T qx_max;
			T qy_max;
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
		pars.prop     = zeros_ND<2, std::complex<T> >({{pars.imageSize[1], pars.imageSize[0]}});
		pars.propBack = zeros_ND<2, std::complex<T> >({{pars.imageSize[1], pars.imageSize[0]}});
		for (auto y = 0; y < pars.qMask.get_dimj(); ++y) {
			for (auto x = 0; x < pars.qMask.get_dimi(); ++x) {
				if (pars.qMask.at(y,x)==1)
				{
					pars.prop.at(y,x)     = exp(-i * pi * complex<T>(pars.lambda, 0) *
					                            complex<T>(pars.meta.sliceThickness, 0) *
					                            complex<T>(q2.at(y, x), 0));
					pars.propBack.at(y,x) = exp(i * pi * complex<T>(pars.lambda, 0) *
					                            complex<T>(pars.meta.cellDim[0], 0) *
					                            complex<T>(q2.at(y, x), 0));
				}
			}
		}

		pars.dr = 2.5 / 1000;
		pars.alphaMax = pars.qMax * pars.lambda;
		vector<T> detectorAngles_d = vecFromRange(pars.dr / 2, pars.dr, pars.alphaMax - pars.dr / 2);
		Array1D<T> detectorAngles(detectorAngles_d, {{detectorAngles_d.size()}});
		pars.detectorAngles = detectorAngles;
		Array2D<T> alpha = q1 * pars.lambda;
//		Array2D<T> alphaInds = (alpha + pars.dr/2) / pars.dr;
//		for (auto& q : alphaInds) q = round(q);

		pars.Ndet = pars.detectorAngles.size();
		if (pars.probeSemiangleArray.size() > 1)throw std::domain_error("Currently only scalar probeSemiangleArray supported. Multiple inputs received.\n");
		T qProbeMax = pars.probeSemiangleArray[0] / pars.lambda; // currently a single semiangle
		Array2D<complex<T> > PsiProbeInit, psi;
		PsiProbeInit = zeros_ND<2, complex<T> >({{q1.get_dimj(), q1.get_dimi()}});
		psi = zeros_ND<2, complex<T> >({{q1.get_dimj(), q1.get_dimi()}});
		pars.dq = (pars.qxa.at(0, 1) + pars.qya.at(1, 0)) / 2;

		transform(PsiProbeInit.begin(), PsiProbeInit.end(),
		          q1.begin(), PsiProbeInit.begin(),
		          [&pars, &qProbeMax](std::complex<T> &a, T &q1_t) {
			          a.real(erf((qProbeMax - q1_t) / (0.5 * pars.dq)) * 0.5 + 0.5);
			          a.imag(0);
			          return a;
		          });


		transform(PsiProbeInit.begin(), PsiProbeInit.end(),
		          q2.begin(), PsiProbeInit.begin(),
		          [&pars, &i, &pi](std::complex<T> &a, T &q2_t) {
			          a = a * exp(-i * pi * pars.lambda * pars.probeDefocusArray[0] * q2_t); // TODO: fix hardcoded length-1 defocus
			          return a;
		          });
		T norm_constant = sqrt(accumulate(PsiProbeInit.begin(), PsiProbeInit.end(),
		                                  0.0, [](T accum, std::complex<T> &a) {
					return accum + abs(a) * abs(a);
				})); // make sure to initialize with 0.0 and NOT 0 or it won't be a float and answer will be wrong
		double a = 0;
		for (auto &i : PsiProbeInit) { a += i.real(); };
		transform(PsiProbeInit.begin(), PsiProbeInit.end(),
		          PsiProbeInit.begin(), [&norm_constant](std::complex<T> &a) {
					return a / norm_constant;
				});

		Array3D<complex<T> > trans = zeros_ND<3, complex<T> >(
				{{pars.pot.get_dimk(), pars.pot.get_dimj(), pars.pot.get_dimi()}});
		{
			auto p = pars.pot.begin();
			for (auto &j:trans)j = exp(i * pars.sigma * (*p++));
		}

		pars.stack = zeros_ND<4, T>({{pars.yp.size(), pars.xp.size(), pars.Ndet, 1}}); // TODO: encapsulate stack creation for 3D/4D output

//		buildMultisliceOutput(pars);
		int debug=0;
	}

}
#endif //PRISM_MULTISLICE_H
