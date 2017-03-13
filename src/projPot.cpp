//
// Created by AJ Pryor on 2/22/17.
//

#include "projPot.h"
#include <vector>
#include "ArrayND.h"
#include "boost/math/special_functions/bessel.hpp"
#include <math.h>
#include <iostream>
#include <algorithm>
#include <numeric>
#include "fparams.h"
#include "defines.h"
#include "configure.h"

namespace PRISM {

	PRISM_FLOAT_PRECISION get_potMin(const Array2D<PRISM_FLOAT_PRECISION>& pot,
	                                        const Array1D<PRISM_FLOAT_PRECISION>& xr,
	                                        const Array1D<PRISM_FLOAT_PRECISION>& yr){
		// I am assuming that xr and yr are symmetric about 0
		const size_t xInd = std::floor(xr.size()/2);
		const size_t yInd = std::floor(yr.size()/2);
		const PRISM_FLOAT_PRECISION dx        = round(sqrt(2*(xInd + 1)- 1));
		const PRISM_FLOAT_PRECISION dy        = round(sqrt(2*(yInd + 1) - 1));
		const PRISM_FLOAT_PRECISION xv[]      = {xInd-dx, xInd+dx, xInd-dx, xInd+dx, 0, 0, (PRISM_FLOAT_PRECISION)xr.size()-1, (PRISM_FLOAT_PRECISION)xr.size()-1};
//		const PRISM_FLOAT_PRECISION xv[]      = {xInd-dx-1, xInd+dx+1, xInd-dx-1, xInd+dx+1, 0, 0, (PRISM_FLOAT_PRECISION)xr.size()-1, (PRISM_FLOAT_PRECISION)xr.size()-1};
		const PRISM_FLOAT_PRECISION yv[]      = {0, 0, (PRISM_FLOAT_PRECISION)yr.size()-1, (PRISM_FLOAT_PRECISION)yr.size()-1, yInd-dy, yInd+dy, yInd-dy, yInd+dy};

		PRISM_FLOAT_PRECISION potMin = 0;
		for (auto i=0; i < 8; ++i)potMin = (pot.at(yv[i],xv[i]) > potMin) ? pot.at(yv[i],xv[i]) : potMin;
		return potMin;
	}

	using namespace std;

	Array2D<PRISM_FLOAT_PRECISION> projPot(const size_t &Z,
	                                              const Array1D<PRISM_FLOAT_PRECISION> &xr,
	                                              const Array1D<PRISM_FLOAT_PRECISION> &yr) {
		static const PRISM_FLOAT_PRECISION pi = std::acos(-1);
		PRISM_FLOAT_PRECISION ss    = 8;
		PRISM_FLOAT_PRECISION a0    = 0.5292;
		PRISM_FLOAT_PRECISION e     = 14.4;
		PRISM_FLOAT_PRECISION term1 = 4*pi*pi*a0*e;
		PRISM_FLOAT_PRECISION term2 = 2*pi*pi*a0*e;
		ArrayND<2, std::vector<PRISM_FLOAT_PRECISION> > result = zeros_ND<2, PRISM_FLOAT_PRECISION>({{yr.size(), xr.size()}});
		const PRISM_FLOAT_PRECISION dx = xr[1] - xr[0];
		const PRISM_FLOAT_PRECISION dy = yr[1] - yr[0];

		PRISM_FLOAT_PRECISION start = -(ss-1)/ss/2;
		const PRISM_FLOAT_PRECISION step  = 1/ss;
		const PRISM_FLOAT_PRECISION end   = -start;
		vector<PRISM_FLOAT_PRECISION> sub_data;
		while (start <= end){
			sub_data.push_back(start);
			start+=step;
		}
		ArrayND<1, std::vector<PRISM_FLOAT_PRECISION> > sub(sub_data,{{sub_data.size()}});

		std::pair<Array2D<PRISM_FLOAT_PRECISION>, Array2D<PRISM_FLOAT_PRECISION> > meshx = meshgrid(xr, sub*dx);
		std::pair<Array2D<PRISM_FLOAT_PRECISION>, Array2D<PRISM_FLOAT_PRECISION> > meshy = meshgrid(yr, sub*dy);

		ArrayND<1, std::vector<PRISM_FLOAT_PRECISION> > xv = zeros_ND<1, PRISM_FLOAT_PRECISION>({{meshx.first.size()}});
		ArrayND<1, std::vector<PRISM_FLOAT_PRECISION> > yv = zeros_ND<1, PRISM_FLOAT_PRECISION>({{meshy.first.size()}});
		{
			auto t_x = xv.begin();
			for (auto i = 0; i < meshx.first.get_dimi(); ++i) {
				for (auto j = 0; j < meshx.first.get_dimj(); ++j) {
					*t_x++ = meshx.first.at(j, i) + meshx.second.at(j, i);
				}
			}
		}

		{
			auto t_y = yv.begin();
			for (auto i = 0; i < meshy.first.get_dimi(); ++i) {
				for (auto j = 0; j < meshy.first.get_dimj(); ++j) {
					*t_y++ = meshy.first.at(j, i) + meshy.second.at(j, i);
				}
			}
		}

		std::pair<Array2D<PRISM_FLOAT_PRECISION>, Array2D<PRISM_FLOAT_PRECISION> > meshxy = meshgrid(xv, yv);
		ArrayND<2, std::vector<PRISM_FLOAT_PRECISION> > r2 = zeros_ND<2, PRISM_FLOAT_PRECISION>({{yv.size(), xv.size()}});
		ArrayND<2, std::vector<PRISM_FLOAT_PRECISION> > r  = zeros_ND<2, PRISM_FLOAT_PRECISION>({{yv.size(), xv.size()}});

		{
			auto t_y = r2.begin();
			for (auto i = 0; i < meshxy.first.get_dimi(); ++i) {
				for (auto j = 0; j < meshxy.first.get_dimj(); ++j) {
					*t_y++ = pow(meshxy.first.at(j,i),2) + pow(meshxy.second.at(j,i),2);
				}
			}
		}


		for (auto i = 0; i < r.size(); ++i)r[i] = sqrt(r2[i]);
		// construct potential
		ArrayND<2, std::vector<PRISM_FLOAT_PRECISION> > potSS  = ones_ND<2, PRISM_FLOAT_PRECISION>({{r2.get_dimj(), r2.get_dimi()}});
		std::vector<PRISM_FLOAT_PRECISION> ap;
		ap.reserve(n_parameters);
		for (auto i = 0; i < n_parameters; ++i){
			ap[i] = fparams[(Z-1)*n_parameters + i];
		}

		using namespace boost::math;
		std::transform(r.begin(), r.end(),
		               r2.begin(), potSS.begin(), [&ap, &term1, &term2](const PRISM_FLOAT_PRECISION& r_t, const PRISM_FLOAT_PRECISION& r2_t){

					return term1*(ap[0] *
					              cyl_bessel_k(0,2*pi*sqrt(ap[1])*r_t)          +
					              ap[2]*cyl_bessel_k(0,2*pi*sqrt(ap[3])*r_t)    +
					              ap[4]*cyl_bessel_k(0,2*pi*sqrt(ap[5])*r_t))   +
					       term2*(ap[6]/ap[7]*exp(-pow(pi,2)/ap[7]*r2_t) +
					              ap[8]/ap[9]*exp(-pow(pi,2)/ap[9]*r2_t)        +
					              ap[10]/ap[11]*exp(-pow(pi,2)/ap[11]*r2_t));
				});

		// integrate
		ArrayND<2, std::vector<PRISM_FLOAT_PRECISION> > pot = zeros_ND<2, PRISM_FLOAT_PRECISION>({{yr.size(), xr.size()}});
		for (auto sy = 0; sy < ss; ++sy){
			for (auto sx = 0; sx < ss; ++sx) {
				for (auto j = 0; j < pot.get_dimj(); ++j) {
					for (auto i = 0; i < pot.get_dimi(); ++i) {
						pot.at(j, i) += potSS.at(j*ss + sy, i*ss + sx);
					}
				}
			}
		}
		pot/=(ss*ss);

		PRISM_FLOAT_PRECISION potMin = get_potMin(pot,xr,yr);
		pot -= potMin;
		transform(pot.begin(),pot.end(),pot.begin(),[](PRISM_FLOAT_PRECISION& a){return a<0?0:a;});

		return pot;
	}
}
