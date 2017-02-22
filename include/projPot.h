//
// Created by AJ Pryor on 2/22/17.
//

#ifndef PRISM_PROJPOT_H
#define PRISM_PROJPOT_H
#include "ArrayND.h"
#include "boost/math/special_functions/bessel.hpp"
#include <math.h>
#include <iostream>
#include <algorithm>
#include <numeric>
#include "fparams.h"
namespace PRISM {
	using namespace std;
	template <class T>
	using Array2D = PRISM::ArrayND<2, std::vector<T> >;
	template <class T>
	using Array1D = PRISM::ArrayND<1, std::vector<T> >;

	template<class T>
	ArrayND<2, std::vector<T> >
	projPot(const size_t &Z, const ArrayND<1, std::vector<T> > &xr, const ArrayND<1, std::vector<T> > &yr) {
		static const T pi = std::acos(-1);
		T ss = 8;
		T a0 = 0.5292;
		T e = 14.4;
		T term1 = 4*pi*pi*a0*e;
		T term2 = 2*pi*pi*a0*e;
		ArrayND<2, std::vector<T> > result = zeros_ND<2, T>({yr.size(), xr.size()});
		const T dx = xr[1] - xr[0];
		const T dy = yr[1] - yr[0];

		T start = -(ss-1)/ss/2;
		const T step  = 1/ss;
		const T end   = -start;
		vector<T> sub_data;
		while (start <= end){
			sub_data.push_back(start);
			start+=step;
		}
		ArrayND<1, std::vector<T> > sub(sub_data,{sub_data.size()});
		for (auto& i : sub)cout << i << endl;

		std::pair<Array2D<T>, Array2D<T> > meshx = meshgrid(xr, sub*dx);
		std::pair<Array2D<T>, Array2D<T> > meshy = meshgrid(yr, sub*dy);

//		for (auto i = 0; i < meshx.first.get_ncols();++i){
//			for (auto j = 0; j < meshx.first.get_nrows(); ++j){
//				cout <<" j,i = " << j << ", " << i << endl;
//				cout << meshx.first.at(j,i) << endl;
//			}
//		}
//
//		for (auto i = 0; i < meshx.first.get_ncols();++i){
//			for (auto j = 0; j < meshx.first.get_nrows(); ++j){
//				cout <<" j,i = " << j << ", " << i << endl;
//				cout << meshx.second.at(j,i) << endl;
//			}
//		}

		ArrayND<1, std::vector<T> > xv = zeros_ND<1, T>({meshx.first.size()});
		ArrayND<1, std::vector<T> > yv = zeros_ND<1, T>({meshy.first.size()});
		{
			auto t_x = xv.begin();
			for (auto i = 0; i < meshx.first.get_ncols(); ++i) {
				for (auto j = 0; j < meshx.first.get_nrows(); ++j) {
					*t_x++ = meshx.first.at(j, i) + meshx.second.at(j, i);
				}
			}
		}

//		cout <<"print xv\n";
//		for (auto &i : xv)cout<<i << endl;



		{
			auto t_y = yv.begin();
			for (auto i = 0; i < meshy.first.get_ncols(); ++i) {
				for (auto j = 0; j < meshy.first.get_nrows(); ++j) {
					*t_y++ = meshy.first.at(j, i) + meshy.second.at(j, i);
				}
			}
		}

//		cout <<"print yv\n";
//		for (auto &i : yv)cout<<i << endl;

//		cout << "yv.at[0] = " << yv[0] << endl;
//		cout << "yv.at[5] = " << yv[5] << endl;
//		cout << "yv.at[10] = " << yv[10] << endl;
//		cout << "yv.at[100] = " << yv[100] << endl;
//		cout << "yv.at[0] = " << yv[0] << endl;



		std::pair<Array2D<T>, Array2D<T> > meshxy = meshgrid(xv, yv);
		ArrayND<2, std::vector<T> > r2 = zeros_ND<2, T>({yv.size(), xv.size()});
		ArrayND<2, std::vector<T> > r = zeros_ND<2, T>({yv.size(), xv.size()});
//		cout << "r.get_ncols = " << r.get_ncols()<<endl;
//		cout << "r.get_nrows = " << r.get_nrows()<<endl;
//		cout << "meshxy.first.size() = " << meshxy.first.size() << endl;
//		cout << "meshxy.first.get_nrows() = " << meshxy.first.get_nrows() << endl;
//		cout << "meshxy.first.get_ncols() = " << meshxy.first.get_ncols() << endl;


//		for (auto j = 0; j < meshxy.first.size(); ++j)r2[j]=pow(meshxy.first[j],2) + pow(meshxy.second[j],2);
		{
			auto t_y = r2.begin();
			for (auto i = 0; i < meshxy.first.get_nrows(); ++i) {
				for (auto j = 0; j < meshxy.first.get_ncols(); ++j) {
					*t_y++ = pow(meshxy.first.at(j,i),2) + pow(meshxy.second.at(j,i),2);
//					r2.at(j,i) = pow(meshxy.first.at(j,i),2) + pow(meshxy.second.at(j,i),2);
				}
			}
		}

//		cout << "r2.at(1,1) = " << r2.at(1,1) << endl;
//		cout << "r2.at(10,10) = " << r2.at(10,10) << endl;
//		cout << "r2.at(100,100) = " << r2.at(100,100) << endl;
//		cout << "r2.at(96,96) = " << r2.at(96,96) << endl;

		for (auto i = 0; i < r.size(); ++i)r[i] = sqrt(r2[i]);
		// construct potential
		ArrayND<2, std::vector<T> > potSS  = ones_ND<2, T>({r2.get_nrows(), r2.get_ncols()});
		std::vector<double> ap;
		ap.reserve(n_parameters);
		for (auto i = 0; i < n_parameters; ++i){
			ap[i] = fparams[(Z-1)*n_parameters + i];
		}

//		for (auto &i:r)cout<<i<<endl;
		using namespace boost::math;
		cout << "2*pi*sqrt(ap[1])*r[0] = " << 2*pi*sqrt(ap[1])*r[0] << endl;
		cout << "cyl_bessel_k(2*pi*sqrt(ap[1])*r[0]) = " << cyl_bessel_k(0,2*pi*sqrt(ap[1])*r[0]) << endl;

		cout << "r2[0] = " << r2[0] << endl;
		cout << "r2[1] = " << r2[1] << endl;
		cout << "r2[2] = " << r2[2] << endl;
		cout << "r[0] = " << r[0] << endl;
		cout << "r[1] = " << r[1] << endl;
		cout << "r[28223] = " << r[28223] << endl;
		cout << "pi = " << pi << endl;
//		potSS = term1 * ( ap[0]*cyl_bessel_k(0,2*pi*sqrt(ap[1]))*r);
		std::transform(r.begin() , r.end() ,
		          r2.begin(),
		          potSS.begin(),[&ap, &term1, &term2](const T& r_t, const T& r2_t){

			return term1*(ap[0] *
	                cyl_bessel_k(0,2*pi*sqrt(ap[1])*r_t)      +
					ap[2]*cyl_bessel_k(0,2*pi*sqrt(ap[3])*r_t)  +
					ap[4]*cyl_bessel_k(0,2*pi*sqrt(ap[5])*r_t)) +
					term2*(ap[6]/ap[7]*exp(-pow(pi,2)/ap[7]*r2_t)    +
					ap[8]/ap[9]*exp(-pow(pi,2)/ap[9]*r2_t)          +
					ap[10]/ap[11]*exp(-pow(pi,2)/ap[11]*r2_t));

		});

		cout << "potSS[0] = " << potSS[0] << endl;
		cout << "potSS[2] = " << potSS[2] << endl;

		// integrate
		ArrayND<2, std::vector<T> > potMid = zeros_ND<2, T>({yr.size(), xr.size()});


#ifndef NDEBUG
#include <iostream>
		using namespace std;
		cout << meshx.first.get_nrows() << endl;
		cout << meshx.second.get_nrows() << endl;
		cout << meshx.first.size() << endl;
		cout << "meshx.first.at(0,1) = " << meshx.first.at(0,1)<< endl;
		cout << "meshx.first.at(0,2) = " << meshx.first.at(0,2)<< endl;
		cout << "meshx.first.at(1,1) = " << meshx.first.at(1,1)<< endl;
		cout << "meshx.second.at(0,1) = " << meshx.second.at(0,1)<< endl;
		cout << "meshx.second.at(0,2) = " << meshx.second.at(0,2)<< endl;
		cout << "meshx.second.at(1,1) = " << meshx.second.at(1,1)<< endl;
		cout << "meshy.first.at(0,1) = " << meshy.first.at(0,1)<< endl;
		cout << "meshy.first.at(0,2) = " << meshy.first.at(0,2)<< endl;
		cout << "meshy.first.at(1,1) = " << meshy.first.at(1,1)<< endl;
		cout << "meshy.second.at(0,1) = " << meshy.second.at(0,1)<< endl;
		cout << "meshy.second.at(0,2) = " << meshy.second.at(0,2)<< endl;
		cout << "meshy.second.at(1,1) = " << meshy.second.at(1,1)<< endl;
		cout << "xv[0] = " << xv[0]<< endl;
		cout << "xv[2] = " << xv[2]<< endl;
		cout << "xv[4] = " << xv[4]<< endl;
		cout << "yv[0] = " << xv[0]<< endl;
		cout << "yv[2] = " << xv[2]<< endl;
		cout << "yv[4] = " << xv[4]<< endl;
//		cout << "meshxy.second.at(0,1) = " << meshxy.second.at(0,1)<< endl;
//		cout << "meshxy.second.at(0,2) = " << meshxy.second.at(0,2)<< endl;
//		cout << "meshxy.second.at(1,1) = " << meshxy.second.at(1,1)<< endl;
		cout << "r2.at(1,1) = " << r2.at(1,1)<< endl;
		cout << "r2.at(3,2) = " << r2.at(3,2)<< endl;
		cout << "a0 = " << a0 << endl;
		cout << "term1 = " << term1 << endl;
		cout << "term2 = " << term2 << endl;
		cout << "term1 = " << term1 << endl;
		cout << "ap terms\n";
		for (auto &i : ap)cout << i << endl;

#endif //NDEBUG
		return result;
	}
}
#endif //PRISM_PROJPOT_H
