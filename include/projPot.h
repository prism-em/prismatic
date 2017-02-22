//
// Created by AJ Pryor on 2/22/17.
//

#ifndef PRISM_PROJPOT_H
#define PRISM_PROJPOT_H
#include "ArrayND.h"
#include <math.h>
#include <iostream>
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
		ArrayND<1, std::vector<T> > xv = zeros_ND<1, T>({meshx.first.size()});
		ArrayND<1, std::vector<T> > yv = zeros_ND<1, T>({meshy.first.size()});
		{
			auto t_x = xv.begin();
			for (auto i = 0; i < meshx.first.get_nrows(); ++i) {
				for (auto j = 0; j < meshx.first.get_ncols(); ++j) {
					*t_x++ = meshx.first.at(j, i) + meshx.second.at(j, i);
				}
			}
		}
		{
			auto t_y = yv.begin();
			for (auto i = 0; i < meshy.first.get_nrows(); ++i) {
				for (auto j = 0; j < meshy.first.get_ncols(); ++j) {
					*t_y++ = meshy.first.at(j, i) + meshy.second.at(j, i);
				}
			}
		}
		std::pair<Array2D<T>, Array2D<T> > meshxy = meshgrid(xv, yv);
		ArrayND<2, std::vector<T> > r2 = zeros_ND<2, T>({yv.size(), xv.size()});
		for (auto j = 0; j < meshxy.first.size(); ++j)r2[j]=pow(meshxy.first[j],2) + pow(meshxy.second[j],2);
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

#endif //NDEBUG
		return result;
	}
}
#endif //PRISM_PROJPOT_H
