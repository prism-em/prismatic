//
// Created by AJ Pryor on 2/22/17.
//

#ifndef PRISM_PROJPOT_H
#define PRISM_PROJPOT_H
#include "ArrayND.h"
#include <math.h>
#include <iostream>
namespace PRISM {
	using namespace std;
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

#ifndef NDEBUG
#include <iostream>
		using namespace std;
		cout << "a0 = " << a0 << endl;
		cout << "term1 = " << term1 << endl;
		cout << "term2 = " << term2 << endl;
		cout << "term1 = " << term1 << endl;

#endif //NDEBUG
		return result;
	}
}
#endif //PRISM_PROJPOT_H
