//
// Created by AJ Pryor on 2/24/17.
//

#include "PRISM02.h"
#include <iostream>
using namespace std;
namespace PRISM {
	template <class T>
	using Array3D = PRISM::ArrayND<3, std::vector<T> >;
	template <class T>
	using Array2D = PRISM::ArrayND<2, std::vector<T> >;
	template <class T>
	using Array1D = PRISM::ArrayND<1, std::vector<T> >;
	template <class T>
	Array1D<T> makeFourierCoords(const size_t& N, const T& pixel_size){
		Array1D<T> result = zeros_ND<1, T>({N});
		long long nc = (size_t)floor( (double)N/2 );

		T dp = 1/(N*pixel_size);
		for (auto i = 0; i < N; ++i){
			result[(nc + (size_t)i) % N] = (i-nc) * dp;
		}
		return result;
};

	template <class T>
	void PRISM02(emdSTEM<T>& pars){


		constexpr double m = 9.109383e-31;
		constexpr double e = 1.602177e-19;
		constexpr double c = 299792458;
		constexpr double h = 6.62607e-34;
		pars.imageSize[0] = pars.pot.get_nlayers();
		pars.imageSize[1] = pars.pot.get_ncols();
		Array1D<T> qx = makeFourierCoords(pars.imageSize[0], pars.pixelSize[0]);
		Array1D<T> qy = makeFourierCoords(pars.imageSize[1], pars.pixelSize[1]);

		pair< Array2D<T>, Array2D<T> > mesh = meshgrid(qy,qx);
		pars.qya = mesh.first;
		pars.qxa = mesh.second;
		Array2D<T> q2(pars.qya);
		transform(pars.qxa.begin(), pars.qxa.end(),
				  pars.qya.begin(), q2.begin(), [](const T& a, const T& b){
				return a*a + b*b;
				});

		// get qMax
		pars.qMax = 0;
		{
			T qx_max;
			T qy_max;
			for (auto i = 0; i < qx.size(); ++i) {
				qx_max = ( abs(qx[i]) > qx_max) ? abs(qx[i]) : qx_max;
				qy_max = ( abs(qy[i]) > qy_max) ? abs(qy[i]) : qy_max;
			}
			pars.qMax = min(qx_max, qy_max) / 2;
		}

		pars.qMask = zeros_ND<2, unsigned int>({pars.imageSize[1], pars.imageSize[1]});
		{
			int offset_x = pars.qMask.get_ncols()/2;
			int offset_y = pars.qMask.get_nrows()/2;

			// Fix this
			for (auto y = 0; y < pars.qMask.get_nrows() / 2; ++y) {
				for (auto x = 0; x < pars.qMask.get_ncols() / 2; ++x) {
					pars.qMask.at((y-offset_y) % pars.qMask.get_nrows(),((x-offset_x) % pars.qMask.get_ncols())) = 1;
				}
			}
		}
		pars.qMask.toMRC_f("/mnt/spareA/clion/PRISM/MATLAB/debug.mrc");
		cout << "pars.qMax = " << pars.qMax << endl;
		cout << "pars.qya.at(1,1) = " << pars.qya.at(1,1) << endl;
		cout << "pars.qya.at(0,1) = " << pars.qya.at(0,1) << endl;
		cout << "q2.at(3,4) = " << q2.at(3,4) << endl;
		cout << "q2.at(5,5) = " << q2.at(5,5) << endl;

		cout << "pars.pixelSize[0] = " << pars.pixelSize[0]<< endl;
		cout << "pars.pixelSize[1] = " << pars.pixelSize[1]<< endl;
		for (auto i = 0; i < 10; ++i){
			cout << "qx[" << i << "] = " << qx[i] << endl;
			cout << "qy[" << i << "] = " << qy[i] << endl;
		}
		cout << "qx[499] = " << qx[499] << endl;
		cout << "qy[499] = " << qy[499] << endl;
		cout << "qx[500] = " << qx[500] << endl;
		cout << "qy[500] = " << qy[500] << endl;

	}
}