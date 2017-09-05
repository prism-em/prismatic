// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// Prismatic is distributed under the GNU General Public License (GPL)
// If you use Prismatic, we kindly ask that you cite the following papers:

// 1. Ophus, C.: A fast image simulation algorithm for scanning
//    transmission electron microscopy. Advanced Structural and
//    Chemical Imaging 3(1), 13 (2017)

// 2. Pryor, Jr., A., Ophus, C., and Miao, J.: A Streaming Multi-GPU
//    Implementation of Image Simulation Algorithms for Scanning
//	  Transmission Electron Microscopy. arXiv:1706.08563 (2017)

#ifndef PRISMATIC_UTILITY_H
#define PRISMATIC_UTILITY_H
#include <vector>
#include <string>
#include <sstream>
#include <mutex>
#include <complex>
#include "defines.h"
#include "fftw3.h"
#include "configure.h"
namespace Prismatic {

	extern std::mutex fftw_plan_lock; // for synchronizing access to shared FFTW resources


	template<class T>
	std::vector<T> vecFromRange(const T &start, const T &step, const T &stop) {
		std::vector<T> result;
		for (auto i = start; i <= stop; i += step) {
			result.push_back(i);
		}
		if (result.empty())result.push_back(start);
		return result;
	};

	template<class T>
	Array1D <T> makeFourierCoords(const size_t &N, const T &pixel_size) {
		Array1D <T> result = zeros_ND<1, T>({{N}});
		long long nc = (size_t) floor((T) N / 2);

		T dp = 1 / (N * pixel_size);
		for (auto i = 0; i < N; ++i) {
			result[(nc + (size_t) i) % N] = (i - nc) * dp;
		}
		return result;
	};

	template<class T>
	Array2D <T> fftshift2(Array2D<T> arr) {
		Array2D<T> result(arr);
		const long sj = std::floor(arr.get_dimj() / 2);
		const long si = std::floor(arr.get_dimi() / 2);
		for (auto j = 0; j < arr.get_dimj(); ++j) {
			for (auto i = 0; i < arr.get_dimi(); ++i) {
				result.at((j + sj) % arr.get_dimj(),(i + si) % arr.get_dimi()) = arr.at(j,i);
			}
		}
		return result;
	};

	template<class T>
	std::string generateFilename(const Parameters <T> &pars, const size_t ay, const size_t ax) {
		std::string result = pars.meta.filenameOutput.substr(0, pars.meta.filenameOutput.find_last_of("."));
		std::stringstream ss;
		ss << "_X" << ax << "_Y" << ay << "_FP" << pars.meta.fpNum;
		//result += "_X" + std::string(ax) + "_Y" + std::string(ay) + "_FP" + std::string(pars.meta.fpNum);
		result += ss.str() + pars.meta.filenameOutput.substr(pars.meta.filenameOutput.find_last_of("."));
		return result;

	}

	std::pair<Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> >, Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> > >
	upsamplePRISMProbe(Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> > probe,
	                   const long dimj, const long dimi, long ys=0, long xs=0);

	PRISMATIC_FLOAT_PRECISION computePearsonCorrelation(Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> > left,
	                                                Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> > right);
	PRISMATIC_FLOAT_PRECISION computeRfactor(Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> > left,
	                                     Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> > right);
}

#endif //PRISMATIC_UTILITY_H
