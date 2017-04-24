// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

#ifndef PRISM_UTILITY_H
#define PRISM_UTILITY_H
#include <vector>
#include <string>
#include <sstream>
namespace PRISM{
	template<class T>
	std::vector<T> vecFromRange(const T &start, const T &step, const T &stop) {
		//std::cout << "vecFromRange start = " << start << std::endl;
		//std::cout << "vecFromRange step = " << step << std::endl;
		//std::cout << "vecFromRange stop = " << stop << std::endl;
		std::vector<T> result;
		for (auto i = start; i <= stop; i += step) {
			result.push_back(i);
		}
		if (result.empty())result.push_back(start);
		return result;
	};
template <class T>
	Array1D<T> makeFourierCoords(const size_t &N, const T &pixel_size) {
		Array1D<T> result = zeros_ND<1, T>({{N}});
		long long nc = (size_t) floor((T) N / 2);

		T dp = 1 / (N * pixel_size);
		for (auto i = 0; i < N; ++i) {
			result[(nc + (size_t) i) % N] = (i - nc) * dp;
		}
		return result;
	};


	template <class T>
	std::string generateFilename(const Parameters<T>& pars, const size_t ay, const size_t ax){
		std::string result = pars.meta.filename_output.substr(0, pars.meta.filename_output.find_last_of("."));
		std::stringstream ss;
        ss << "_X"  << ax << "_Y" << ay << "_FP" << pars.meta.fpNum;
        //result += "_X" + std::string(ax) + "_Y" + std::string(ay) + "_FP" + std::string(pars.meta.fpNum);
		result += ss.str() + pars.meta.filename_output.substr(pars.meta.filename_output.find_last_of("."));
		return result;

	}

}
#endif //PRISM_UTILITY_H
