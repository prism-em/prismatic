//
// Created by AJ Pryor on 3/6/17.
//

#ifndef PRISM_UTILITY_H
#define PRISM_UTILITY_H
#include <vector>
namespace PRISM{
	template<class T>
	std::vector<T> vecFromRange(const T &start, const T &step, const T &stop) {
		std::vector<T> result;
		for (auto i = start; i <= stop; i += step) {
			result.push_back(i);
		}
		return result;
	};
}
#endif //PRISM_UTILITY_H
