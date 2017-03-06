//
// Created by AJ Pryor on 3/6/17.
//

#ifndef PRISM_UTILITY_H
#define PRISM_UTILITY_H
namespace PRISM{
	template<class T>
	vector<T> vecFromRange(const T &start, const T &step, const T &stop) {
		vector<T> result;
		for (auto i = start; i <= stop; i += step) {
			result.push_back(i);
		}
		return result;
	};
}
#endif //PRISM_UTILITY_H
