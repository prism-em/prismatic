#ifndef IO_TESTS_H
#define IO_TESTS_H
#include <random>
#include "ArrayND.h"
#include <vector>
#include "params.h"

namespace Prismatic{

template <size_t N>
void assignRandomValues(ArrayND<N, std::vector<PRISMATIC_FLOAT_PRECISION>> &arr, std::default_random_engine &de)
{
    std::normal_distribution<PRISMATIC_FLOAT_PRECISION> randn(0,1);
    for (auto i = 0; i < arr.size(); i++) arr[i] = randn(de);
        
};

template <size_t N>
void assignRandomValues(ArrayND<N, std::vector<std::complex<PRISMATIC_FLOAT_PRECISION>>> &arr, std::default_random_engine &de)
{
    std::normal_distribution<PRISMATIC_FLOAT_PRECISION> randn(0,1);
    for (auto i = 0; i < arr.size(); i++)  arr[i] = std::complex<PRISMATIC_FLOAT_PRECISION>(randn(de), randn(de));
};

template <size_t N, class T>
bool compareSize(ArrayND<N, T> &ref, ArrayND<N, T> &test)
{  
    std::array<size_t, N> ref_dims = ref.get_dimarr();
    std::array<size_t, N> test_dims = test.get_dimarr();
    return ref_dims == test_dims;
};

template <size_t N, class T>
PRISMATIC_FLOAT_PRECISION compareValues(ArrayND<N, T> &ref, ArrayND<N, T> &test)
{
    PRISMATIC_FLOAT_PRECISION errorSum = 0.0;
    for(auto i = 0; i < ref.size(); i++) errorSum += std::abs(ref[i]-test[i]);
    return errorSum;
};

} //namespace Prismatic
#endif