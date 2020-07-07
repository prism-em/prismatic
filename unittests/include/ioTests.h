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

template <size_t N, class T>
ArrayND<N-1, T> subspace(ArrayND<N, T> &orig,  const size_t &index)
{
    //returns a reduced array of all values in array at index along slowest dimension
    std::array<size_t, N> inDims = orig.get_dimarr();
    std::array<size_t, N-1> outDims;
    size_t strides = 1;
    for(auto i = 1; i < N; i++)
    {
        strides *= inDims[i];
        outDims[i-1] = inDims[i];
    }

    ArrayND<N-1, T> output = zeros_ND<N-1, typename T::value_type>(outDims);
    std::copy(&orig[index*strides], &orig[(index+1)*strides], output.begin());

    return output;
};

template <size_t N>
ArrayND<N, std::vector<PRISMATIC_FLOAT_PRECISION>> getAmplitude(ArrayND<N, std::vector<std::complex<PRISMATIC_FLOAT_PRECISION>>> &input)
{
    ArrayND<N, std::vector<PRISMATIC_FLOAT_PRECISION>> output = zeros_ND<N, PRISMATIC_FLOAT_PRECISION>(input.get_dimarr());

    auto in_ptr = input.begin();
    for(auto& i : output) i = pow(std::abs(*in_ptr++), 2);

    return output;
};

// herr_t myOperator(void *elem, hid_t type_id, unsigned ndim, 
//                 const hsize_t *point, void *operator_data)
// {
//     herr_t a = 0;
//     return a;
// };

} //namespace Prismatic
#endif