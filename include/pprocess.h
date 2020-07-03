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

#ifndef PRISMATIC_PROCESS_H
#define PRISMATIC_PROCESS_H
#include <vector>
#include <string>
#include <sstream>
#include <complex>
#include <ctime>
#include <iomanip>
#include "fftw3.h"
#include "ArrayND.h"
#include "utility.h"
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>

namespace Prismatic
{

boost::random::mt19937 urng(std::time(0)); //uniform random number generator
static const PRISMATIC_FLOAT_PRECISION pi = acos(-1);

//array utilities
template <size_t N, class T>
void scaleArray(ArrayND<N, std::vector<T>> &arr, PRISMATIC_FLOAT_PRECISION &scale)
{
    arr *= scale;
};

//dstribution sampling

PRISMATIC_FLOAT_PRECISION gaussian_sample(const PRISMATIC_FLOAT_PRECISION mu, 
                                    const PRISMATIC_FLOAT_PRECISION sigma, 
                                    const PRISMATIC_FLOAT_PRECISION sample)
{
    return exp(-(sample-mu)*(sample-mu)/(2.0*sigma*sigma)) / (sigma*sqrt(2*pi));	
};

Array1D<PRISMATIC_FLOAT_PRECISION> gaussian1D(const PRISMATIC_FLOAT_PRECISION mu, 
                                            const PRISMATIC_FLOAT_PRECISION sigma, 
                                            const Array1D<PRISMATIC_FLOAT_PRECISION> samples)
{
    Array1D<PRISMATIC_FLOAT_PRECISION> output(samples);
    for(auto i = 0; i < samples.get_dimi(); i++)
    {
        output.at(i) = gaussian_sample(mu, sigma, samples.at(i));
    }
    return output;
};

Array2D<PRISMATIC_FLOAT_PRECISION> gaussian2D(const PRISMATIC_FLOAT_PRECISION mu_x, 
                                            const PRISMATIC_FLOAT_PRECISION sigma_x,
                                            const PRISMATIC_FLOAT_PRECISION mu_y,
                                            const PRISMATIC_FLOAT_PRECISION sigma_y,
                                            const Array2D<PRISMATIC_FLOAT_PRECISION> samples_x,
                                            const Array2D<PRISMATIC_FLOAT_PRECISION> samples_y)
{
    Array2D<PRISMATIC_FLOAT_PRECISION> output(samples_x);
    for(auto j = 0; j < samples_x.get_dimj(); j++)
    {
        for(auto i = 0; i < samples_x.get_dimi(); i++)
        {
            output.at(j,i) = gaussian_sample(mu_x, sigma_x, samples_x.at(j, i))+gaussian_sample(mu_y, sigma_y, samples_y.at(j, i));
        }
    }

    return output;
};

PRISMATIC_FLOAT_PRECISION poisson_sample(const PRISMATIC_FLOAT_PRECISION lambda)
{
    //get distribution
    boost::random::poisson_distribution<int, PRISMATIC_FLOAT_PRECISION> pd(lambda);
    return (PRISMATIC_FLOAT_PRECISION) pd(urng);
}

template <size_t N>
void poisson_array(ArrayND<N, std::vector<PRISMATIC_FLOAT_PRECISION>> &arr)
{
    //assume that arr is scaled to array of scalar lambdas
    for (auto i = 0; i < arr.size(); i++) arr[i] = poisson_sample(arr[i]);
};

template <size_t N>
void applyPoisson(ArrayND<N, std::vector<PRISMATIC_FLOAT_PRECISION>> &arr, PRISMATIC_FLOAT_PRECISION &scale)
{
    scaleArray(arr,scale);
    poisson_array(arr);
};

template <size_t N>
void applyPoisson_norm(ArrayND<N, std::vector<PRISMATIC_FLOAT_PRECISION>> &arr, PRISMATIC_FLOAT_PRECISION &scale)
{
    scaleArray(arr,scale);
    poisson_array(arr);
    scaleArray(arr, 1.0/scale);
};

template <size_t N, class T>
ArrayND<N, T> subarray(ArrayND<N, T> &arr, std::array<size_t, N> starts, std::array<size_t, N> stops)
{
    //return an array of same dimension but cropped
    std::array<size_t, N-1> strides = arr.get_strides();
    std::array<int, N> dims;
    std::array<int, N> pdims;
    for(auto i = 0; i < N; i++) dims[i] = stops[i]-starts[i];
    pdims[N-1] = dims[N-1];
    for(auto i = N-2; i >= 0; i--) pdims[i] = pdims[i+1]*dims[i];

    int num = 1;
    for(auto i = 0; i < N; i++) num*=dims[i];

    ArrayND<N, T> sub = zeros_ND<N, typename T::value_type>(dims);

    int idx;
    std::array<int, N> out_idx;
    std::array<int, N> in_idx;
    for(auto i = 0; i < num; i++)
    {
        //calculate relative indices
        out_idx[N-1] = i % dims[N-1];
        in_idx[N-1] = (i % dims[N-1])+dims[N-1];
        for(auto j = N-2; j >=0; j--)
        {
            out_idx[j] = i / pdims[j];
            in_idx[j] = out_idx[j] + dims[j];
        }

        //use rel indices and offset to calcualte spot in 1D array
        idx = 0;
        for(auto j = 0; j < N-1; j++) idx += in_idx[j]*strides[j];
        idx += in_idx[N-1];
        sub[i] = arr[idx];
    }

    return sub;
};

template <size_t N, class T>
ArrayND<N, T> subslice(ArrayND<N, T> &arr)
{
    //return an array of N-1 dimensions
    return arr;
};

} //namespace Prismatic

#endif //PRISMATIC_PPROCESS_H