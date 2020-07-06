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

template<typename T>
Array1D<T> subarray(Array1D<T> &arr, size_t start, size_t stop)
{
    Array1D<T> output = zeros_ND<1, T>({{stop-start}});
    for(auto i = start; i < stop; i++)
    {
        output.at(i-start) = arr.at(i);
    }
    return output;
};

template<typename T>
Array2D<T> subarray(Array2D<T> &arr, std::array<size_t, 2> start, std::array<size_t, 2> stop)
{
    std::array<size_t, 2> dims;
    for(auto i = 0; i < 2; i++) dims[i] = stop[i] - start[i];

    Array2D<T> output = zeros_ND<2, T>(dims);
    for(auto j = start[0]; j < stop[0]; j++)
    {
        for(auto i = start[1]; i < stop[1]; i++)
        {
            output.at(j-start[0], i-start[1]) = arr.at(j, i);
        }
    }
    return output;
};

template<typename T>
Array3D<T> subarray(Array3D<T> &arr, std::array<size_t, 3> start, std::array<size_t, 3> stop)
{
    std::array<size_t, 3> dims;
    for(auto i = 0; i < 3; i++) dims[i] = stop[i] - start[i];

    Array3D<T> output = zeros_ND<3, T>(dims);
    for(auto k = start[0]; k < stop[0]; k++)
    {
        for(auto j = start[1]; j < stop[1]; j++)
        {
            for(auto i = start[2]; i < stop[2]; i++)
            {
                output.at(k-start[0], j-start[1], i-start[2]) = arr.at(k, j, i);
            }
        }
    }
    return output;
};

template<typename T>
Array4D<T> subarray(Array4D<T> &arr, std::array<size_t, 4> start, std::array<size_t, 4> stop)
{
    std::array<size_t, 4> dims;
    for(auto i = 0; i < 4; i++) dims[i] = stop[i] - start[i];

    Array4D<T> output = zeros_ND<4, T>(dims);
    for(auto l = start[0]; l < stop[0]; l++)
    {
        for(auto k = start[1]; k < stop[1]; k++)
        {
            for(auto j = start[2]; j < stop[2]; j++)
            {
                for(auto i = start[3]; i < stop[3]; i++)
                {
                    output.at(l-start[0], k-start[1], j-start[2], i-start[3]) = arr.at(l, k, j, i);
                }
            }
        }
    }
    return output;
};

template<typename T>
Array1D<T> subslice(Array2D<T> &arr, size_t dim, size_t idx)
{
    //return an array of N-1 dimensions where arr is indexed along a single dimension
    Array1D<T> output;
    switch (dim)
    {
        case 1: //freeze dim j
            output = zeros_ND<1, T>({arr.get_dimi()});
            for(auto i = 0; i < arr.get_dimi(); i++)
                output.at(i) = arr.at(idx, i);
            break;
        case 0: //freeze dim i
            output = zeros_ND<1, T>({arr.get_dimj()});
            for(auto i = 0; i < arr.get_dimj(); i++)
                output.at(i) = arr.at(i, idx);
            break;
    }
    return output;
};

template<typename T>
Array2D<T> subslice(Array3D<T> &arr, size_t dim, size_t idx)
{
    //return an array of N-1 dimensions where arr is indexed along a single dimension
    //maintain order of dimensions otherwise
    Array2D<T> output;
    switch (dim)
    {
        case 2: //freeze dim k
            output = zeros_ND<2, T>({arr.get_dimj(), arr.get_dimi()});
            for(auto j = 0; j < arr.get_dimj(); j++)
            {
                for(auto i = 0; i < arr.get_dimi(); i++)
                {
                    output.at(j, i) = arr.at(idx, j, i);
                }
            }
            break;
        case 1: //freeze dim j
            output = zeros_ND<2, T>({arr.get_dimk(), arr.get_dimi()});
            for(auto j = 0; j < arr.get_dimk(); j++)
            {
                for(auto i = 0; i < arr.get_dimi(); i++)
                {
                    output.at(j, i) = arr.at(j, idx, i);
                }
            }
            break;
        case 0: //freeze dim i
            output = zeros_ND<2, T>({arr.get_dimk(), arr.get_dimj()});
            for(auto j = 0; j < arr.get_dimk(); j++)
            {
                for(auto i = 0; i < arr.get_dimj(); i++)
                {
                    output.at(j, i) = arr.at(j, i, idx);
                }
            }
            break;
    }
    return output;
};

template<typename T>
Array3D<T> subslice(Array4D<T> &arr, size_t dim, size_t idx)
{
    //return an array of N-1 dimensions where arr is indexed along a single dimension
    //maintain order of dimensions otherwise
    Array3D<T> output;
    switch (dim)
    {
        case 3: //freeze dim l
            output = zeros_ND<3, T>({arr.get_dimk(), arr.get_dimj(), arr.get_dimi()});
            for(auto k = 0; k < arr.get_dimk(); k++)
            {
                for(auto j = 0; j < arr.get_dimj(); j++)
                {
                    for(auto i = 0; i < arr.get_dimi(); i++)
                    {
                        output.at(k, j, i) = arr.at(idx, k, j, i);
                    }
                }
            }
            break;
        case 2: //freeze dim k
            output = zeros_ND<3, T>({arr.get_diml(), arr.get_dimj(), arr.get_dimi()});
            for(auto k = 0; k < arr.get_diml(); k++)
            {
                for(auto j = 0; j < arr.get_dimj(); j++)
                {
                    for(auto i = 0; i < arr.get_dimi(); i++)
                    {
                        output.at(k, j, i) = arr.at(k, idx, j, i);
                    }
                }
            }
            break;
        case 1: //freeze dim j
            output = zeros_ND<3, T>({arr.get_diml(), arr.get_dimk(), arr.get_dimi()});
            for(auto k = 0; k < arr.get_diml(); k++)
            {
                for(auto j = 0; j < arr.get_dimk(); j++)
                {
                    for(auto i = 0; i < arr.get_dimi(); i++)
                    {
                        output.at(k, j, i) = arr.at(k, j, idx, i);
                    }
                }
            }
            break;
        case 0: //freeze dim i
            output = zeros_ND<3, T>({arr.get_diml(), arr.get_dimk(), arr.get_dimj()});
            for(auto k = 0; k < arr.get_diml(); k++)
            {
                for(auto j = 0; j < arr.get_dimk(); j++)
                {
                    for(auto i = 0; i < arr.get_dimj(); i++)
                    {
                        output.at(k, j, i) = arr.at(k, j, i, idx);
                    }
                }
            }
            break;
    }
    return output;
};

} //namespace Prismatic

#endif //PRISMATIC_PPROCESS_H