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
#include <mutex>

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

Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> gaussian2D_k(const PRISMATIC_FLOAT_PRECISION mu_x, 
                                            const PRISMATIC_FLOAT_PRECISION sigma_x,
                                            const PRISMATIC_FLOAT_PRECISION mu_y,
                                            const PRISMATIC_FLOAT_PRECISION sigma_y,
                                            const Array2D<PRISMATIC_FLOAT_PRECISION> samples_x,
                                            const Array2D<PRISMATIC_FLOAT_PRECISION> samples_y)
{
    //return a fourier transformed gaussian kernel
    Array2D<PRISMATIC_FLOAT_PRECISION> tmp(samples_x);
    for(auto j = 0; j < samples_x.get_dimj(); j++)
    {
        for(auto i = 0; i < samples_x.get_dimi(); i++)
        {
            tmp.at(j,i) = gaussian_sample(mu_x, sigma_x, samples_x.at(j, i))+gaussian_sample(mu_y, sigma_y, samples_y.at(j, i));
        }
    }

    //prepare FFT calls
    extern std::mutex fftw_plan_lock;
    Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> output = zeros_ND<2, std::complex<PRISMATIC_FLOAT_PRECISION>>({{samples_x.get_dimj(), samples_x.get_dimi()}});

	//create FFT plans 
	PRISMATIC_FFTW_INIT_THREADS();
	PRISMATIC_FFTW_PLAN_WITH_NTHREADS(1);
	
	std::unique_lock<std::mutex> gatekeeper(fftw_plan_lock);
	PRISMATIC_FFTW_PLAN plan_forward = PRISMATIC_FFTW_PLAN_DFT_2D(output.get_dimj(), output.get_dimi(),
															reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&output[0]),
															reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&output[0]),
															FFTW_FORWARD,
															FFTW_ESTIMATE);

	gatekeeper.unlock();

    //copy data to transform
    for(auto i = 0; i < output.size(); i++){
         output[i] = tmp[i];
    }
    
    //transform, multiply, transform
    PRISMATIC_FFTW_EXECUTE(plan_forward);

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

template<typename T>
Array2D<T> bin(Array2D<T> &arr, size_t fi, size_t fj)
{
    //check for divisibility
    size_t extent_i = (arr.get_dimi() % fi) ? arr.get_dimi() - arr.get_dimi() % fi : arr.get_dimi();
    size_t extent_j = (arr.get_dimj() % fj) ? arr.get_dimj() - arr.get_dimj() % fj : arr.get_dimj();

    Array2D<T> output = zeros_ND<2, T>({{extent_j / fj, extent_i / fi}});
    for(auto j = 0; j < extent_j; j++)
    {
        for(auto i =0; i < extent_i; i++)
        {
            output.at(j/fj, i/fi) += arr.at(j,i);
        }
    }

    return output;
};

template<typename T>
Array3D<T> bin(Array3D<T> &arr, size_t fi, size_t fj, size_t fk)
{
    //check for divisibility
    size_t extent_i = (arr.get_dimi() % fi) ? arr.get_dimi() - arr.get_dimi() % fi : arr.get_dimi();
    size_t extent_j = (arr.get_dimj() % fj) ? arr.get_dimj() - arr.get_dimj() % fj : arr.get_dimj();
    size_t extent_k = (arr.get_dimk() % fk) ? arr.get_dimk() - arr.get_dimk() % fk : arr.get_dimk();

    Array3D<T> output = zeros_ND<3, T>({{extent_k / fk, extent_j / fj, extent_i / fi}});
    for(auto k = 0; k < extent_k; k++)
    {
        for(auto j = 0; j < extent_j; j++)
        {
            for(auto i =0; i < extent_i; i++)
            {
                output.at(k/fk, j/fj, i/fi) += arr.at(k,j,i);
            }
        }
    }

    return output;
};

template<typename T>
Array4D<T> bin(Array4D<T> &arr, size_t fi, size_t fj, size_t fk, size_t fl)
{
    //check for divisibility
    size_t extent_i = (arr.get_dimi() % fi) ? arr.get_dimi() - arr.get_dimi() % fi : arr.get_dimi();
    size_t extent_j = (arr.get_dimj() % fj) ? arr.get_dimj() - arr.get_dimj() % fj : arr.get_dimj();
    size_t extent_k = (arr.get_dimk() % fk) ? arr.get_dimk() - arr.get_dimk() % fk : arr.get_dimk();
    size_t extent_l = (arr.get_diml() % fl) ? arr.get_diml() - arr.get_diml() % fl : arr.get_diml();

    Array4D<T> output = zeros_ND<4, T>({{extent_l / fl, extent_k / fk, extent_j / fj, extent_i / fi}});
    for(auto l = 0; l < extent_l; l++)
    {
        for(auto k = 0; k < extent_k; k++)
        {
            for(auto j = 0; j < extent_j; j++)
            {
                for(auto i =0; i < extent_i; i++)
                {
                    output.at(l/fl, k/fk, j/fj, i/fi) += arr.at(l,k,j,i);
                }
            }
        }
    }

    return output;
};

template<size_t N>
ArrayND<N, std::vector<PRISMATIC_FLOAT_PRECISION>> getAmp(ArrayND<N, std::vector<std::complex<PRISMATIC_FLOAT_PRECISION>>> &arr)
{
    ArrayND<N, std::vector<PRISMATIC_FLOAT_PRECISION>> output = zeros_ND<N, PRISMATIC_FLOAT_PRECISION>(arr.get_dimarr());

    for(auto i = 0; i < arr.size(); i++) output[i] = pow(std::abs(arr[i]), 2.0);

    return output;
};

template<size_t N>
ArrayND<N, std::vector<PRISMATIC_FLOAT_PRECISION>> getPhase(ArrayND<N, std::vector<std::complex<PRISMATIC_FLOAT_PRECISION>>> &arr)
{
    ArrayND<N, std::vector<PRISMATIC_FLOAT_PRECISION>> output = zeros_ND<N, PRISMATIC_FLOAT_PRECISION>(arr.get_dimarr());

    for(auto i = 0; i < arr.size(); i++) output[i] = std::arg(arr[i]);
    return output;
};

Array2D<PRISMATIC_FLOAT_PRECISION> fourierDownsample(Array2D<PRISMATIC_FLOAT_PRECISION> &arr, int Ni, int Nj)
{
    extern std::mutex fftw_plan_lock;

    Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> fstore = zeros_ND<2,std::complex<PRISMATIC_FLOAT_PRECISION>>({{arr.get_dimj(), arr.get_dimi()}});
	Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> bstore = zeros_ND<2,std::complex<PRISMATIC_FLOAT_PRECISION>>({{Nj, Ni}});
	Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> farr = zeros_ND<2,std::complex<PRISMATIC_FLOAT_PRECISION>>({{arr.get_dimj(),arr.get_dimi()}});
	Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> barr = zeros_ND<2,std::complex<PRISMATIC_FLOAT_PRECISION>>({{Nj, Ni}});
    Array2D<PRISMATIC_FLOAT_PRECISION> result = zeros_ND<2, PRISMATIC_FLOAT_PRECISION>({{Nj, Ni}});

	//create FFT plans 
	PRISMATIC_FFTW_INIT_THREADS();
	PRISMATIC_FFTW_PLAN_WITH_NTHREADS(1);
	
	std::unique_lock<std::mutex> gatekeeper(fftw_plan_lock);
	PRISMATIC_FFTW_PLAN plan_forward = PRISMATIC_FFTW_PLAN_DFT_2D(fstore.get_dimj(), fstore.get_dimi(),
															reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&farr[0]),
															reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&fstore[0]),
															FFTW_FORWARD,
															FFTW_ESTIMATE);

	PRISMATIC_FFTW_PLAN plan_inverse = PRISMATIC_FFTW_PLAN_DFT_2D(bstore.get_dimj(), bstore.get_dimi(),
															reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&bstore[0]),
															reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&barr[0]),
															FFTW_BACKWARD,
															FFTW_ESTIMATE);
	gatekeeper.unlock();

    //calculate indices for downsampling in fourier space
	int nyqi = std::floor(Ni/2) + 1;
	int nyqj = std::floor(Nj/2) + 1;

    //copy data to forward transform
    for(auto i = 0; i < farr.size(); i++) farr[i] = arr[i];
    
    //forward transform 
    PRISMATIC_FFTW_EXECUTE(plan_forward);

    //copy relevant quadrants to backward store
    //manual looping through quadrants
    for(auto j = 0; j < nyqj; j++)
    {
        for(auto i = 0; i < nyqi; i++)
        {
            bstore.at(j, i) = fstore.at(j, i);
        }
    }

    for(auto j = nyqj-Nj; j < 0; j++)
    {
        for(auto i = 0; i < nyqi; i++)
        {
            bstore.at(Nj + j, i) = fstore.at(fstore.get_dimj() + j, i);
        }
    }

    for(auto j = 0; j < nyqj; j++)
    {
        for(auto i = nyqi-Ni; i < 0; i++)
        {
            bstore.at(j, Ni + i) = fstore.at(j, fstore.get_dimi() + i);
        }
    }

    for(auto j = nyqj-Nj; j < 0; j++)
    {
        for(auto i = nyqi-Ni; i < 0; i++)
        {
            bstore.at(Nj + j, Ni + i) = fstore.at(fstore.get_dimj() + j, fstore.get_dimi() + i);
        }
    }

    //inverse transform
    PRISMATIC_FFTW_EXECUTE(plan_inverse);

    //store slice in potential
    for(auto i = 0; i < barr.size(); i++) result[i] = barr[i].real();

	PRISMATIC_FLOAT_PRECISION orig_x = arr.get_dimi();
	PRISMATIC_FLOAT_PRECISION orig_y = arr.get_dimj();
	PRISMATIC_FLOAT_PRECISION new_x = Ni;
	PRISMATIC_FLOAT_PRECISION new_y = Nj;
	result /= Ni*Nj;
	result *= (new_x/orig_x)*(new_y/orig_y);

    return result;
};


void convolve2D(Array2D<PRISMATIC_FLOAT_PRECISION> &arr, Array2D<PRISMATIC_FLOAT_PRECISION> &kernel)
{   
    //convolve two arrays using fourier transform method
    //result is stored in arr
    //enforce size equivalence
    if(arr.get_dimi() != kernel.get_dimi() || arr.get_dimj() != kernel.get_dimj()) return;

    //prepare FFT calls
    extern std::mutex fftw_plan_lock;

    Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> karr = zeros_ND<2,std::complex<PRISMATIC_FLOAT_PRECISION>>({{arr.get_dimj(), arr.get_dimi()}});
	Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> kkern = zeros_ND<2,std::complex<PRISMATIC_FLOAT_PRECISION>>({{kernel.get_dimj(),kernel.get_dimi()}});

	//create FFT plans 
	PRISMATIC_FFTW_INIT_THREADS();
	PRISMATIC_FFTW_PLAN_WITH_NTHREADS(1);
	
	std::unique_lock<std::mutex> gatekeeper(fftw_plan_lock);
	PRISMATIC_FFTW_PLAN plan_forward_arr = PRISMATIC_FFTW_PLAN_DFT_2D(karr.get_dimj(), karr.get_dimi(),
															reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&karr[0]),
															reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&karr[0]),
															FFTW_FORWARD,
															FFTW_ESTIMATE);

	PRISMATIC_FFTW_PLAN plan_inv_arr = PRISMATIC_FFTW_PLAN_DFT_2D(karr.get_dimj(), karr.get_dimi(),
															reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&karr[0]),
															reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&karr[0]),
															FFTW_BACKWARD,
															FFTW_ESTIMATE);

	PRISMATIC_FFTW_PLAN plan_forward_kern = PRISMATIC_FFTW_PLAN_DFT_2D(kkern.get_dimj(), kkern.get_dimi(),
															reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&kkern[0]),
															reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&kkern[0]),
															FFTW_BACKWARD,
															FFTW_ESTIMATE);
	gatekeeper.unlock();

    //copy data to transform
    for(auto i = 0; i < arr.size(); i++){
         karr[i] = arr[i];
         kkern[i] = kernel[i];
    }
    
    //transform, multiply, transform
    PRISMATIC_FFTW_EXECUTE(plan_forward_arr);
    PRISMATIC_FFTW_EXECUTE(plan_forward_kern);

    for(auto i = 0; i < karr.size(); i++) karr[i] *= kkern[i];

    PRISMATIC_FFTW_EXECUTE(plan_inv_arr);

    //copy data back and scale
    for(auto i = 0; i < arr.size(); i++) arr[i] = karr[i].real();

    arr/=arr.get_dimi()*arr.get_dimj();
};

void convolve2D(Array2D<PRISMATIC_FLOAT_PRECISION> &arr, Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> &kkernel)
{   
    //convolve two arrays using fourier transform method
    //assumes input kernel has already been transformed
    //result is stored in arr
    //enforce size equivalence
    if(arr.get_dimi() != kkernel.get_dimi() || arr.get_dimj() != kkernel.get_dimj()) return;

    //prepare FFT calls
    extern std::mutex fftw_plan_lock;

    Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> karr = zeros_ND<2,std::complex<PRISMATIC_FLOAT_PRECISION>>({{arr.get_dimj(), arr.get_dimi()}});

	//create FFT plans 
	PRISMATIC_FFTW_INIT_THREADS();
	PRISMATIC_FFTW_PLAN_WITH_NTHREADS(1);
	
	std::unique_lock<std::mutex> gatekeeper(fftw_plan_lock);
	PRISMATIC_FFTW_PLAN plan_forward_arr = PRISMATIC_FFTW_PLAN_DFT_2D(karr.get_dimj(), karr.get_dimi(),
															reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&karr[0]),
															reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&karr[0]),
															FFTW_FORWARD,
															FFTW_ESTIMATE);

	PRISMATIC_FFTW_PLAN plan_inv_arr = PRISMATIC_FFTW_PLAN_DFT_2D(karr.get_dimj(), karr.get_dimi(),
															reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&karr[0]),
															reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&karr[0]),
															FFTW_BACKWARD,
															FFTW_ESTIMATE);

	gatekeeper.unlock();

    //copy data to transform
    for(auto i = 0; i < arr.size(); i++){
         karr[i] = arr[i];
    }
    
    //transform, multiply, transform
    PRISMATIC_FFTW_EXECUTE(plan_forward_arr);

    for(auto i = 0; i < karr.size(); i++) karr[i] *= kkernel[i];

    PRISMATIC_FFTW_EXECUTE(plan_inv_arr);

    //copy data back and scale
    for(auto i = 0; i < arr.size(); i++) arr[i] = karr[i].real();

    arr/=arr.get_dimi()*arr.get_dimj();
};

} //namespace Prismatic

#endif //PRISMATIC_PPROCESS_H