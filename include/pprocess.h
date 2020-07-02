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

namespace Prismatic
{

static const PRISMATIC_FLOAT_PRECISION pi = acos(-1);

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



} //namespace Prismatic

#endif //PRISMATIC_PPROCESS_H