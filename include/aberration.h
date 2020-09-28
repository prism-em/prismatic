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

#ifndef PRISM_ABERRATION_H
#define PRISM_ABERRATION_H
#include <string>
#include <vector>
#include <array>
#include <iostream>
#include "ArrayND.h"
#include "defines.h"

struct aberration
{
    int m;
    int n;
    PRISMATIC_FLOAT_PRECISION mag;
    PRISMATIC_FLOAT_PRECISION angle;
    void to_string() const
    {
        std::cout << "m = " << m << std::endl;
        std::cout << "n = " << n << std::endl;
        std::cout << "mag = " << mag << std::endl;
        std::cout << "angle = " << angle << std::endl;
    };

    bool operator==(const aberration &a) const
    {
        return (m == a.m) && (n == a.n);
    };
};


namespace Prismatic
{

template <class T>
using Array2D = Prismatic::ArrayND<2, std::vector<T> >;

std::vector<aberration> readAberrations(const std::string &filename);

Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> getChi(Array2D<PRISMATIC_FLOAT_PRECISION> &q,
                                                        Array2D<PRISMATIC_FLOAT_PRECISION> &qTheta,
                                                        PRISMATIC_FLOAT_PRECISION &lambda, 
                                                        std::vector<aberration> &ab);

std::vector<aberration> updateAberrations(std::vector<aberration> ab, 
                                        PRISMATIC_FLOAT_PRECISION C1, 
                                        PRISMATIC_FLOAT_PRECISION C3, 
                                        PRISMATIC_FLOAT_PRECISION C5,
                                        const PRISMATIC_FLOAT_PRECISION &lambda);
                        
} // namespace Prismatic
#endif //PRISM_ABERRATION_H