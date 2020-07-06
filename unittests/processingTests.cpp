#include <boost/test/unit_test.hpp>
#include "ArrayND.h"
#include <iostream>
#include <vector>
#include <stdio.h>
#include <random>
#include "utility.h"
#include "fileIO.h"
#include "pprocess.h"
#include "ioTests.h"

namespace Prismatic{

BOOST_AUTO_TEST_SUITE(processingTests);

BOOST_AUTO_TEST_CASE(poissonNoise)
{
    Array4D<PRISMATIC_FLOAT_PRECISION> testArr = zeros_ND<4,PRISMATIC_FLOAT_PRECISION>({{2,3,5,7}});
    for(auto i =0; i < testArr.size(); i++) testArr[i] = i+1;

    PRISMATIC_FLOAT_PRECISION scale = 1.0;
    applyPoisson(testArr, scale);
    
}

BOOST_AUTO_TEST_CASE(subindexing)
{
    //prepare test arrays
    int seed = 10101;
    srand(seed);
    std::default_random_engine de(seed);

    Array2D<PRISMATIC_FLOAT_PRECISION> testArr2D = zeros_ND<2,PRISMATIC_FLOAT_PRECISION>({{3,5}});
    Array3D<PRISMATIC_FLOAT_PRECISION> testArr3D = zeros_ND<3,PRISMATIC_FLOAT_PRECISION>({{3,5,7}});
    Array4D<PRISMATIC_FLOAT_PRECISION> testArr4D = zeros_ND<4,PRISMATIC_FLOAT_PRECISION>({{3,5,7,11}});
    assignRandomValues(testArr2D, de);
    assignRandomValues(testArr3D, de);
    assignRandomValues(testArr4D, de);

    //test subarray
    std::array<size_t, 2> starts2D = {1,2};
    std::array<size_t, 2> stops2D = {2,3};
    Array2D<PRISMATIC_FLOAT_PRECISION> sub2D = subarray(testArr2D, starts2D, stops2D);

    std::array<size_t, 3> starts3D = {1,2,3};
    std::array<size_t, 3> stops3D = {2,3,5};
    Array3D<PRISMATIC_FLOAT_PRECISION> sub3D = subarray(testArr3D, starts3D, stops3D);

    std::array<size_t, 4> starts4D = {1,2,3,4};
    std::array<size_t, 4> stops4D = {2,3,5,7};
    Array4D<PRISMATIC_FLOAT_PRECISION> sub4D = subarray(testArr4D, starts4D, stops4D);

    PRISMATIC_FLOAT_PRECISION err = 0.0;
    PRISMATIC_FLOAT_PRECISION tol = 0.001;
    for(auto j = starts2D[0]; j < stops2D[0]; j++)
    {
        for(auto i = starts2D[1]; i < stops2D[1]; i++)
        {
            err += std::abs(testArr2D.at(j,i) - sub2D.at(j-starts2D[0], i-starts2D[1]));
        }
    }
    BOOST_TEST(err < tol);
    err = 0.0;

    for(auto k = starts3D[0]; k < stops3D[0]; k++)
    {
        for(auto j = starts3D[1]; j < stops3D[1]; j++)
        {
            for(auto i = starts3D[2]; i < stops3D[2]; i++)
            {
                err += std::abs(testArr3D.at(k,j,i) - sub3D.at(k-starts3D[0], j-starts3D[1], i-starts3D[2]));
            }
        }
    }
    BOOST_TEST(err < tol);
    err = 0.0;

    for(auto l = starts4D[0]; l < stops4D[0]; l++)
    {
        for(auto k = starts4D[1]; k < stops4D[1]; k++)
        {
            for(auto j = starts4D[2]; j < stops4D[2]; j++)
            {
                for(auto i = starts4D[3]; i < stops4D[3]; i++)
                {
                    err += std::abs(testArr4D.at(l, k,j,i) - sub4D.at(l-starts4D[0], k-starts4D[1], j-starts4D[2], i-starts4D[3]));
                }
            }
        }
    }
    BOOST_TEST(err < tol);
    err = 0.0;

    //test subslicing
    size_t dim = 0;
    size_t idx = 3;

    Array1D<PRISMATIC_FLOAT_PRECISION> slice2D = subslice(testArr2D, dim, idx);
    Array2D<PRISMATIC_FLOAT_PRECISION> slice3D = subslice(testArr3D, dim, idx);
    Array3D<PRISMATIC_FLOAT_PRECISION> slice4D = subslice(testArr4D, dim, idx);
    for(auto j = 0; j < testArr2D.get_dimj(); j++)
    {
        err += std::abs(testArr2D.at(j,idx) - slice2D.at(j));
    }
    BOOST_TEST(err < tol);
    err = 0.0;

    for(auto k = 0; k < testArr3D.get_dimk(); k++)
    {
        for(auto j = 0; j < testArr3D.get_dimj(); j++)
        {
            err += std::abs(testArr3D.at(k,j,idx) - slice3D.at(k,j));
        }
    }
    BOOST_TEST(err < tol);
    err = 0.0;

    for(auto l = 0; l < testArr4D.get_diml(); l++)
    {
        for(auto k = 0; k < testArr4D.get_dimk(); k++)
        {
            for(auto j = 0; j < testArr4D.get_dimj(); j++)
            {
                err += std::abs(testArr4D.at(l,k,j,idx) - slice4D.at(l,k,j));
            }
        }
    }
    BOOST_TEST(err < tol);
    err = 0.0;
}

BOOST_AUTO_TEST_CASE(binning)
{
    //prepare test arrays
    int seed = 10101;
    srand(seed);
    std::default_random_engine de(seed);

    //square case
    size_t Rx = 16; size_t Ry = 16; size_t Qx = 64; size_t Qy = 64;
    Array2D<PRISMATIC_FLOAT_PRECISION> testArr2D = zeros_ND<2,PRISMATIC_FLOAT_PRECISION>({{Rx,Ry}});
    Array3D<PRISMATIC_FLOAT_PRECISION> testArr3D = zeros_ND<3,PRISMATIC_FLOAT_PRECISION>({{Rx,Ry,Qx}});
    Array4D<PRISMATIC_FLOAT_PRECISION> testArr4D = zeros_ND<4,PRISMATIC_FLOAT_PRECISION>({{Rx,Ry,Qx,Qy}});
    assignRandomValues(testArr2D, de);
    assignRandomValues(testArr3D, de);

    assignRandomValues(testArr4D, de);
    Array2D<PRISMATIC_FLOAT_PRECISION> bin2D = bin(testArr2D, 4, 4);
    Array3D<PRISMATIC_FLOAT_PRECISION> bin3D = bin(testArr3D, 8, 4, 4);
    Array4D<PRISMATIC_FLOAT_PRECISION> bin4D = bin(testArr4D, 8, 8, 4, 4);

    //check for arr size
    std::array<size_t, 2> bin_dims2D = {Rx/4, Ry/4};
    std::array<size_t, 3> bin_dims3D = {Rx/4, Ry/4, Qx/8};
    std::array<size_t, 4> bin_dims4D = {Rx/4, Ry/4, Qx/8, Qy/8};
    BOOST_TEST(bin_dims2D == bin2D.get_dimarr());
    BOOST_TEST(bin_dims3D == bin3D.get_dimarr());
    BOOST_TEST(bin_dims4D == bin4D.get_dimarr());

    //check for loss 
    PRISMATIC_FLOAT_PRECISION test2D_sum = 0.0;
    PRISMATIC_FLOAT_PRECISION test3D_sum = 0.0;
    PRISMATIC_FLOAT_PRECISION test4D_sum = 0.0;
    for(auto i = 0; i < Rx*Ry; i++) test2D_sum += testArr2D[i];
    for(auto i = 0; i < Rx*Ry*Qx; i++) test3D_sum += testArr3D[i];
    for(auto i = 0; i < Rx*Ry*Qx*Qy; i++) test4D_sum += testArr4D[i];

    PRISMATIC_FLOAT_PRECISION bin2D_sum = 0.0;
    PRISMATIC_FLOAT_PRECISION bin3D_sum = 0.0;
    PRISMATIC_FLOAT_PRECISION bin4D_sum = 0.0;
    for(auto i = 0; i < (Rx/4)*(Ry/4); i++) bin2D_sum += bin2D[i];
    for(auto i = 0; i < (Rx/4)*(Ry/4)*(Qx/8); i++) bin3D_sum += bin3D[i];
    for(auto i = 0; i < (Rx/4)*(Ry/4)*(Qx/8)*(Qy/8); i++) bin4D_sum += bin4D[i];
    
    PRISMATIC_FLOAT_PRECISION tol = 0.001;
    BOOST_TEST(std::abs(test2D_sum-bin2D_sum) < tol);
    BOOST_TEST(std::abs(test3D_sum-bin3D_sum) < tol);
    BOOST_TEST(std::abs(test4D_sum-bin4D_sum) < tol);

    //uneven rectangle
    Rx = 17; Ry = 19; Qx = 65; Qy = 67;
    testArr2D = zeros_ND<2,PRISMATIC_FLOAT_PRECISION>({{Rx,Ry}});
    testArr3D = zeros_ND<3,PRISMATIC_FLOAT_PRECISION>({{Rx,Ry,Qx}});
    testArr4D = zeros_ND<4,PRISMATIC_FLOAT_PRECISION>({{Rx,Ry,Qx,Qy}});
    assignRandomValues(testArr2D, de);
    assignRandomValues(testArr3D, de);
    assignRandomValues(testArr4D, de);

    bin2D = bin(testArr2D, 4, 4);
    bin3D = bin(testArr3D, 8, 4, 4);
    bin4D = bin(testArr4D, 8, 8, 4, 4);

    //only check size
    BOOST_TEST(bin_dims2D == bin2D.get_dimarr());
    BOOST_TEST(bin_dims3D == bin3D.get_dimarr());
    BOOST_TEST(bin_dims4D == bin4D.get_dimarr());

}

BOOST_AUTO_TEST_CASE(complexDisplay)
{
    Array1D<std::complex<PRISMATIC_FLOAT_PRECISION>> testArr = zeros_ND<1,std::complex<PRISMATIC_FLOAT_PRECISION>>({5});
    testArr[0] = {1.0, 0.0};
    testArr[1] = {0.0, 1.0};
    testArr[2] = {-1.0, 0.0};
    testArr[3] = {0.0, -1.0};
    testArr[4] = {1.0, 1.0};

    Array1D<PRISMATIC_FLOAT_PRECISION> testArr_amp = getAmp(testArr);
    Array1D<PRISMATIC_FLOAT_PRECISION> testArr_phase = getPhase(testArr);

    PRISMATIC_FLOAT_PRECISION tol = 0.00001;
    BOOST_TEST(std::abs(testArr_amp[0] - 1.0) < tol);
    BOOST_TEST(std::abs(testArr_amp[4] - 2.0) < tol);

    PRISMATIC_FLOAT_PRECISION pi = acos(-1);
    BOOST_TEST(testArr_phase[0] == 0.0);
    BOOST_TEST(testArr_phase[1] == pi/2);
    BOOST_TEST(testArr_phase[2] == pi);
    BOOST_TEST(testArr_phase[3] == -pi/2);
    BOOST_TEST(testArr_phase[4] == pi/4);

}

BOOST_AUTO_TEST_CASE(downsample)
{
    //only need to test that the downsampling runs; method is a direct lift from tested potential downsampling
    size_t Rx = 100; size_t Ry = 100;
    size_t Tx = 80; size_t Ty = 90;
    Array2D<PRISMATIC_FLOAT_PRECISION> testArr2D = zeros_ND<2,PRISMATIC_FLOAT_PRECISION>({{Rx,Ry}});

    Array2D<PRISMATIC_FLOAT_PRECISION> smallArr = fourierDownsample(testArr2D, Tx, Ty);
    BOOST_TEST(smallArr.get_dimi() == Tx);
    BOOST_TEST(smallArr.get_dimj() == Ty);
}

BOOST_AUTO_TEST_SUITE_END();

} //namespace Prismatic