#include <boost/test/unit_test.hpp>
#include "ArrayND.h"
#include <iostream>
#include <vector>
#include <stdio.h>
#include <random>
#include "utility.h"
#include "fileIO.h"
#include "pprocess.h"

namespace Prismatic{

BOOST_AUTO_TEST_SUITE(processingTests);

BOOST_AUTO_TEST_CASE(poissonNoise)
{
    Array4D<PRISMATIC_FLOAT_PRECISION> testArr = zeros_ND<4,PRISMATIC_FLOAT_PRECISION>({{2,3,5,7}});
    for(auto i =0; i < testArr.size(); i++) testArr[i] = i+1;

    PRISMATIC_FLOAT_PRECISION scale = 1.0;
    applyPoisson(testArr, scale);
    for(auto i =0; i < 45; i++) std::cout << testArr[i] << std::endl;
    
}

BOOST_AUTO_TEST_SUITE_END();

} //namespace Prismatic