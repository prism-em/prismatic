#include "potentialTests.h"
#include "PRISM01_calcPotential.h"
#include <boost/test/unit_test.hpp>
#include "ArrayND.h"
#include <iostream>

namespace Prismatic{
    
BOOST_AUTO_TEST_SUITE(potentialTests);

BOOST_AUTO_TEST_CASE(kirklandPotential){
    Array3D<PRISMATIC_FLOAT_PRECISION> radius = ones_ND<3,PRISMATIC_FLOAT_PRECISION>({{10,10,10}});
    Array1D<PRISMATIC_FLOAT_PRECISION> factors = ones_ND<1, PRISMATIC_FLOAT_PRECISION>({{12}});
    PRISMATIC_FLOAT_PRECISION tol = 0.000001;
    PRISMATIC_FLOAT_PRECISION error_sum = 0;
    //test first term
    for(auto i = 6; i < 12; i+=2) factors.at(i) = 0;
    Array3D<PRISMATIC_FLOAT_PRECISION> val = kirklandPotential3D(factors, radius);
    PRISMATIC_FLOAT_PRECISION expected = 0.8423963; //all terms = 1 at radius 1
    for(auto z = 0; z < radius.get_dimk(); z++){
        for(auto y = 0; y < radius.get_dimj(); y++ ){
            for(auto x = 0; x < radius.get_dimi(); x++){
                error_sum += std::abs(val.at(z,y,x)-expected);
            }
        }
    }
    error_sum/=1000; //take mean
    BOOST_TEST(error_sum<tol);

    //test second term
    for(auto i = 0; i < 12; i++) factors.at(i) = 1;
    for(auto i = 0; i < 6; i++) factors.at(i) = 0;
    val = kirklandPotential3D(factors, radius);
    expected = 0.041355;
    error_sum = 0.0;
    for(auto z = 0; z < radius.get_dimk(); z++){
        for(auto y = 0; y < radius.get_dimj(); y++ ){
            for(auto x = 0; x < radius.get_dimi(); x++){
                error_sum += std::abs(val.at(z,y,x)-expected);
            }
        }
    }
    error_sum/=1000; //take mean
    BOOST_TEST(error_sum<tol);

    //test both
    for(auto i = 0; i < 12; i++) factors.at(i) = 1;
    val = kirklandPotential3D(factors, radius);
    expected = 0.883751;
    error_sum = 0.0;
    for(auto z = 0; z < radius.get_dimk(); z++){
        for(auto y = 0; y < radius.get_dimj(); y++ ){
            for(auto x = 0; x < radius.get_dimi(); x++){
                error_sum += std::abs(val.at(z,y,x)-expected);
            }
        }
    }
    error_sum/=1000; //take mean
    BOOST_TEST(error_sum<tol);
};

BOOST_AUTO_TEST_SUITE_END();

}