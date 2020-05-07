#include "potentialTests.h"
#include "projectedPotential.h"
#include <boost/test/unit_test.hpp>
#include "ArrayND.h"
#include <iostream>
#include "kirkland_params.h"

namespace Prismatic{

    static const PRISMATIC_FLOAT_PRECISION pi = std::acos(-1);
	PRISMATIC_FLOAT_PRECISION a0 = 0.529; //bohr radius
	PRISMATIC_FLOAT_PRECISION e = 14.4; //electron charge in Volt-Angstoms
	PRISMATIC_FLOAT_PRECISION term1 =  2*pi*pi*a0*e;
	PRISMATIC_FLOAT_PRECISION term2 = 2*pow(pi,5.0/2.0)*a0*e;
    
BOOST_AUTO_TEST_SUITE(potentialTests);

BOOST_AUTO_TEST_CASE(potential3D){
    Array3D<PRISMATIC_FLOAT_PRECISION> radius = ones_ND<3,PRISMATIC_FLOAT_PRECISION>({{10,10,10}});
    Array1D<PRISMATIC_FLOAT_PRECISION> xr = ones_ND<1,PRISMATIC_FLOAT_PRECISION>({{10}});
    Array1D<PRISMATIC_FLOAT_PRECISION> yr = ones_ND<1,PRISMATIC_FLOAT_PRECISION>({{10}});
    Array1D<PRISMATIC_FLOAT_PRECISION> zr = ones_ND<1,PRISMATIC_FLOAT_PRECISION>({{10}});

    //normalize radius to 1 at all points in grid
    xr /= sqrt(3);
    yr /= sqrt(3);
    zr /= sqrt(3);
    
    const size_t Z = 1;
    std::vector<PRISMATIC_FLOAT_PRECISION> parameters;
    parameters.resize(NUM_PARAMETERS);
    for (auto i =0; i < NUM_PARAMETERS; i++) parameters[i] = fparams[(Z-1)*NUM_PARAMETERS + i];
    
    PRISMATIC_FLOAT_PRECISION expected = term1*(parameters[0]*exp(-2*pi*sqrt(parameters[1]))
                                        + parameters[2]*exp(-2*pi*sqrt(parameters[3]))
                                        + parameters[4]*exp(-2*pi*sqrt(parameters[5])))
                                + term2*(parameters[6]*pow(parameters[7],-3.0/2.0)*exp(-pi*pi/parameters[7])
                                        + parameters[8]*pow(parameters[9],-3.0/2.0)*exp(-pi*pi/parameters[9])
                                        + parameters[10]*pow(parameters[11],-3.0/2.0)*exp(-pi*pi/parameters[11])); 
  
    Array3D<PRISMATIC_FLOAT_PRECISION> pot = kirklandPotential3D(Z, xr, yr, zr);
    PRISMATIC_FLOAT_PRECISION tol = 0.000001;
    PRISMATIC_FLOAT_PRECISION error_sum = 0;
    for(auto z = 0; z < radius.get_dimk(); z++){
        for(auto y = 0; y < radius.get_dimj(); y++ ){
            for(auto x = 0; x < radius.get_dimi(); x++){
                error_sum += std::abs(pot.at(z,y,x)-expected);
            }
        }
    }
    error_sum/=1000; //take mean
    BOOST_TEST(error_sum<tol);
};

BOOST_AUTO_TEST_SUITE_END();

}