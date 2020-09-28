#include "potentialTests.h"
#include "projectedPotential.h"
#include "PRISM01_calcPotential.h"
#include <boost/test/unit_test.hpp>
#include "ArrayND.h"
#include <iostream>
#include "kirkland_params.h"
#include <vector>
#include "params.h"
#include "atom.h"
#include "go.h"

namespace Prismatic{

static const PRISMATIC_FLOAT_PRECISION pi = std::acos(-1);
PRISMATIC_FLOAT_PRECISION a0 = 0.529; //bohr radius
PRISMATIC_FLOAT_PRECISION e = 14.4; //electron charge in Volt-Angstoms
PRISMATIC_FLOAT_PRECISION term1 =  2*pi*pi*a0*e;
PRISMATIC_FLOAT_PRECISION term2 = 2*pow(pi,5.0/2.0)*a0*e;

class basicCell {

    public:
    basicCell()     {setupCell(),BOOST_TEST_MESSAGE( "Setting up fixture");}
    ~basicCell()    {BOOST_TEST_MESSAGE( "Tearing down fixture");}
    PRISMATIC_FLOAT_PRECISION cellDim = 20.0;
    PRISMATIC_FLOAT_PRECISION potBound = 2.0;
    PRISMATIC_FLOAT_PRECISION sliceThickness = 2.0;
    PRISMATIC_FLOAT_PRECISION sliceThicknessFactor = 4;
    PRISMATIC_FLOAT_PRECISION xPixel = 0.05;
    PRISMATIC_FLOAT_PRECISION yPixel = 0.05;
    PRISMATIC_FLOAT_PRECISION dzPot = sliceThickness/sliceThicknessFactor;
    Array1D<size_t> imageSize = zeros_ND<1, size_t>({{2}});
    std::vector<PRISMATIC_FLOAT_PRECISION> pixelSize;
    PRISMATIC_FLOAT_PRECISION yleng;
    PRISMATIC_FLOAT_PRECISION xleng;
    PRISMATIC_FLOAT_PRECISION zleng;
    Array1D<long> xvec;
    Array1D<long> yvec;
    Array1D<PRISMATIC_FLOAT_PRECISION> zvec;
    Array1D<PRISMATIC_FLOAT_PRECISION> xr;
    Array1D<PRISMATIC_FLOAT_PRECISION> yr;
    Array1D<PRISMATIC_FLOAT_PRECISION> zr;
    std::vector<atom> atoms;
    atom new_atom;
    std::vector<PRISMATIC_FLOAT_PRECISION> tiledCellDim;
    Parameters<PRISMATIC_FLOAT_PRECISION> pars;

    void setupCell()
    {
        imageSize[0] = (round(cellDim/(yPixel*16))*16);
        imageSize[1] = (round(cellDim/(xPixel*16))*16);
        pixelSize.push_back(cellDim/imageSize[0]);
        pixelSize.push_back(cellDim/imageSize[1]);
        pixelSize.push_back(dzPot);
        yleng = std::ceil(potBound/pixelSize[0]);
        xleng = std::ceil(potBound/pixelSize[1]);
        zleng = std::ceil(potBound/dzPot);
        xvec = zeros_ND<1,long>({{2*xleng+1}});
        yvec = zeros_ND<1,long>({{2*yleng+1}});
        zvec = zeros_ND<1,PRISMATIC_FLOAT_PRECISION>({{2*zleng}});
        xr = zeros_ND<1,PRISMATIC_FLOAT_PRECISION>({{2*xleng+1}});
        yr = zeros_ND<1,PRISMATIC_FLOAT_PRECISION>({{2*yleng+1}});
        zr = zeros_ND<1,PRISMATIC_FLOAT_PRECISION>({{2*zleng}});

        {
            PRISMATIC_FLOAT_PRECISION tmpx = -xleng;
            PRISMATIC_FLOAT_PRECISION tmpy = -yleng;
            for (auto &i : xvec)
                i = tmpx++;
            for (auto &j : yvec)
                j = tmpy++;
        }

        for (auto j = -zleng; j < zleng; j++)
		{
			zvec[j+zleng] = (PRISMATIC_FLOAT_PRECISION) j + 0.5;
		}

        for (auto i = 0; i < xr.size(); ++i) xr[i] = (PRISMATIC_FLOAT_PRECISION)xvec[i] * pixelSize[1];
        for (auto j = 0; j < yr.size(); ++j) yr[j] = (PRISMATIC_FLOAT_PRECISION)yvec[j] * pixelSize[0];
        for (auto j = 0; j < zr.size(); ++j) zr[j] = (PRISMATIC_FLOAT_PRECISION)zvec[j] * dzPot;

        // create list of atoms
        tiledCellDim.push_back(cellDim);
        tiledCellDim.push_back(cellDim);
        tiledCellDim.push_back(cellDim);
        pars.tiledCellDim = tiledCellDim;

        for(auto i = 0; i < 2; i++)
        {
            for(auto j = 0; j < 2; j++)
            {
                for(auto k = 0; k < 2; k++)
                {
                    new_atom.x = (1.0+2.0*i)/4.0;
                    new_atom.y = (1.0+2.0*j)/4.0;
                    new_atom.z = (1.0+2.0*k)/4.0;
                    new_atom.sigma = 0.1;
                    new_atom.species = 79;
                    new_atom.occ = 1;
                    atoms.push_back(new_atom);
                }
            }
        }

        // setting up parameter class
        pars.meta.sliceThickness = sliceThickness;
        pars.meta.zSampling = sliceThicknessFactor;
        pars.meta.randomSeed = 11111;
        pars.meta.numThreads = 8;
        pars.atoms = atoms;
        pars.pixelSize = pixelSize;
        pars.imageSize = imageSize;
        pars.dzPot = dzPot;
    }
};


//misc helper functions
PRISMATIC_FLOAT_PRECISION checkFaces(Array3D<PRISMATIC_FLOAT_PRECISION> &pot)
{
    
    size_t xFace[] = {0, pot.get_dimi()-1};
    size_t yFace[] = {0, pot.get_dimj()-1};
    size_t zFace[] = {0, pot.get_dimk()-1};
    PRISMATIC_FLOAT_PRECISION errSum = 0;
    //+-x
    for(auto i = 0; i < 2; i++)
    {
        for(auto j = 0; j < pot.get_dimj(); j++)
        {
            for(auto k = 0; k < pot.get_dimk(); k++)
            {
                errSum += pot.at(k,j,xFace[i]);
            }
        }
    }

    //+-y
    for(auto j = 0; j < 2; j++)
    {
        for(auto i = 0; i < pot.get_dimi(); i++)
        {
            for(auto k = 0; k < pot.get_dimk(); k++)
            {
                errSum += pot.at(k,yFace[j],i);
            }
        }
    }

    //+-z
    for(auto k = 0; k < 2; k++)
    {
        for(auto i = 0; i < pot.get_dimi(); i++)
        {
            for(auto j = 0; j < pot.get_dimj(); j++)
            {
                errSum += pot.at(zFace[k],j,i);
            }
        }
    }
    return errSum;
};

Array3D<PRISMATIC_FLOAT_PRECISION> generateInterpolationIdentity(size_t xSize,size_t ySize, size_t zSize,
                                                                            PRISMATIC_FLOAT_PRECISION dx,
                                                                            PRISMATIC_FLOAT_PRECISION dy,
                                                                            PRISMATIC_FLOAT_PRECISION dz)
{
    Array3D<PRISMATIC_FLOAT_PRECISION> identity = zeros_ND<3,PRISMATIC_FLOAT_PRECISION>({{zSize,ySize,xSize}});

    std::vector<size_t> xind = (dx < 0) ? std::vector<size_t>{0, identity.get_dimi()-2} : std::vector<size_t>{1, identity.get_dimi()-1};
    std::vector<size_t> yind = (dy < 0) ? std::vector<size_t>{0, identity.get_dimj()-2} : std::vector<size_t>{1, identity.get_dimj()-1};
    std::vector<size_t> zind = (dz < 0) ? std::vector<size_t>{0, identity.get_dimk()-2} : std::vector<size_t>{1, identity.get_dimk()-1};

    PRISMATIC_FLOAT_PRECISION dxm = (dx == 0) ? 1 : dx;
    PRISMATIC_FLOAT_PRECISION dym = (dy == 0) ? 1 : dy;
    PRISMATIC_FLOAT_PRECISION dzm = (dz == 0) ? 1 : dz;

    for(auto k = zind[0]; k <= zind[1]; k++)
    {
        for(auto j = yind[0]; j <= yind[1]; j++)
        {
            for(auto i = xind[0]; i <= xind[1]; i++)
            {
                identity.at(k,j,i) = 1;
                if(i == xind[0]) identity.at(k,j,i)*=std::abs(dxm);
                if(i == xind[1]) identity.at(k,j,i)*=std::abs(dx);

                if(j == yind[0]) identity.at(k,j,i)*=std::abs(dym);
                if(j == yind[1]) identity.at(k,j,i)*=std::abs(dy);
                
                if(k == zind[0]) identity.at(k,j,i)*=std::abs(dzm);
                if(k == zind[1]) identity.at(k,j,i)*=std::abs(dz);
            }
        }
    }

    return identity;
};

void printArray3D(Array3D<PRISMATIC_FLOAT_PRECISION> &arr)
{
    for(auto k = 0; k < arr.get_dimk(); k++)
    {
        for(auto j = 0; j < arr.get_dimj(); j++)
        {
            for(auto i = 0; i < arr.get_dimi(); i++)
            {
                std::cout<<arr.at(k,j,i) << " ";
            }
            std::cout<<std::endl;
        }
        std::cout << "---" << std::endl;
    }
    std::cout << "-----------------------" << std::endl;
};

void addMismatchedArray(Array3D<PRISMATIC_FLOAT_PRECISION> &big,
                        Array3D<PRISMATIC_FLOAT_PRECISION> &small,
                        const size_t &x_offset,
                        const size_t &y_offset,
                        const size_t &z_offset)
{//adds a smaller array to a region of larger area

    for(auto k = 0; k < small.get_dimk(); k++)
    {
        for(auto j = 0; j < small.get_dimj(); j++)
        {
            for(auto i = 0; i < small.get_dimi(); i++)
            {
                if(i+x_offset>= big.get_dimi()) std::cout <<"huhx" << std::endl;
                if(j+y_offset>= big.get_dimj()) std::cout <<"huhy" << std::endl;
                if(k+z_offset>= big.get_dimk()) std::cout <<"huhz" << std::endl;

                big.at(k+z_offset,j+y_offset,i+x_offset) += small.at(k,j,i);
            }
        }
    }

};

BOOST_AUTO_TEST_SUITE(potentialTests);

BOOST_AUTO_TEST_CASE(pot3DFunction)
{
    Array3D<PRISMATIC_FLOAT_PRECISION> radius = ones_ND<3,PRISMATIC_FLOAT_PRECISION>({{10,10,10}});
    Array1D<PRISMATIC_FLOAT_PRECISION> xr = ones_ND<1,PRISMATIC_FLOAT_PRECISION>({{10}});
    Array1D<PRISMATIC_FLOAT_PRECISION> yr = ones_ND<1,PRISMATIC_FLOAT_PRECISION>({{10}});
    Array1D<PRISMATIC_FLOAT_PRECISION> zr = ones_ND<1,PRISMATIC_FLOAT_PRECISION>({{10}});

    //seed corner with lower radius so that we can test function with potMin taken into acct
    xr.at(0) = 0.1;
    yr.at(0) = 0.1;
    zr.at(0) = 0.1;

    PRISMATIC_FLOAT_PRECISION r2 = pow(xr.at(0),2) + pow(yr.at(0),2) + pow(zr.at(0),2);
    PRISMATIC_FLOAT_PRECISION r = sqrt(r2);
    
    const size_t Z = 1;
    std::vector<PRISMATIC_FLOAT_PRECISION> parameters;
    parameters.resize(NUM_PARAMETERS);
    for (auto i =0; i < NUM_PARAMETERS; i++) parameters[i] = fparams[(Z-1)*NUM_PARAMETERS + i];
    
    PRISMATIC_FLOAT_PRECISION dz = zr.at(1) - zr.at(0);
    PRISMATIC_FLOAT_PRECISION pot_at_1 = term1*(parameters[0]*exp(-2*pi*sqrt(parameters[1]))
                                        + parameters[2]*exp(-2*pi*sqrt(parameters[3]))
                                        + parameters[4]*exp(-2*pi*sqrt(parameters[5])))
                                + term2*(parameters[6]*pow(parameters[7],-3.0/2.0)*exp(-pi*pi/parameters[7])
                                        + parameters[8]*pow(parameters[9],-3.0/2.0)*exp(-pi*pi/parameters[9])
                                        + parameters[10]*pow(parameters[11],-3.0/2.0)*exp(-pi*pi/parameters[11]));

    PRISMATIC_FLOAT_PRECISION pot_at_0 = term1*(parameters[0]*exp(-2*pi*r*sqrt(parameters[1]))/r
                                        + parameters[2]*exp(-2*pi*r*sqrt(parameters[3]))/r
                                        + parameters[4]*exp(-2*pi*r*sqrt(parameters[5]))/r)
                                + term2*(parameters[6]*pow(parameters[7],-3.0/2.0)*exp(-pi*pi*r2/parameters[7])
                                        + parameters[8]*pow(parameters[9],-3.0/2.0)*exp(-pi*pi*r2/parameters[9])
                                        + parameters[10]*pow(parameters[11],-3.0/2.0)*exp(-pi*pi*r2/parameters[11]));

    PRISMATIC_FLOAT_PRECISION expected = dz*(pot_at_0 - pot_at_1); //since potMin only checks 3 points
  
    Array3D<PRISMATIC_FLOAT_PRECISION> pot = kirklandPotential3D(Z, xr, yr, zr);
    PRISMATIC_FLOAT_PRECISION tol = 0.01;
    PRISMATIC_FLOAT_PRECISION error = std::abs(pot.at(0,0,0) - expected);

    BOOST_TEST(error<tol);
};

/*
BOOST_AUTO_TEST_CASE(potMin)
{
    //after potMin correction, potential should be non-negative with all faces of 3D potential prism = 0; 
    //test cubic
    PRISMATIC_FLOAT_PRECISION yleng = 20;
	PRISMATIC_FLOAT_PRECISION xleng = 20;
    PRISMATIC_FLOAT_PRECISION zleng = 20;

	ArrayND<1, std::vector<long>> xvec(std::vector<long>(2 * (size_t)xleng + 1, 0), {{2 * (size_t)xleng + 1}});
	ArrayND<1, std::vector<long>> yvec(std::vector<long>(2 * (size_t)yleng + 1, 0), {{2 * (size_t)yleng + 1}});
	ArrayND<1, std::vector<long>> zvec(std::vector<long>(2 * (size_t)zleng + 1, 0), {{2 * (size_t)zleng + 1}});
	{
		PRISMATIC_FLOAT_PRECISION tmpx = -xleng;
		PRISMATIC_FLOAT_PRECISION tmpy = -yleng;
		PRISMATIC_FLOAT_PRECISION tmpz = -zleng;
		for (auto &i : xvec)
			i = tmpx++;
		for (auto &j : yvec)
			j = tmpy++;
		for (auto &k : zvec)
			k = tmpz++;
	}

    PRISMATIC_FLOAT_PRECISION xPixel = 0.1;
    PRISMATIC_FLOAT_PRECISION yPixel = 0.1;
    PRISMATIC_FLOAT_PRECISION zPixel = 0.1;

    Array1D<PRISMATIC_FLOAT_PRECISION> xr(std::vector<PRISMATIC_FLOAT_PRECISION>(2 * (size_t)xleng + 1, 0), {{2 * (size_t)xleng + 1}});
	Array1D<PRISMATIC_FLOAT_PRECISION> yr(std::vector<PRISMATIC_FLOAT_PRECISION>(2 * (size_t)yleng + 1, 0), {{2 * (size_t)yleng + 1}});
	Array1D<PRISMATIC_FLOAT_PRECISION> zr(std::vector<PRISMATIC_FLOAT_PRECISION>(2 * (size_t)zleng + 1, 0), {{2 * (size_t)zleng + 1}});

	for (auto i = 0; i < xr.size(); ++i) xr[i] = (PRISMATIC_FLOAT_PRECISION)xvec[i] * xPixel;
	for (auto j = 0; j < yr.size(); ++j) yr[j] = (PRISMATIC_FLOAT_PRECISION)yvec[j] * yPixel;
	for (auto j = 0; j < zr.size(); ++j) zr[j] = (PRISMATIC_FLOAT_PRECISION)zvec[j] * zPixel;

    const size_t Z = 1; //Hydrogen!
    Array3D<PRISMATIC_FLOAT_PRECISION> pot = kirklandPotential3D(Z, xr, yr, zr);

    PRISMATIC_FLOAT_PRECISION tol = 0.0001;
    //test faces
    PRISMATIC_FLOAT_PRECISION errSum = checkFaces(pot);
    PRISMATIC_FLOAT_PRECISION minVal = pow(2,10); //check for nonnegativity
    for(auto i = 0; i < pot.size(); i++) minVal = (pot[i] < minVal) ? pot[i] : minVal;
    BOOST_TEST(errSum < tol);
    BOOST_TEST(minVal >= 0);

    //test rectangular in 1 direction
    xPixel = 0.05;
    for (auto i = 0; i < xr.size(); ++i) xr[i] = (PRISMATIC_FLOAT_PRECISION)xvec[i] * xPixel;
	for (auto j = 0; j < yr.size(); ++j) yr[j] = (PRISMATIC_FLOAT_PRECISION)yvec[j] * yPixel;
	for (auto j = 0; j < zr.size(); ++j) zr[j] = (PRISMATIC_FLOAT_PRECISION)zvec[j] * zPixel;
    pot = kirklandPotential3D(Z, xr, yr, zr);
    errSum = checkFaces(pot);
    minVal = pow(2,10);
    for(auto i = 0; i < pot.size(); i++) minVal = (pot[i] < minVal) ? pot[i] : minVal;
    BOOST_TEST(errSum < tol);
    BOOST_TEST(minVal >= 0);

    xPixel = 0.1;
    yPixel = 0.05;
	for (auto i = 0; i < xr.size(); ++i) xr[i] = (PRISMATIC_FLOAT_PRECISION)xvec[i] * xPixel;
	for (auto j = 0; j < yr.size(); ++j) yr[j] = (PRISMATIC_FLOAT_PRECISION)yvec[j] * yPixel;
	for (auto j = 0; j < zr.size(); ++j) zr[j] = (PRISMATIC_FLOAT_PRECISION)zvec[j] * zPixel;
    pot = kirklandPotential3D(Z, xr, yr, zr);
    errSum = checkFaces(pot);
    minVal = pow(2,10);
    for(auto i = 0; i < pot.size(); i++) minVal = (pot[i] < minVal) ? pot[i] : minVal;
    BOOST_TEST(errSum < tol);
    BOOST_TEST(minVal >= 0);

    yPixel = 0.1;
    zPixel = 0.05;
	for (auto i = 0; i < xr.size(); ++i) xr[i] = (PRISMATIC_FLOAT_PRECISION)xvec[i] * xPixel;
	for (auto j = 0; j < yr.size(); ++j) yr[j] = (PRISMATIC_FLOAT_PRECISION)yvec[j] * yPixel;
	for (auto j = 0; j < zr.size(); ++j) zr[j] = (PRISMATIC_FLOAT_PRECISION)zvec[j] * zPixel;
    pot = kirklandPotential3D(Z, xr, yr, zr);
    errSum = checkFaces(pot);
    minVal = pow(2,10);
    for(auto i = 0; i < pot.size(); i++) minVal = (pot[i] < minVal) ? pot[i] : minVal;
    BOOST_TEST(errSum < tol);
    BOOST_TEST(minVal >= 0);

    //test rectangulars in 2 directions
    xPixel = 0.05;
    yPixel = 0.05;
    zPixel = 0.1;
	for (auto i = 0; i < xr.size(); ++i) xr[i] = (PRISMATIC_FLOAT_PRECISION)xvec[i] * xPixel;
	for (auto j = 0; j < yr.size(); ++j) yr[j] = (PRISMATIC_FLOAT_PRECISION)yvec[j] * yPixel;
	for (auto j = 0; j < zr.size(); ++j) zr[j] = (PRISMATIC_FLOAT_PRECISION)zvec[j] * zPixel;
    pot = kirklandPotential3D(Z, xr, yr, zr);
    errSum = checkFaces(pot);
    minVal = pow(2,10);
    for(auto i = 0; i < pot.size(); i++) minVal = (pot[i] < minVal) ? pot[i] : minVal;
    BOOST_TEST(errSum < tol);
    BOOST_TEST(minVal >= 0);

    xPixel = 0.1;
    zPixel = 0.05;
	for (auto i = 0; i < xr.size(); ++i) xr[i] = (PRISMATIC_FLOAT_PRECISION)xvec[i] * xPixel;
	for (auto j = 0; j < yr.size(); ++j) yr[j] = (PRISMATIC_FLOAT_PRECISION)yvec[j] * yPixel;
	for (auto j = 0; j < zr.size(); ++j) zr[j] = (PRISMATIC_FLOAT_PRECISION)zvec[j] * zPixel;
    pot = kirklandPotential3D(Z, xr, yr, zr);
    errSum = checkFaces(pot);
    minVal = pow(2,10);
    for(auto i = 0; i < pot.size(); i++) minVal = (pot[i] < minVal) ? pot[i] : minVal;
    BOOST_TEST(errSum < tol);
    BOOST_TEST(minVal >= 0);

    yPixel = 0.1;
    xPixel = 0.05;
	for (auto i = 0; i < xr.size(); ++i) xr[i] = (PRISMATIC_FLOAT_PRECISION)xvec[i] * xPixel;
	for (auto j = 0; j < yr.size(); ++j) yr[j] = (PRISMATIC_FLOAT_PRECISION)yvec[j] * yPixel;
	for (auto j = 0; j < zr.size(); ++j) zr[j] = (PRISMATIC_FLOAT_PRECISION)zvec[j] * zPixel;
    pot = kirklandPotential3D(Z, xr, yr, zr);
    errSum = checkFaces(pot);
    minVal = pow(2,10);
    for(auto i = 0; i < pot.size(); i++) minVal = (pot[i] < minVal) ? pot[i] : minVal;
    BOOST_TEST(errSum < tol);
    BOOST_TEST(minVal >= 0);
};
*/

/*
BOOST_FIXTURE_TEST_CASE(projPotGeneration, basicCell)
{
    PRISMATIC_FLOAT_PRECISION tol = 0.001;

    //preparing reference potential for single atom potential
    const size_t Z = 79; //Gold!
    Array3D<PRISMATIC_FLOAT_PRECISION> refPot = zeros_ND<3,PRISMATIC_FLOAT_PRECISION>({{zr.size()*2,yr.size()*2,xr.size()*2}});
    Array3D<PRISMATIC_FLOAT_PRECISION> oneH = kirklandPotential3D(Z, xr, yr, zr);
    //loop through quadrants
    size_t offset_x;
    size_t offset_y;
    size_t offset_z;
    for(auto k = 0; k < 2; k++)
    {
        for(auto j = 0; j < 2; j++)
        {
            for(auto i = 0; i < 2; i++)
            {
                offset_x = i*oneH.get_dimi();
                offset_y = j*oneH.get_dimj();
                offset_z = k*oneH.get_dimk();
                addMismatchedArray(refPot, oneH, offset_x, offset_y, offset_z);
            }
        }
    }
    //expand dims to match lookup table
    Array4D<PRISMATIC_FLOAT_PRECISION> oneH_table = zeros_ND<4, PRISMATIC_FLOAT_PRECISION>({{1,oneH.get_dimk(),oneH.get_dimj(),oneH.get_dimi()}});
    for(auto k = 0; k < oneH.get_dimk(); k++)
    {
        for(auto j = 0; j < oneH.get_dimj(); j++)
        {
            for(auto i = 0; i < oneH.get_dimi(); i++)
            {
                oneH_table.at(0,k,j,i) = oneH.at(k,j,i);
            }
        }
    }

	std::vector<size_t> unique_species = get_unique_atomic_species(pars);
    generateProjectedPotentials3D(pars, oneH_table, unique_species, xvec, yvec, zvec); //H is the only lookup
    Array3D<PRISMATIC_FLOAT_PRECISION> testPot = pars.pot;
    PRISMATIC_FLOAT_PRECISION refPotSum = 0;
    PRISMATIC_FLOAT_PRECISION testPotSum = 0;
    PRISMATIC_FLOAT_PRECISION refPotSum2 = 0;
    for(auto i = 0; i < refPot.size(); i++) refPotSum += refPot[i];
    for(auto i = 0; i < oneH.size(); i++) refPotSum2 += 8*oneH[i];
    for(auto i = 0; i < testPot.size(); i++) testPotSum += testPot[i];

    // std::cout << "\% error (test vs ref): " << 100*std::abs(refPotSum-testPotSum)/refPotSum << std::endl;
    // std::cout << "\% error (test vs ref2): " << 100*std::abs(refPotSum2-testPotSum)/refPotSum2 << std::endl;
    // std::cout << "\% error (ref vs ref2): " << 100*std::abs(refPotSum-refPotSum2)/refPotSum2 << std::endl;
    BOOST_TEST(std::abs(refPotSum-testPotSum)/refPotSum<tol);
    BOOST_TEST(std::abs(refPotSum2-testPotSum)/refPotSum2<tol);

    //print min val in testPot
    PRISMATIC_FLOAT_PRECISION minVal = pow(2,10); //check for nonnegativity
    for(auto i = 0; i < testPot.size(); i++) minVal = (testPot[i] < minVal) ? testPot[i] : minVal;
    // std::cout << "minVal: " << minVal << std::endl;
    BOOST_TEST(minVal >=0);
    
    // testPot.toMRC_f("testPot.mrc");
    // refPot.toMRC_f("refPot.mrc");
};
*/

BOOST_FIXTURE_TEST_CASE(PRISM01_integration, basicCell)
{
    PRISMATIC_FLOAT_PRECISION tol = 0.001;
    //preparing reference potential for single atom potential
    const size_t Z = 79; //Gold!
    Array3D<PRISMATIC_FLOAT_PRECISION> refPot = zeros_ND<3,PRISMATIC_FLOAT_PRECISION>({{zr.size()*2,yr.size()*2,xr.size()*2}});
    Array3D<PRISMATIC_FLOAT_PRECISION> oneH = kirklandPotential3D(Z, xr, yr, zr);
    //loop through quadrants
    size_t offset_x;
    size_t offset_y;
    size_t offset_z;
    for(auto k = 0; k < 2; k++)
    {
        for(auto j = 0; j < 2; j++)
        {
            for(auto i = 0; i < 2; i++)
            {
                offset_x = i*oneH.get_dimi();
                offset_y = j*oneH.get_dimj();
                offset_z = k*oneH.get_dimk();
                addMismatchedArray(refPot, oneH, offset_x, offset_y, offset_z);
            }
        }
    }

    pars.meta.potential3D = true;
    pars.meta.numThreads = 1;
    PRISM01_calcPotential(pars);
    Array3D<PRISMATIC_FLOAT_PRECISION> testPot = pars.pot;

    PRISMATIC_FLOAT_PRECISION refPotSum = 0;
    PRISMATIC_FLOAT_PRECISION testPotSum = 0;
    PRISMATIC_FLOAT_PRECISION refPotSum2 = 0;
    for(auto i = 0; i < refPot.size(); i++) refPotSum += refPot[i];
    for(auto i = 0; i < oneH.size(); i++) refPotSum2 += 8*oneH[i];
    for(auto i = 0; i < testPot.size(); i++) testPotSum += testPot[i];

    BOOST_TEST(std::abs(refPotSum-testPotSum)/refPotSum<tol);
    BOOST_TEST(std::abs(refPotSum2-testPotSum)/refPotSum2<tol);

    std::cout << "min refpot1: " << *std::min_element(refPot.begin(), refPot.end()) << " max refpot1: " << *std::max_element(refPot.begin(), refPot.end()) << std::endl;
    std::cout << "min oneH: " << *std::min_element(oneH.begin(), oneH.end()) << " max oneH: " << *std::max_element(oneH.begin(), oneH.end()) << std::endl;
    std::cout << "min testPot: " << *std::min_element(testPot.begin(), testPot.end()) << " max testPot: " << *std::max_element(testPot.begin(), testPot.end()) << std::endl;
    //print min val in testPot
    PRISMATIC_FLOAT_PRECISION minVal = pow(2,10); //check for nonnegativity
    for(auto i = 0; i < testPot.size(); i++) minVal = (testPot[i] < minVal) ? testPot[i] : minVal;
    // std::cout << "minVal: " << minVal << std::endl;
    BOOST_TEST(minVal >=0);

    H5::H5File testFile = H5::H5File("../test/newpot3D.h5", H5F_ACC_TRUNC);
    hsize_t mdims[3] = {pars.pot.get_dimi(), pars.pot.get_dimj(), pars.pot.get_dimk()};
    H5::DataSpace mspace(3, mdims);
    H5::DataSet testPot_ds = testFile.createDataSet("pot_3D", PFP_TYPE, mspace);
    H5::DataSpace fspace = testPot_ds.getSpace();

    Array3D<PRISMATIC_FLOAT_PRECISION> strided_pot = zeros_ND<3,PRISMATIC_FLOAT_PRECISION>({{pars.pot.get_dimi(), pars.pot.get_dimj(), pars.pot.get_dimk()}});
    for(auto i = 0; i < pars.pot.get_dimi(); i++)
    {
        for(auto j = 0; j < pars.pot.get_dimj(); j++)
        {
            for(auto k = 0; k < pars.pot.get_dimk(); k++)
            {
                strided_pot.at(i,j,k) = pars.pot.at(k,j,i);
            }
        }
    }
    testPot_ds.write(&strided_pot[0], PFP_TYPE, mspace, fspace);

    fspace.close();
    mspace.close();
    testPot_ds.close();
    testFile.close();


};

BOOST_AUTO_TEST_CASE(pot_comparison)
{
    Metadata<PRISMATIC_FLOAT_PRECISION> meta;
    meta.filenameAtoms = "../unittests/pfiles/center_Au.xyz";
    // meta.filenameAtoms = "../SI100.XYZ";
    meta.filenameOutput = "../test/new_pot_ref.h5";
    meta.realspacePixelSize[0] = 0.05;
    meta.realspacePixelSize[1] = 0.05;
    meta.algorithm = Algorithm::Multislice;
    meta.probeStepX = 40;
    meta.probeStepY = 40;
    meta.savePotentialSlices = true;
    meta.potBound = 3.0;
    meta.numThreads = 1;
    meta.sliceThickness = 4.0;
    meta.potential3D = false;
    meta.includeThermalEffects = false;


    go(meta);

    std::cout << "#######################################################\n##############################################" << std::endl;
    meta.filenameOutput = "../test/new_pot_test.h5";
    meta.zSampling = 16;
    meta.potential3D = true;

    go(meta);
};

BOOST_AUTO_TEST_SUITE_END();

}