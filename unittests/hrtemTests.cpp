#include "ioTests.h"
#include <boost/test/unit_test.hpp>
#include "ArrayND.h"
#include <iostream>
#include <vector>
#include "go.h"
#include "meta.h"
#include "params.h"
#include <stdio.h>
#include <random>
#include "fileIO.h"
#include "H5Cpp.h"

namespace Prismatic{

class basicSim{

    public:
    basicSim()     {setupSim(),BOOST_TEST_MESSAGE( "Setting up fixture");}
    ~basicSim()    {BOOST_TEST_MESSAGE( "Tearing down fixture");}
    Metadata<PRISMATIC_FLOAT_PRECISION> meta;
    Parameters<PRISMATIC_FLOAT_PRECISION> pars;
    std::string logPath = "ioTests.log";
    int fd;
    fpos_t pos;

    void setupSim()
    {
        //running from build directory
        meta.filenameAtoms = "../SI100.XYZ";
        meta.filenameOutput = "../test/fileIOtests.h5";
        meta.includeThermalEffects = 0;
        meta.save2DOutput = true;
        meta.save3DOutput = true;
        meta.save4DOutput = true;
        meta.saveDPC_CoM  = true;
        meta.savePotentialSlices = true;
        pars = Parameters<PRISMATIC_FLOAT_PRECISION>(meta);
    }
    

};

class logFile{
    public:
    logFile()       {setupLog(), BOOST_TEST_MESSAGE("Setting up ioTests.log file.");}
    ~logFile()      {BOOST_TEST_MESSAGE("Releasing ioTests.log file.");}
    std::string logPath;

    void setupLog()
    {
        logPath = "ioTests.log";
        FILE *fp = fopen(logPath.c_str(),"w");
        fprintf(fp,"########## BEGIN TEST SUITE: ioTests ##########\n");
        fclose(fp);
    }
};

void divertOutput(fpos_t &pos, int &fd, const std::string &file);
void revertOutput(const int &fd, fpos_t &pos);
void removeFile(const std::string &filepath);

BOOST_GLOBAL_FIXTURE(logFile);

BOOST_AUTO_TEST_SUITE(hrtemTests);

BOOST_FIXTURE_TEST_CASE(planeWave, basicSim)
{
    meta.algorithm = Algorithm::HRTEM;
    meta.filenameOutput = "../test/planeWave.h5";
    meta.filenameAtoms = "../test/au_np.xyz";
    meta.saveSMatrix = false;
    meta.savePotentialSlices = false;
    meta.saveComplexOutputWave = true;
    meta.potential3D = false;
    meta.maxXtilt = 0.3 / 1000;
    meta.maxYtilt = 0.3 / 1000;
    meta.xTiltOffset = 0.0 / 1000;
    meta.yTiltOffset = 0.0 / 1000;

    divertOutput(pos, fd, logPath);
    std::cout << "\n######### BEGIN TEST CASE: planeWave ##########\n";
    go(meta);
    std::cout << "########### END TEST CASE: planeWave ##########\n";
    revertOutput(fd, pos);


    H5::H5File testFile = H5::H5File(meta.filenameOutput.c_str(), H5F_ACC_RDONLY);
    H5::Group realslices = testFile.openGroup("4DSTEM_simulation/data/realslices");
    bool dataCheck = realslices.nameExists("HRTEM");
    BOOST_TEST(dataCheck);

    // check to see scaling is right
    double errSum = 0.0;
    double tol = 0.05; //5 % tolerance since can lose some electrons
    if(meta.saveComplexOutputWave)
    {
        Array3D<std::complex<PRISMATIC_FLOAT_PRECISION>> output;
        readComplexDataSet(output, meta.filenameOutput, "4DSTEM_simulation/data/realslices/HRTEM/realslice");
        for(auto &i : output) errSum += pow(std::abs(i), 2.0);
        errSum /= (PRISMATIC_FLOAT_PRECISION) output.size();
        BOOST_TEST(output.get_dimk() == 1);
    }
    else
    {
        Array3D<PRISMATIC_FLOAT_PRECISION> output = readDataSet3D(meta.filenameOutput, "4DSTEM_simulation/data/realslices/HRTEM/realslice");
        for(auto &i : output) errSum += i;
        errSum /= (PRISMATIC_FLOAT_PRECISION) output.size();
        BOOST_TEST(output.get_dimk() == 1);
    }
    BOOST_TEST(std::abs(1-errSum) < tol);

    removeFile(meta.filenameOutput);
}

BOOST_FIXTURE_TEST_CASE(imageTilts, basicSim)
{
    meta.algorithm = Algorithm::HRTEM;
    meta.filenameOutput = "../test/imageTilts.h5";
    meta.filenameAtoms = "../test/au_np.xyz";
    meta.saveSMatrix = false;
    meta.savePotentialSlices = false;
    meta.saveComplexOutputWave = false;
    meta.potential3D = false;
    meta.numGPUs = 1;
    meta.batchSizeCPU = 1;
    meta.realspacePixelSize[1] = 0.1;
    meta.realspacePixelSize[0] = 0.1;

    meta.minXtilt = 0 / 1000.0;
    meta.minYtilt = 0 / 1000.0;
    meta.maxXtilt = 20 / 1000.0; //mrads
    meta.maxYtilt = 20 / 1000.0;
    meta.xTiltOffset = 0.0 / 1000.0;
    meta.yTiltOffset = 0.0 / 1000.0;

    divertOutput(pos, fd, logPath);
    std::cout << "\n######### BEGIN TEST CASE: imageTilts #########\n";
    go(meta);
    std::cout << "########### END TEST CASE: imageTilts #########\n";
    revertOutput(fd, pos);

    // check to see scaling is right
    double errSum = 0.0;
    double tol = 0.05; //5 % tolerance since can lose some electrons
    if(meta.saveComplexOutputWave)
    {
        Array3D<std::complex<PRISMATIC_FLOAT_PRECISION>> output;
        readComplexDataSet(output, meta.filenameOutput, "4DSTEM_simulation/data/realslices/HRTEM/realslice");
        for(auto &i : output) errSum += pow(std::abs(i), 2.0);
        errSum /= (PRISMATIC_FLOAT_PRECISION) output.size();
        std::cout << "Number of tilts: " << output.get_dimk() << std::endl;
    }
    else
    {
        Array3D<PRISMATIC_FLOAT_PRECISION> output = readDataSet3D(meta.filenameOutput, "4DSTEM_simulation/data/realslices/HRTEM/realslice");
        for(auto &i : output) errSum += i;
        errSum /= (PRISMATIC_FLOAT_PRECISION) output.size();
        std::cout << "Number of tilts: " << output.get_dimk() << std::endl;
    }
    BOOST_TEST(std::abs(1-errSum) < tol);
    removeFile(meta.filenameOutput);
}

BOOST_FIXTURE_TEST_CASE(datacube, basicSim)
{
    meta.algorithm = Algorithm::HRTEM;
    meta.filenameOutput = "../test/hrtem_datacube.h5";
    meta.saveSMatrix = true;
    meta.savePotentialSlices = false;
    meta.saveComplexOutputWave = false;
    meta.potential3D = false;
    meta.maxXtilt = 2 / 1000.0; //mrads
    meta.maxYtilt = 2 / 1000.0;

    divertOutput(pos, fd, logPath);
    std::cout << "\n########## BEGIN TEST CASE: datacube ##########\n";
    go(meta);
    std::cout << "############ END TEST CASE: datacube ##########\n";
    revertOutput(fd, pos);

    // check to see that data exists in datacube location
    H5::H5File output = H5::H5File(meta.filenameOutput.c_str(), H5F_ACC_RDONLY);
    H5::Group datacubes = output.openGroup("4DSTEM_simulation/data/datacubes");
    BOOST_TEST(datacubes.nameExists("HRTEM"));

    // check to see that dim3 and dim4 are created properly; i.e., x and y are monotonically increasing
    std::array<size_t, 4> dim3 = {0,1,2,3};
    std::array<size_t, 4> dim4 = {0,2,1,3};
    bool dim3_check = true;
    bool dim4_check = true;
    for(auto i = 1; i < dim3.size(); i++)
    {
        if(dim3[i] < dim3[i-1]) dim3_check = false;
    }
    for(auto i = 1; i < dim4.size(); i++)
    {
        if(dim4[i] < dim4[i-1]) dim4_check = false;
    }
    
    BOOST_TEST(dim3_check);
    BOOST_TEST(dim4_check);
    //check to see N_dim3 x N_dim4 = N_beams in smatrix
    removeFile(meta.filenameOutput);
}

BOOST_AUTO_TEST_CASE(sortTest)
{
    //tests sorting along two dimensions for the purpose of ordering beams 
    int seed = 12345;
    srand(seed);
    std::default_random_engine de(seed);
    size_t N = 6;
    Array1D<PRISMATIC_FLOAT_PRECISION> xdata = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{N}});
    Array1D<PRISMATIC_FLOAT_PRECISION> ydata = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{N}});
    assignRandomValues(xdata,de);
    assignRandomValues(ydata,de);

    //we'll be working with vectors in the core code so copy data over into final pair form
    std::vector<std::pair<PRISMATIC_FLOAT_PRECISION, PRISMATIC_FLOAT_PRECISION>> tilts(N);
    for(auto i = 0; i < N; i++) tilts[i] = std::make_pair(xdata[i], ydata[i]);
    tilts[4].first = tilts[3].first;

    for(auto i = 0; i < N; i++) std::cout << tilts[i].first << " " << tilts[i].second << std::endl;
    std::cout << "---------------------------" << std::endl;

    //populate indices vector
    std::vector<size_t> indices(N);
    for(auto i = 0; i < N; i++) indices[i] = i;

    //sort the indices using tilt pairs as compartor
    std::sort(indices.begin(), indices.end(), [&](int i, int j){return tilts[i]<tilts[j];} );
    for(auto i = 0; i < N; i++) std::cout << tilts[i].first << " " << tilts[i].second << std::endl;

    std::cout << "---------------------------" << std::endl;
    for(auto i = 0; i < N; i++) std::cout << indices[i] << std::endl;

}

BOOST_AUTO_TEST_CASE(hdfEmpty)
{
    std::array<size_t, 4> test1 = {0,1,2,3};
    std::array<size_t, 4> test2 = {0,0,0,0};
    std::cout << &test1[1] << std::endl;
    std::cout << &test1[2] << std::endl;
    std::cout << &test1[3] << std::endl;
    std::cout << &test1[4] << std::endl;
    std::copy(&test1[0], &test1[4], &test2[0]);

    // removeFile(fname);
}
BOOST_AUTO_TEST_SUITE_END();

}