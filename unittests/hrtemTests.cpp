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

    meta.minXtilt = 104.8 / 1000.0;
    meta.minYtilt = 104.8 / 1000.0;
    meta.maxXtilt = 106 / 1000.0; //mrads
    meta.maxYtilt = 106 / 1000.0;
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
    // removeFile(meta.filenameOutput);
}

BOOST_AUTO_TEST_SUITE_END();

}