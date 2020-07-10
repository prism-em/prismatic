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
#include "utility.h"

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

BOOST_AUTO_TEST_SUITE(aberrationsTests);

BOOST_FIXTURE_TEST_CASE(probeIO_M, basicSim)
{
    meta.filenameOutput = "../test/probeIO_M.h5";
    meta.potential3D = false;
    meta.save2DOutput = false;
    meta.save3DOutput = true;
    meta.save4DOutput = false;
    meta.savePotentialSlices = false;
    meta.saveSMatrix = false;
    meta.saveDPC_CoM = false;
    meta.saveProbe = true;
    meta.algorithm = Algorithm::Multislice;

    divertOutput(pos, fd, logPath);
    std::cout << "\n######### BEGIN TEST CASE: probeIO_M ##########\n";

    go(meta);

    std::cout << "########### END TEST CASE: probeIO_M ##########\n";
    revertOutput(fd, pos);

    //check for array
    Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> probeArr;
    Array1D<PRISMATIC_FLOAT_PRECISION> dim1;
    Array1D<PRISMATIC_FLOAT_PRECISION> dim2;

    std::string datapath = "4DSTEM_simulation/data/diffractionslices/probe/data";
    std::string dim1path = "4DSTEM_simulation/data/diffractionslices/probe/dim1";
    std::string dim2path = "4DSTEM_simulation/data/diffractionslices/probe/dim2";
    std::vector<size_t> order2D = {0,1};
    std::vector<size_t> order1D = {0};
    bool dataCheck;
    try
    {
        readComplexDataSet(probeArr, meta.filenameOutput, datapath, order2D);
        dataCheck = true;
    }
    catch(...)
    {
        dataCheck = false;
    }

    bool dimCheck;
    try
    {
        readRealDataSet(dim1, meta.filenameOutput, dim1path, order1D);
        readRealDataSet(dim2, meta.filenameOutput, dim2path, order1D);
        dimCheck = true;
    }
    catch(...)
    {
        dimCheck = false;
    }

    BOOST_TEST(dataCheck);
    BOOST_TEST(dimCheck);
    BOOST_TEST(probeArr.get_dimi() = dim1.get_dimi());
    BOOST_TEST(probeArr.get_dimj() = dim2.get_dimi());
    removeFile(meta.filenameOutput);
};

BOOST_FIXTURE_TEST_CASE(probeIO_P, basicSim)
{
    meta.filenameOutput = "../test/probeIO_P.h5";
    meta.potential3D = false;
    meta.save2DOutput = false;
    meta.save3DOutput = true;
    meta.save4DOutput = false;
    meta.savePotentialSlices = false;
    meta.saveSMatrix = false;
    meta.saveDPC_CoM = false;
    meta.saveProbe = true;
    meta.algorithm = Algorithm::PRISM;

    divertOutput(pos, fd, logPath);
    std::cout << "\n######### BEGIN TEST CASE: probeIO_P ##########\n";

    go(meta);

    std::cout << "########### END TEST CASE: probeIO_P ##########\n";
    revertOutput(fd, pos);

    //check for array
    Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> probeArr;
    Array1D<PRISMATIC_FLOAT_PRECISION> dim1;
    Array1D<PRISMATIC_FLOAT_PRECISION> dim2;

    std::string datapath = "4DSTEM_simulation/data/diffractionslices/probe/data";
    std::string dim1path = "4DSTEM_simulation/data/diffractionslices/probe/dim1";
    std::string dim2path = "4DSTEM_simulation/data/diffractionslices/probe/dim2";
    std::vector<size_t> order2D = {0,1};
    std::vector<size_t> order1D = {0};
    bool dataCheck;
    try
    {
        readComplexDataSet(probeArr, meta.filenameOutput, datapath, order2D);
        dataCheck = true;
    }
    catch(...)
    {
        dataCheck = false;
    }

    bool dimCheck;
    try
    {
        readRealDataSet(dim1, meta.filenameOutput, dim1path, order1D);
        readRealDataSet(dim2, meta.filenameOutput, dim2path, order1D);
        dimCheck = true;
    }
    catch(...)
    {
        dimCheck = false;
    }

    BOOST_TEST(dataCheck);
    BOOST_TEST(dimCheck);
    BOOST_TEST(probeArr.get_dimi() = dim1.get_dimi());
    BOOST_TEST(probeArr.get_dimj() = dim2.get_dimi());
    removeFile(meta.filenameOutput);
};

BOOST_AUTO_TEST_SUITE_END();

}