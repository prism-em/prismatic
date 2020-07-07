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

BOOST_AUTO_TEST_SUITE(probeTests);

BOOST_FIXTURE_TEST_CASE(rectGrid, basicSim)
{
    std::string refname = "../test/rectGridRef.h5";
    std::string testname = "../test/rectGridTest.h5";
    meta.filenameOutout = refname;
    divertOutput(pos, fd, logPath);
    std::cout << "\n######### BEGIN TEST CASE: planeWave ##########\n";
    go(meta);
    std::cout << "--------------------------------------------------\n";


    meta.probes_x = vecFromRange((PRISMATIC_FLOAT_PRECISION) 0.0, meta.probeStepX, (PRISMATIC_FLOAT_PRECISION) 0.99999*meta.cellDim[2]);
    meta.probes_y = vecFromRange((PRISMATIC_FLOAT_PRECISION) 0.0, meta.probeStepY, (PRISMATIC_FLOAT_PRECISION) 0.99999*meta.cellDim[1]);
    meta.filenameOutout = testname;
    meta.arbitraryProbes = true;
    go(meta);
    std::cout << "########### END TEST CASE: planeWave ##########\n";
    revertOutput(fd, pos);

    
    BOOST_TEST(1 == 0);

    removeFile(refname);
    removeFile(testname);
}

BOOST_AUTO_TEST_SUITE_END();

}