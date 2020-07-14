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
#include "utility.h"
#include "probe.h"

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

BOOST_AUTO_TEST_SUITE(seriesTests);

BOOST_AUTO_TEST_CASE(first)
{
    BOOST_TEST(1==0);
}

BOOST_AUTO_TEST_SUITE_END();

} //namespace Prismatic