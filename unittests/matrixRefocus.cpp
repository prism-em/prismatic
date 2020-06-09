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
    logFile()       {setupLog(), BOOST_TEST_MESSAGE("Setting up matrixRefocus.log file.");}
    ~logFile()      {BOOST_TEST_MESSAGE("Releasing matrixRefocus.log file.");}
    std::string logPath;

    void setupLog()
    {
        logPath = "matrixRefocus.log";
        FILE *fp = fopen(logPath.c_str(),"w");
        fprintf(fp,"####### BEGIN TEST SUITE: matrixRefocus #######\n");
        fclose(fp);
    }
};

void divertOutput(fpos_t &pos, int &fd, const std::string &file);
void revertOutput(const int &fd, fpos_t &pos);
void removeFile(const std::string &filepath);

BOOST_GLOBAL_FIXTURE(logFile);

BOOST_AUTO_TEST_SUITE(matrixRefocus);

BOOST_FIXTURE_TEST_CASE(refocus_test, basicSim)
{   

    std::string fname = "../test/matrixRefocus.h5";
    meta.filenameOutput = fname;
    // meta.matrixRefocus = false;
    meta.numGPUs = 1;
    divertOutput(pos, fd, logPath);
    std::cout << "\n####### BEGIN TEST CASE: refocus_test #########\n";

    go(meta);
    std::cout << "######### END TEST CASE: refocus_test #########\n";
    revertOutput(fd, pos);

    removeFile(fname);
}

BOOST_AUTO_TEST_SUITE_END();
}