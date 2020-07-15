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

BOOST_FIXTURE_TEST_CASE(CC_series_M, basicSim)
{
    meta.algorithm = Algorithm::Multislice;
    meta.simSeries = true;
    meta.seriesVals = {{-10.0, 0.0, 10.0}};
    meta.seriesKeys = {"probeDefocus"};
    meta.seriesTags = {"_df0000", "_df0001", "_df0002"};
    meta.filenameOutput = "../test/CC_series.h5";
    meta.save3DOutput = true;
    meta.save2DOutput = false;
    meta.save4DOutput = false;
    meta.savePotentialSlices = false;
    meta.saveDPC_CoM = false;

    divertOutput(pos, fd, logPath);
    std::cout << "\n######## BEGIN TEST CASE: CC_series_M ##########\n";
    go(meta);
    std::cout << "########## END TEST CASE: CC_series_M ##########\n";
    revertOutput(fd, pos);

    H5::H5File output = H5::H5File(meta.filenameOutput.c_str(), H5F_ACC_RDONLY);
    std::string basename = "virtual_detector_depth0000";
    H5::Group realslices = output.openGroup("4DSTEM_simulation/data/realslices");
    for(auto i = 0; i < meta.seriesTags.size(); i++)
    {
        std::string cur_name = basename + meta.seriesTags[i];
        bool nameCheck = realslices.nameExists(cur_name.c_str());
        BOOST_TEST(nameCheck);
    }
}

BOOST_AUTO_TEST_SUITE_END();

} //namespace Prismatic