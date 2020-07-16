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
    meta.save2DOutput = true;
    meta.save4DOutput = true;
    meta.savePotentialSlices = false;
    meta.saveDPC_CoM = true;
    meta.probeStepX = 1;
    meta.probeStepY = 1;

    divertOutput(pos, fd, logPath);
    std::cout << "\n######## BEGIN TEST CASE: CC_series_M ##########\n";
    go(meta);
    std::cout << "########## END TEST CASE: CC_series_M ##########\n";
    revertOutput(fd, pos);

    H5::H5File output = H5::H5File(meta.filenameOutput.c_str(), H5F_ACC_RDONLY);
    std::string basename_3D = "virtual_detector_depth0000";
    std::string basename_2D = "annular_detector_depth0000";
    std::string basename_4D = "CBED_array_depth0000";
    std::string basename_DPC = "DPC_CoM_depth0000";

    H5::Group realslices = output.openGroup("4DSTEM_simulation/data/realslices");
    H5::Group datacubes = output.openGroup("4DSTEM_simulation/data/datacubes");

    for(auto i = 0; i < meta.seriesTags.size(); i++)
    {
        std::string cur_name_2D = basename_2D + meta.seriesTags[i];
        bool nameCheck_2D = realslices.nameExists(cur_name_2D.c_str());
        BOOST_TEST(nameCheck_2D);
        
        std::string cur_name_3D = basename_3D + meta.seriesTags[i];
        bool nameCheck_3D = realslices.nameExists(cur_name_3D.c_str());
        BOOST_TEST(nameCheck_3D);

        std::string cur_name_4D = basename_4D + meta.seriesTags[i];
        bool nameCheck_4D = datacubes.nameExists(cur_name_4D.c_str());
        BOOST_TEST(nameCheck_4D);
        
        std::string cur_name_DPC = basename_DPC + meta.seriesTags[i];
        bool nameCheck_DPC = realslices.nameExists(cur_name_DPC.c_str());
        BOOST_TEST(nameCheck_DPC);
    }
}

BOOST_AUTO_TEST_SUITE_END();

} //namespace Prismatic