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
#include "fftw3.h"
#include "ioTests.h"
#include "utility.h"

namespace Prismatic{

class basicSim{

    public:
    basicSim()     {setupSim(),BOOST_TEST_MESSAGE( "Setting up fixture");}
    ~basicSim()    {BOOST_TEST_MESSAGE( "Tearing down fixture");}
    Metadata<PRISMATIC_FLOAT_PRECISION> meta;
    Parameters<PRISMATIC_FLOAT_PRECISION> pars;
    std::string logPath = "prismatic-tests.log";
    int fd;
    fpos_t pos;

    void setupSim()
    {
        //running from build directory
        meta.filenameAtoms = "../SI100.XYZ";
        meta.filenameOutput = "../unittests/outputs/fileIOtests.h5";
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
        fprintf(fp,"####### BEGIN TEST SUITE: refocusTests #######\n");
        fclose(fp);
    }
};

void divertOutput(fpos_t &pos, int &fd, const std::string &file);
void revertOutput(const int &fd, fpos_t &pos);
void removeFile(const std::string &filepath);


BOOST_GLOBAL_FIXTURE(logFile);

BOOST_AUTO_TEST_SUITE(refocusTests);

BOOST_FIXTURE_TEST_CASE(matrixRefocus, basicSim)
{   
    /*
    Test case is designed as follows.
    Matrix refocusing is only relevant for the PRISM method; i.e, there is no error to correct for in multislice
    Therefore, ideally we can compare to a multislice simulation.

    tile the Z dimension for the SI example cell 10-20 times, so it is decently thick and refocusing matters

    In multislice, we choose three defocii - C1 = -100, 0, 100 and simulate the complex output wave 4D output
    apply a backpropagation operator to each image so that probes can compared on equivalent bases

    In PRISM, simulate the same defocii with f=1 refocusing the matrix (and using the C1 abberation in the probe formation)

    Compare values of complex output probes and ensure equality within tolerance

    Could be expanded to f > 1 but would need carefulc ropping and downsampling for proper checking


    */

    std::string fname_m = "../unittests/outputs/matrixRefocus_m.h5";
    std::string fname_p = "../unittests/outputs/matrixRefocus_p.h5";
    meta.tileX = 5;
    meta.tileY = 5;
    meta.tileZ = 1;
    meta.filenameOutput = fname_m;
    meta.probeDefocus = 10.0;
    meta.matrixRefocus = false;
    meta.algorithm = Algorithm::Multislice;
    meta.saveComplexOutputWave = false;
    meta.savePotentialSlices = true;
    meta.probeStepX = 10;
    meta.probeStepY = 10;
    meta.alphaBeamMax = 50 / 1000.0;
    meta.potBound = 3.0;
    meta.numGPUs = 1;
    meta.numThreads = 12;
    meta.realspacePixelSize[0] = 0.1;
    meta.realspacePixelSize[1] = 0.1;

    divertOutput(pos, fd, logPath);
    std::cout << "\n####### BEGIN TEST CASE: refocus_test #########\n";
    go(meta);

    std::cout << "\n--------------------------------------------\n";
    meta.matrixRefocus = true;
    meta.algorithm = Algorithm::PRISM;
    meta.interpolationFactorX = 1;
    meta.interpolationFactorY = 1;
    meta.filenameOutput = fname_p;
    go(meta);

    std::cout << "######### END TEST CASE: refocus_test #########\n";
    revertOutput(fd, pos);


    //read and compare probes    
    std::string dataPath4D = "4DSTEM_simulation/data/datacubes/CBED_array_depth0000/data";
    Array4D<std::complex<PRISMATIC_FLOAT_PRECISION>> refProbes;
    Array4D<std::complex<PRISMATIC_FLOAT_PRECISION>> testProbes;
    readComplexDataSet_inOrder(refProbes, fname_m, dataPath4D);
    readComplexDataSet_inOrder(testProbes, fname_p, dataPath4D);

    PRISMATIC_FLOAT_PRECISION tol = 1e-7;
    BOOST_TEST(compareSize(refProbes, testProbes));
    if(compareSize(refProbes, testProbes))
    {
        PRISMATIC_FLOAT_PRECISION error = 0.0;
        for(auto i = 0; i < refProbes.size(); i++)
        {
            error += std::abs(pow(std::abs(refProbes[i]), 2.0) - pow(std::abs(testProbes[i]),2.0));
        }
        PRISMATIC_FLOAT_PRECISION meanError = error / refProbes.size();
        BOOST_TEST(meanError < tol);
    }
    else
    {
        std::array<size_t, 4> refDims = refProbes.get_dimarr();
        std::array<size_t, 4> testDims = testProbes.get_dimarr();
        for(auto i = 0; i < 4; i++)
        {
            std::cout << refDims[i] << " " << testDims[i] << std::endl;
        }
    }


    removeFile(fname_m);
    removeFile(fname_p);
}

BOOST_AUTO_TEST_CASE(boolstream)
{
    bool check = true;
    double two = 2.4;
    char test = two;
}

BOOST_AUTO_TEST_SUITE_END();
}