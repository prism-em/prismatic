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

BOOST_AUTO_TEST_SUITE(probeTests);

BOOST_FIXTURE_TEST_CASE(rectGrid_M, basicSim)
{
    std::string refname = "../test/rectGridRef.h5";
    std::string testname = "../test/rectGridTest.h5";
    meta.filenameOutput = refname;
    meta.algorithm = Algorithm::Multislice;
    divertOutput(pos, fd, logPath);
    std::cout << "\n######## BEGIN TEST CASE: rectGrid_M ##########\n";
    go(meta);
    std::cout << "--------------------------------------------------\n";

    //prepare identical rectangular grid as lisdt
    PRISMATIC_FLOAT_PRECISION r0 = 0.0;
    PRISMATIC_FLOAT_PRECISION r1 = 0.99999*5.43;
    std::vector<PRISMATIC_FLOAT_PRECISION> px = vecFromRange(r0, meta.probeStepX, r1);
    std::vector<PRISMATIC_FLOAT_PRECISION> py = vecFromRange(r0, meta.probeStepY, r1);
    meta.probes_x = {};
    meta.probes_y = {};
    for(auto i = 0; i < px.size(); i++)
    {
        for(auto j = 0; j < py.size(); j++)
        {
            meta.probes_x.push_back(px[i]);
            meta.probes_y.push_back(py[j]);
        }
    }
    meta.filenameOutput = testname;
    meta.arbitraryProbes = true;
    go(meta);
    std::cout << "########## END TEST CASE: rectGrid_M ##########\n";
    revertOutput(fd, pos);

    std::string datapath = "4DSTEM_simulation/data/realslices/virtual_detector_depth0000/realslice";
    Array3D<PRISMATIC_FLOAT_PRECISION> testArr;
    Array3D<PRISMATIC_FLOAT_PRECISION> refArr;
    std::vector<size_t> order = {0,1,2};
    readRealDataSet(refArr, refname, datapath, order);
    readRealDataSet(testArr, testname, datapath, order);

    //check for proper reshaping of arrays
    BOOST_TEST(testArr.get_dimj() == 1);
    BOOST_TEST(testArr.get_dimi() == (refArr.get_dimi()*refArr.get_dimj()));

    //since equivalent probe positions were run, check total error
    PRISMATIC_FLOAT_PRECISION err = 0.0;
    PRISMATIC_FLOAT_PRECISION tol = 0.0001;

    for(auto i = 0; i < testArr.get_dimi(); i++)
    {
        size_t ay = i % refArr.get_dimi();
        size_t ax = i / refArr.get_dimi();
        for(auto b = 0; b < testArr.get_dimk(); b++)
        {
            err += std::abs(testArr.at(b,0,i) - refArr.at(b,ay,ax));
        }
    }

    BOOST_TEST(err < tol);
    removeFile(refname);
    removeFile(testname);
}

BOOST_FIXTURE_TEST_CASE(rectGrid_P, basicSim)
{
    std::string refname = "../test/rectGridRef.h5";
    std::string testname = "../test/rectGridTest.h5";
    meta.filenameOutput = refname;
    meta.algorithm = Algorithm::PRISM;
    divertOutput(pos, fd, logPath);
    std::cout << "\n######## BEGIN TEST CASE: rectGrid_P ##########\n";
    go(meta);
    std::cout << "--------------------------------------------------\n";

    //prepare identical rectangular grid as lisdt
    PRISMATIC_FLOAT_PRECISION r0 = 0.0;
    PRISMATIC_FLOAT_PRECISION r1 = 0.99999*5.43;
    std::vector<PRISMATIC_FLOAT_PRECISION> px = vecFromRange(r0, meta.probeStepX, r1);
    std::vector<PRISMATIC_FLOAT_PRECISION> py = vecFromRange(r0, meta.probeStepY, r1);
    meta.probes_x = {};
    meta.probes_y = {};
    for(auto i = 0; i < px.size(); i++)
    {
        for(auto j = 0; j < py.size(); j++)
        {
            meta.probes_x.push_back(px[i]);
            meta.probes_y.push_back(py[j]);
        }
    }
    meta.filenameOutput = testname;
    meta.arbitraryProbes = true;
    go(meta);
    std::cout << "########## END TEST CASE: rectGrid_P ##########\n";
    revertOutput(fd, pos);

    std::string datapath = "4DSTEM_simulation/data/realslices/virtual_detector_depth0000/realslice";
    Array3D<PRISMATIC_FLOAT_PRECISION> testArr;
    Array3D<PRISMATIC_FLOAT_PRECISION> refArr;
    std::vector<size_t> order = {0,1,2};
    readRealDataSet(refArr, refname, datapath, order);
    readRealDataSet(testArr, testname, datapath, order);

    //check for proper reshaping of arrays
    BOOST_TEST(testArr.get_dimj() == 1);
    BOOST_TEST(testArr.get_dimi() == (refArr.get_dimi()*refArr.get_dimj()));

    //since equivalent probe positions were run, check total error
    PRISMATIC_FLOAT_PRECISION err = 0.0;
    PRISMATIC_FLOAT_PRECISION tol = 0.0001;

    for(auto i = 0; i < testArr.get_dimi(); i++)
    {
        size_t ay = i % refArr.get_dimi();
        size_t ax = i / refArr.get_dimi();
        for(auto b = 0; b < testArr.get_dimk(); b++)
        {
            err += std::abs(testArr.at(b,0,i) - refArr.at(b,ay,ax));
        }
    }

    BOOST_TEST(err < tol);
    removeFile(refname);
    removeFile(testname);
}

BOOST_AUTO_TEST_CASE(parser)
{

    std::vector<PRISMATIC_FLOAT_PRECISION> xprobes;
    std::vector<PRISMATIC_FLOAT_PRECISION> yprobes;
    std::string fname = "../unittests/pfiles/probes_real";

    std::vector<PRISMATIC_FLOAT_PRECISION> xcheck = {0.0, 1.0, 2.0, 3.145, 7.6};
    std::vector<PRISMATIC_FLOAT_PRECISION> ycheck = {5.0, 3.145, 7.6, 1.0, 0.32};
    readProbes(fname, xprobes, yprobes);
    BOOST_TEST(xprobes == xcheck);
    BOOST_TEST(yprobes == ycheck);

}

BOOST_AUTO_TEST_SUITE_END();

} //namespace Prismatic