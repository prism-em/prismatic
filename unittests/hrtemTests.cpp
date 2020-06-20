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
    meta.saveComplexOutputWave = false;
    meta.potential3D = false;

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
    PRISMATIC_FLOAT_PRECISION errSum = 0.0;
    PRISMATIC_FLOAT_PRECISION tol = 0.05; //5 % tolerance since can lose some electrons
    if(meta.saveComplexOutputWave)
    {
        Array3D<std::complex<PRISMATIC_FLOAT_PRECISION>> output;
        readComplexDataSet(output, meta.filenameOutput, "4DSTEM_simulation/data/realslices/HRTEM/realslice");
        for(auto &i : output) errSum += pow(std::abs(i), 2.0);
        errSum /= (PRISMATIC_FLOAT_PRECISION) output.size();
    }
    else
    {
        Array3D<PRISMATIC_FLOAT_PRECISION> output = readDataSet3D(meta.filenameOutput, "4DSTEM_simulation/data/realslices/HRTEM/realslice");
        for(auto &i : output) errSum += i;
        errSum /= (PRISMATIC_FLOAT_PRECISION) output.size();
    }
    BOOST_TEST(std::abs(1-errSum) < tol);
    removeFile(meta.filenameOutput);
}

BOOST_FIXTURE_TEST_CASE(imageTilts, basicSim)
{
    meta.algorithm = Algorithm::HRTEM;
    meta.filenameOutput = "../test/imageTilts.h5";
    meta.filenameAtoms = "../test/au_np.xyz";
    meta.saveSMatrix = true;
    meta.savePotentialSlices = true;
    meta.saveComplexOutputWave = false;
    meta.potential3D = false;
    meta.maxXtilt = 0.5 / 1000.0; //mrads
    meta.maxYtilt = 0.3 / 1000.0;
    meta.numGPUs = 0;
    meta.batchSizeCPU = 1;

    divertOutput(pos, fd, logPath);
    std::cout << "\n######### BEGIN TEST CASE: imageTilts #########\n";
    go(meta);
    std::cout << "########### END TEST CASE: imageTilts #########\n";
    revertOutput(fd, pos);

    // check to see scaling is right
    PRISMATIC_FLOAT_PRECISION errSum = 0.0;
    PRISMATIC_FLOAT_PRECISION tol = 0.05; //5 % tolerance since can lose some electrons
    if(meta.saveComplexOutputWave)
    {
        Array3D<std::complex<PRISMATIC_FLOAT_PRECISION>> output;
        readComplexDataSet(output, meta.filenameOutput, "4DSTEM_simulation/data/realslices/HRTEM/realslice");
        for(auto &i : output) errSum += pow(std::abs(i), 2.0);
        errSum /= (PRISMATIC_FLOAT_PRECISION) output.size();
    }
    else
    {
        Array3D<PRISMATIC_FLOAT_PRECISION> output = readDataSet3D(meta.filenameOutput, "4DSTEM_simulation/data/realslices/HRTEM/realslice");
        for(auto &i : output) errSum += i;
        errSum /= (PRISMATIC_FLOAT_PRECISION) output.size();
    }
    BOOST_TEST(std::abs(1-errSum) < tol);

    // std::vector<std::array<size_t, 3>> orders;
    // orders.push_back(std::array<size_t,3>{0,1,2});
    // orders.push_back(std::array<size_t,3>{0,2,1});
    // orders.push_back(std::array<size_t,3>{1,0,2});
    // orders.push_back(std::array<size_t,3>{1,2,0});
    // orders.push_back(std::array<size_t,3>{2,0,1});
    // orders.push_back(std::array<size_t,3>{2,1,0});
    // std::array<size_t, 3> dims = output.get_dimarr();
    // std::array<size_t, 3> cur_dims = {0,0,0};
    // Array3D<PRISMATIC_FLOAT_PRECISION> tmp_output = zeros_ND<3, PRISMATIC_FLOAT_PRECISION>(dims);
    // for(auto i = 0; i < 6; i++)
    // {
    //     for(auto j = 0; j < 6; j++)
    //     {
    //         cur_dims = {dims[orders[i][0]], dims[orders[i][1]], dims[orders[i][2]]};
    //         tmp_output = restride(output, cur_dims, orders[j]);
    //         std::string fname = "dump/"+std::to_string(i)+"_"+std::to_string(j)+".mrc";
    //         tmp_output.toMRC_f(fname.c_str());
    //     }
    // }
    // removeFile(meta.filenameOutput);
}

BOOST_AUTO_TEST_SUITE_END();

}