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
#include "aberration.h"

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
    logFile()       {setupLog(), BOOST_TEST_MESSAGE("Setting up prismatic-tests.log file.");}
    ~logFile()      {BOOST_TEST_MESSAGE("Releasing prismatic-tests.log file.");}
    std::string logPath;

    void setupLog()
    {
        logPath = "prismatic-tests.log";
        FILE *fp = fopen(logPath.c_str(),"w");
        fprintf(fp,"########## BEGIN TEST SUITE: aberrationTests ##########\n");
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
    meta.filenameOutput = "../unittests/outputs/probeIO_M.h5";
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
    meta.filenameOutput = "../unittests/outputs/probeIO_P.h5";
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

BOOST_AUTO_TEST_CASE(parser)
{
    std::string fname = "../unittests/pfiles/ab_list";
    std::vector<aberration> abberations = readAberrations(fname);
    PRISMATIC_FLOAT_PRECISION tol = 0.0001;
    BOOST_TEST(abberations[0].m == 1);
    BOOST_TEST(abberations[0].n == 2);
    BOOST_TEST(std::abs(abberations[0].mag - 3.14) < tol);
    BOOST_TEST(std::abs(abberations[0].angle - 23) < tol);

};

BOOST_AUTO_TEST_CASE(astig)
{

    std::string fname = "../unittests/pfiles/astig";
    std::vector<aberration> abberations = readAberrations(fname);

    size_t imsize = 256;
    PRISMATIC_FLOAT_PRECISION pixelSize = 0.25;
    Array1D<PRISMATIC_FLOAT_PRECISION> qx = makeFourierCoords(imsize, pixelSize);
    Array1D<PRISMATIC_FLOAT_PRECISION> qy = makeFourierCoords(imsize, pixelSize);
    
    std::pair< Array2D<PRISMATIC_FLOAT_PRECISION>, Array2D<PRISMATIC_FLOAT_PRECISION> > mesh = meshgrid(qy,qx);
    Array2D<PRISMATIC_FLOAT_PRECISION> qya = mesh.first;
    Array2D<PRISMATIC_FLOAT_PRECISION> qxa = mesh.second;
    Array2D<PRISMATIC_FLOAT_PRECISION> q2(qya);
    std::transform(qxa.begin(), qxa.end(),
                qya.begin(), q2.begin(), [](const PRISMATIC_FLOAT_PRECISION& a, const PRISMATIC_FLOAT_PRECISION& b){
                return a*a + b*b;
            });
    Array2D<PRISMATIC_FLOAT_PRECISION> q1(q2);
    q2 = q2;
    q1 = q1;
    for (auto& q : q1)q=sqrt(q);

    Array2D<PRISMATIC_FLOAT_PRECISION> qTheta(q1);
    std::transform(qxa.begin(), qxa.end(),
                   qya.begin(), qTheta.begin(), [](const PRISMATIC_FLOAT_PRECISION& a, const PRISMATIC_FLOAT_PRECISION& b){
                    return atan2(b,a);
                });

    PRISMATIC_FLOAT_PRECISION lambda;
    constexpr double m = 9.109383e-31;
    constexpr double e = 1.602177e-19;
    constexpr double c = 299792458;
    constexpr double h = 6.62607e-34;
    PRISMATIC_FLOAT_PRECISION E0 = 80e3;
    lambda = (PRISMATIC_FLOAT_PRECISION) (h / sqrt(2 * m * e * E0) / sqrt(1 + e * E0 / 2 / m / c / c) * 1e10);

    Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> chi = getChi(q1, qTheta, lambda, abberations);
    for(auto i = 0; i < chi.size(); i++) chi[i] *= (q1[i] < 1.0) ? 1.0 : 0.0;
    H5::H5File testFile = H5::H5File("../unittests/outputs/astig.h5", H5F_ACC_TRUNC);
    H5::Group group = testFile.createGroup("test");
    hsize_t mdims[2] = {chi.get_dimi(), chi.get_dimj()};
    std::vector<size_t> order = {0,1};
    writeComplexDataSet(group, "data", &chi[0], mdims, 2, order);
    writeRealDataSet(group, "q1", &q1[0], mdims, 2, order);
    writeRealDataSet(group, "qtheta", &qTheta[0], mdims, 2, order);

};

BOOST_AUTO_TEST_CASE(abb1)
{

    std::string fname = "../unittests/pfiles/abb1";
    std::vector<aberration> abberations = readAberrations(fname);

    size_t imsize = 256;
    PRISMATIC_FLOAT_PRECISION pixelSize = 0.25;
    Array1D<PRISMATIC_FLOAT_PRECISION> qx = makeFourierCoords(imsize, pixelSize);
    Array1D<PRISMATIC_FLOAT_PRECISION> qy = makeFourierCoords(imsize, pixelSize);
    
    std::pair< Array2D<PRISMATIC_FLOAT_PRECISION>, Array2D<PRISMATIC_FLOAT_PRECISION> > mesh = meshgrid(qy,qx);
    Array2D<PRISMATIC_FLOAT_PRECISION> qya = mesh.first;
    Array2D<PRISMATIC_FLOAT_PRECISION> qxa = mesh.second;
    Array2D<PRISMATIC_FLOAT_PRECISION> q2(qya);
    std::transform(qxa.begin(), qxa.end(),
                qya.begin(), q2.begin(), [](const PRISMATIC_FLOAT_PRECISION& a, const PRISMATIC_FLOAT_PRECISION& b){
                return a*a + b*b;
            });
    Array2D<PRISMATIC_FLOAT_PRECISION> q1(q2);
    q2 = q2;
    q1 = q1;
    for (auto& q : q1)q=sqrt(q);

    Array2D<PRISMATIC_FLOAT_PRECISION> qTheta(q1);
    std::transform(qxa.begin(), qxa.end(),
                   qya.begin(), qTheta.begin(), [](const PRISMATIC_FLOAT_PRECISION& a, const PRISMATIC_FLOAT_PRECISION& b){
                    return atan2(b,a);
                });

    PRISMATIC_FLOAT_PRECISION lambda;
    constexpr double m = 9.109383e-31;
    constexpr double e = 1.602177e-19;
    constexpr double c = 299792458;
    constexpr double h = 6.62607e-34;
    PRISMATIC_FLOAT_PRECISION E0 = 80e3;
    lambda = (PRISMATIC_FLOAT_PRECISION) (h / sqrt(2 * m * e * E0) / sqrt(1 + e * E0 / 2 / m / c / c) * 1e10);

    Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> chi = getChi(q1, qTheta, lambda, abberations);
    for(auto i = 0; i < chi.size(); i++) chi[i] *= (q1[i] < 1.0) ? 1.0 : 0.0;
    H5::H5File testFile = H5::H5File("../unittests/outputs/abb1.h5", H5F_ACC_TRUNC);
    H5::Group group = testFile.createGroup("test");
    hsize_t mdims[2] = {chi.get_dimi(), chi.get_dimj()};
    std::vector<size_t> order = {0,1};
    writeComplexDataSet(group, "data", &chi[0], mdims, 2, order);
    writeRealDataSet(group, "q1", &q1[0], mdims, 2, order);
    writeRealDataSet(group, "qtheta", &qTheta[0], mdims, 2, order);
};

BOOST_AUTO_TEST_SUITE_END();

}