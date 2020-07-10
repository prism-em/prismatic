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
#include <thread>

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

void divertOutput(fpos_t &pos, int &fd, const std::string &file)
{
    std::cout << "Opening log file " << file << " to capture simulation output." << std::endl;
    fgetpos(stdout, &pos);
    fd = dup(fileno(stdout));
    freopen(file.c_str(),"a",stdout);
};

void revertOutput(const int &fd, fpos_t &pos)
{
    //clean up files
    fflush(stdout);
    dup2(fd,fileno(stdout));
    close(fd);
    clearerr(stdout);
    fsetpos(stdout, &pos);
    std::cout << "Log file closed. Output returning to terminal." << std::endl;
};

void removeFile(const std::string &filepath)
{
    if( remove( filepath.c_str() ) != 0 )
        perror( "Error deleting file" );
    else
        puts( "Test file successfully deleted" );
};

void testFunc(void *elem)
{
    int* ip = (int *) elem;
    *ip += *ip;
};

void intFunc(int x, int y)
{
    std::cout << x+y << std::endl;
};

herr_t CBED_process(void *elem, hid_t type_id, unsigned ndim, const hsize_t *point, void *operator_data)
{
    try
    {
        PRISMATIC_FLOAT_PRECISION* tmp_elem = (PRISMATIC_FLOAT_PRECISION *) elem;
        PRISMATIC_FLOAT_PRECISION* tmp_data = (PRISMATIC_FLOAT_PRECISION *) operator_data;
        *tmp_elem += *tmp_data;
        operator_data++;
        return 0;
    }
    catch(...)
    {
        //if fails
        return -1;
    }

}

BOOST_GLOBAL_FIXTURE(logFile);

BOOST_AUTO_TEST_SUITE(ioTests);

BOOST_FIXTURE_TEST_CASE(operationReorganization, basicSim)
{
    //make sure nothing is broken by moving all file IO operations to their own source
    divertOutput(pos, fd, logPath);

    std::cout << "\n## BEGIN TEST CASE: operationReorganization ###\n";

    go(meta);

    std::cout << "\n--------------------------------------------\n";
    meta.algorithm = Algorithm::Multislice;
    go(meta);
    std::cout << "### END TEST CASE: operationReorganization ####\n";

    revertOutput(fd, pos);
    
    removeFile(meta.filenameOutput);

};

BOOST_FIXTURE_TEST_CASE(readH5, basicSim)
{
    //setting up test arrays
    int seed = 10101;
    srand(seed);
    std::default_random_engine de(seed);

    Array2D<PRISMATIC_FLOAT_PRECISION> testArr2D = zeros_ND<2,PRISMATIC_FLOAT_PRECISION>({{2,7}});
    Array3D<PRISMATIC_FLOAT_PRECISION> testArr3D = zeros_ND<3,PRISMATIC_FLOAT_PRECISION>({{2,7,5}});
    Array4D<PRISMATIC_FLOAT_PRECISION> testArr4D = zeros_ND<4,PRISMATIC_FLOAT_PRECISION>({{2,7,5,3}});
    Array2D<PRISMATIC_FLOAT_PRECISION> readBuffer = zeros_ND<2,PRISMATIC_FLOAT_PRECISION>({{5,3}});
    assignRandomValues(testArr2D, de);
    assignRandomValues(testArr3D, de);
    assignRandomValues(testArr4D, de);

    //setting up output file to read from
    pars.meta.filenameOutput = "../test/fileIOreadTest.h5";
    pars.outputFile = H5::H5File(pars.meta.filenameOutput.c_str(), H5F_ACC_TRUNC);
    setupOutputFile(pars);

    //create datasets
	H5::Group datacubes = pars.outputFile.openGroup("4DSTEM_simulation/data/datacubes");
    
    //lower scope to destroy variables
    {
        std::string base_name = "test4Darray";
        H5::Group test4D_group(datacubes.createGroup(base_name.c_str()));
        hsize_t data_dims[4];
        data_dims[0] = {testArr4D.get_dimk()};
        data_dims[1] = {testArr4D.get_diml()};
        data_dims[2] = {testArr4D.get_dimi()};
        data_dims[3] = {testArr4D.get_dimj()};

        //create dataset
        H5::DataSpace mspace(4, data_dims); //rank is 4
        H5::DataSet test4D_data = test4D_group.createDataSet("datacube", PFP_TYPE, mspace);
        mspace.close();
        test4D_group.close();
    }
    datacubes.close();

    //write 4D array to HDF5 file
    std::string nameString = "4DSTEM_simulation/data/datacubes/test4Darray";
    PRISMATIC_FLOAT_PRECISION numFP = 1;

    {
        hsize_t mdims[4];
        mdims[0] = {1};
        mdims[1] = {1};
        mdims[2] = {testArr4D.get_dimi()};
        mdims[3] = {testArr4D.get_dimj()};
        hsize_t offset[4] = {0,0,0,0};
        for(auto l = 0; l < testArr4D.get_diml(); l++)
        {
            for(auto k = 0; k < testArr4D.get_dimk(); k++)
            {
                offset[0] = k;
                offset[1] = l;
                size_t arrayStart = l*testArr4D.get_dimk()*testArr4D.get_dimj()*testArr4D.get_dimi()+k*testArr4D.get_dimj()*testArr4D.get_dimi();
                writeDatacube4D(pars,&testArr4D[arrayStart],&readBuffer[0],mdims,offset,numFP,nameString.c_str());

            }
        }
    }

    H5::Group realslices = pars.outputFile.openGroup("4DSTEM_simulation/data/realslices");
    {
        std::string base_name = "test3Darray";
        H5::Group test3D_group(realslices.createGroup(base_name.c_str()));
        hsize_t data_dims[3];
        data_dims[0] = {testArr3D.get_dimi()};
        data_dims[1] = {testArr3D.get_dimj()};
        data_dims[2] = {testArr3D.get_dimk()};

        //create dataset
        H5::DataSpace mspace(3, data_dims);
        H5::DataSet test3D_data = test3D_group.createDataSet("realslice", PFP_TYPE, mspace);
        mspace.close();
        test3D_group.close();
    }

    {
        std::string base_name = "test2Darray";
        H5::Group test2D_group(realslices.createGroup(base_name.c_str()));
        hsize_t data_dims[2];
        data_dims[0] = {testArr2D.get_dimi()};
        data_dims[1] = {testArr2D.get_dimj()};

        //create dataset
        H5::DataSpace mspace(2, data_dims);
        H5::DataSet test2D_data = test2D_group.createDataSet("realslice", PFP_TYPE, mspace);
        mspace.close();
        test2D_group.close();
    }

    //write 2D and 3D arrays
    {
        H5::DataSet test3D_data = realslices.openDataSet("test3Darray/realslice");
        hsize_t mdims[3];
        mdims[0] = {testArr3D.get_dimi()};
        mdims[1] = {testArr3D.get_dimj()};
        mdims[2] = {testArr3D.get_dimk()};
        writeDatacube3D(test3D_data, &testArr3D[0],mdims);
        test3D_data.close();
    }

    {
        H5::DataSet test2D_data = realslices.openDataSet("test2Darray/realslice");
        hsize_t mdims[2];
        mdims[0] = {testArr2D.get_dimi()};
        mdims[1] = {testArr2D.get_dimj()};
        writeRealSlice(test2D_data, &testArr2D[0],mdims);
        test2D_data.close();
    }
    realslices.close();
    pars.outputFile.close();

    std::string dataPath2D = "4DSTEM_simulation/data/realslices/test2Darray/realslice";
    std::string dataPath3D = "4DSTEM_simulation/data/realslices/test3Darray/realslice";
    std::string dataPath4D = "4DSTEM_simulation/data/datacubes/test4Darray/datacube";

    Array2D<PRISMATIC_FLOAT_PRECISION> read2D = readDataSet2D(pars.meta.filenameOutput, dataPath2D);
    Array3D<PRISMATIC_FLOAT_PRECISION> read3D = readDataSet3D(pars.meta.filenameOutput, dataPath3D);
    Array4D<PRISMATIC_FLOAT_PRECISION> read4D = readDataSet4D(pars.meta.filenameOutput, dataPath4D);

    //check for array sizing equivalance
    BOOST_TEST(compareSize(read2D, testArr2D));
    BOOST_TEST(compareSize(read3D, testArr3D));
    BOOST_TEST(compareSize(read4D, testArr4D));

    //check value equivalence
    PRISMATIC_FLOAT_PRECISION tol = 0.00001;
    BOOST_TEST(compareValues(read2D, testArr2D) < tol);
    BOOST_TEST(compareValues(read3D, testArr3D) < tol);
    BOOST_TEST(compareValues(read4D, testArr4D) < tol);

    removeFile(pars.meta.filenameOutput);

};

BOOST_FIXTURE_TEST_CASE(importPotential2D_P, basicSim)
{
    //run simulations

    meta.potential3D = false;

    divertOutput(pos, fd, logPath);
    std::cout << "\n#### BEGIN TEST CASE: importPotential2D_P #####\n";

    std::string importFile = "../test/potentialImport.h5";
    meta.filenameOutput = "../test/potentialImport.h5";
    go(meta);

    std::cout << "\n--------------------------------------------\n";

    meta.filenameOutput = "../test/potentialRerun.h5";
    meta.importFile     = "../test/potentialImport.h5";
    meta.importPath     = "4DSTEM_simulation/data/realslices/ppotential_fp0000/realslice";
    meta.importPotential = true;
    go(meta);
    std::cout << "###### END TEST CASE: importPotential2D_P #####\n";

    revertOutput(fd, pos);

    //read in output arrays and compare
    std::string dataPath2D = "4DSTEM_simulation/data/realslices/annular_detector_depth0000/realslice";
    std::string dataPathDPC = "4DSTEM_simulation/data/realslices/DPC_CoM_depth0000/realslice";
    std::string dataPath3D = "4DSTEM_simulation/data/realslices/virtual_detector_depth0000/realslice";
    std::string dataPath4D = "4DSTEM_simulation/data/datacubes/CBED_array_depth0000/datacube";
    std::string dataPathPS = "4DSTEM_simulation/data/realslices/ppotential_fp0000/realslice";

    Array2D<PRISMATIC_FLOAT_PRECISION> refAnnular = readDataSet2D(importFile, dataPath2D);
    Array3D<PRISMATIC_FLOAT_PRECISION> refPS = readDataSet3D(importFile, dataPathPS);
    Array3D<PRISMATIC_FLOAT_PRECISION> refDPC = readDataSet3D(importFile, dataPathDPC);
    Array3D<PRISMATIC_FLOAT_PRECISION> refVD = readDataSet3D(importFile, dataPath3D);
    Array4D<PRISMATIC_FLOAT_PRECISION> refCBED = readDataSet4D(importFile, dataPath4D);

    Array2D<PRISMATIC_FLOAT_PRECISION> testAnnular = readDataSet2D(meta.filenameOutput, dataPath2D);
    Array3D<PRISMATIC_FLOAT_PRECISION> testPS = readDataSet3D(meta.filenameOutput, dataPathPS);
    Array3D<PRISMATIC_FLOAT_PRECISION> testDPC = readDataSet3D(meta.filenameOutput, dataPathDPC);
    Array3D<PRISMATIC_FLOAT_PRECISION> testVD = readDataSet3D(meta.filenameOutput, dataPath3D);
    Array4D<PRISMATIC_FLOAT_PRECISION> testCBED = readDataSet4D(meta.filenameOutput, dataPath4D);

    PRISMATIC_FLOAT_PRECISION tol = 0.001;
    PRISMATIC_FLOAT_PRECISION errorSum = 0.0;

    BOOST_TEST(compareSize(refAnnular, testAnnular));
    BOOST_TEST(compareSize(refPS, testPS));
    BOOST_TEST(compareSize(refVD, testVD));
    BOOST_TEST(compareSize(refDPC, testDPC));
    BOOST_TEST(compareSize(refCBED, testCBED));

    BOOST_TEST(compareValues(refAnnular, testAnnular) < tol);
    BOOST_TEST(compareValues(refPS, testPS) < tol);
    BOOST_TEST(compareValues(refDPC, testDPC) < tol);
    BOOST_TEST(compareValues(refVD, testVD) < tol);
    BOOST_TEST(compareValues(refCBED, testCBED) < tol);

    removeFile(importFile);
    removeFile(meta.filenameOutput);
};

BOOST_FIXTURE_TEST_CASE(importPotential3D_P, basicSim)
{
    //run simulations

    meta.potential3D = true;

    divertOutput(pos, fd, logPath);
    std::cout << "\n#### BEGIN TEST CASE: importPotential3D_P #####\n";

    std::string importFile = "../test/potentialImport.h5";
    meta.filenameOutput = "../test/potentialImport.h5";
    go(meta);

    std::cout << "\n--------------------------------------------\n";

    meta.filenameOutput = "../test/potentialRerun.h5";
    meta.importFile     = "../test/potentialImport.h5";
    meta.importPath     = "4DSTEM_simulation/data/realslices/ppotential_fp0000/realslice";
    meta.importPotential = true;
    go(meta);
    std::cout << "###### END TEST CASE: importPotential3D_P #####\n";

    revertOutput(fd, pos);

    //read in output arrays and compare
    std::string dataPath2D = "4DSTEM_simulation/data/realslices/annular_detector_depth0000/realslice";
    std::string dataPathDPC = "4DSTEM_simulation/data/realslices/DPC_CoM_depth0000/realslice";
    std::string dataPath3D = "4DSTEM_simulation/data/realslices/virtual_detector_depth0000/realslice";
    std::string dataPath4D = "4DSTEM_simulation/data/datacubes/CBED_array_depth0000/datacube";
    std::string dataPathPS = "4DSTEM_simulation/data/realslices/ppotential_fp0000/realslice";

    Array2D<PRISMATIC_FLOAT_PRECISION> refAnnular = readDataSet2D(importFile, dataPath2D);
    Array3D<PRISMATIC_FLOAT_PRECISION> refPS = readDataSet3D(importFile, dataPathPS);
    Array3D<PRISMATIC_FLOAT_PRECISION> refDPC = readDataSet3D(importFile, dataPathDPC);
    Array3D<PRISMATIC_FLOAT_PRECISION> refVD = readDataSet3D(importFile, dataPath3D);
    Array4D<PRISMATIC_FLOAT_PRECISION> refCBED = readDataSet4D(importFile, dataPath4D);

    Array2D<PRISMATIC_FLOAT_PRECISION> testAnnular = readDataSet2D(meta.filenameOutput, dataPath2D);
    Array3D<PRISMATIC_FLOAT_PRECISION> testPS = readDataSet3D(meta.filenameOutput, dataPathPS);
    Array3D<PRISMATIC_FLOAT_PRECISION> testDPC = readDataSet3D(meta.filenameOutput, dataPathDPC);
    Array3D<PRISMATIC_FLOAT_PRECISION> testVD = readDataSet3D(meta.filenameOutput, dataPath3D);
    Array4D<PRISMATIC_FLOAT_PRECISION> testCBED = readDataSet4D(meta.filenameOutput, dataPath4D);

    PRISMATIC_FLOAT_PRECISION tol = 0.001;
    PRISMATIC_FLOAT_PRECISION errorSum = 0.0;

    BOOST_TEST(compareSize(refAnnular, testAnnular));
    BOOST_TEST(compareSize(refPS, testPS));
    BOOST_TEST(compareSize(refVD, testVD));
    BOOST_TEST(compareSize(refDPC, testDPC));
    BOOST_TEST(compareSize(refCBED, testCBED));

    BOOST_TEST(compareValues(refAnnular, testAnnular) < tol);
    BOOST_TEST(compareValues(refPS, testPS) < tol);
    BOOST_TEST(compareValues(refDPC, testDPC) < tol);
    BOOST_TEST(compareValues(refVD, testVD) < tol);
    BOOST_TEST(compareValues(refCBED, testCBED) < tol);

    removeFile(importFile);
    removeFile(meta.filenameOutput);
};

BOOST_FIXTURE_TEST_CASE(importPotential2D_M, basicSim)
{
    //run simulations
    meta.potential3D = false;
    meta.algorithm = Algorithm::Multislice;
    meta.numGPUs =0;

    divertOutput(pos, fd, logPath);
    std::cout << "\n#### BEGIN TEST CASE: importPotential2D_M #####\n";

    std::string importFile = "../test/potentialImport.h5";
    meta.filenameOutput = "../test/potentialImport.h5";
    go(meta);

    std::cout << "\n--------------------------------------------\n";

    meta.filenameOutput = "../test/potentialRerun.h5";
    meta.importFile     = "../test/potentialImport.h5";
    meta.importPath     = "4DSTEM_simulation/data/realslices/ppotential_fp0000/realslice";
    meta.importPotential = true;
    go(meta);
    std::cout << "###### END TEST CASE: importPotential2D_M #####\n";

    revertOutput(fd, pos);

    //read in output arrays and compare
    std::string dataPath2D = "4DSTEM_simulation/data/realslices/annular_detector_depth0000/realslice";
    std::string dataPathDPC = "4DSTEM_simulation/data/realslices/DPC_CoM_depth0000/realslice";
    std::string dataPath3D = "4DSTEM_simulation/data/realslices/virtual_detector_depth0000/realslice";
    std::string dataPath4D = "4DSTEM_simulation/data/datacubes/CBED_array_depth0000/datacube";
    std::string dataPathPS = "4DSTEM_simulation/data/realslices/ppotential_fp0000/realslice";

    Array2D<PRISMATIC_FLOAT_PRECISION> refAnnular = readDataSet2D(importFile, dataPath2D);
    Array3D<PRISMATIC_FLOAT_PRECISION> refPS = readDataSet3D(importFile, dataPathPS);
    Array3D<PRISMATIC_FLOAT_PRECISION> refDPC = readDataSet3D(importFile, dataPathDPC);
    Array3D<PRISMATIC_FLOAT_PRECISION> refVD = readDataSet3D(importFile, dataPath3D);
    Array4D<PRISMATIC_FLOAT_PRECISION> refCBED = readDataSet4D(importFile, dataPath4D);

    Array2D<PRISMATIC_FLOAT_PRECISION> testAnnular = readDataSet2D(meta.filenameOutput, dataPath2D);
    Array3D<PRISMATIC_FLOAT_PRECISION> testPS = readDataSet3D(meta.filenameOutput, dataPathPS);
    Array3D<PRISMATIC_FLOAT_PRECISION> testDPC = readDataSet3D(meta.filenameOutput, dataPathDPC);
    Array3D<PRISMATIC_FLOAT_PRECISION> testVD = readDataSet3D(meta.filenameOutput, dataPath3D);
    Array4D<PRISMATIC_FLOAT_PRECISION> testCBED = readDataSet4D(meta.filenameOutput, dataPath4D);

    PRISMATIC_FLOAT_PRECISION tol = 0.001;
    PRISMATIC_FLOAT_PRECISION errorSum = 0.0;

    BOOST_TEST(compareSize(refAnnular, testAnnular));
    BOOST_TEST(compareSize(refPS, testPS));
    BOOST_TEST(compareSize(refVD, testVD));
    BOOST_TEST(compareSize(refDPC, testDPC));
    BOOST_TEST(compareSize(refCBED, testCBED));

    BOOST_TEST(compareValues(refAnnular, testAnnular) < tol);
    BOOST_TEST(compareValues(refPS, testPS) < tol);
    BOOST_TEST(compareValues(refDPC, testDPC) < tol);
    BOOST_TEST(compareValues(refVD, testVD) < tol);
    BOOST_TEST(compareValues(refCBED, testCBED) < tol);

    removeFile(importFile);
    removeFile(meta.filenameOutput);
};

BOOST_FIXTURE_TEST_CASE(importPotential3D_M, basicSim)
{
    //run simulations

    meta.potential3D = true;
    meta.algorithm = Algorithm::Multislice;

    divertOutput(pos, fd, logPath);
    std::cout << "\n#### BEGIN TEST CASE: importPotential3D_M #####\n";

    std::string importFile = "../test/potentialImport.h5";
    meta.filenameOutput = "../test/potentialImport.h5";
    go(meta);

    std::cout << "\n--------------------------------------------\n";

    meta.filenameOutput = "../test/potentialRerun.h5";
    meta.importFile     = "../test/potentialImport.h5";
    meta.importPath     = "4DSTEM_simulation/data/realslices/ppotential_fp0000/realslice";
    meta.importPotential = true;
    go(meta);
    std::cout << "###### END TEST CASE: importPotential3D_M #####\n";

    revertOutput(fd, pos);

    //read in output arrays and compare
    std::string dataPath2D = "4DSTEM_simulation/data/realslices/annular_detector_depth0000/realslice";
    std::string dataPathDPC = "4DSTEM_simulation/data/realslices/DPC_CoM_depth0000/realslice";
    std::string dataPath3D = "4DSTEM_simulation/data/realslices/virtual_detector_depth0000/realslice";
    std::string dataPath4D = "4DSTEM_simulation/data/datacubes/CBED_array_depth0000/datacube";
    std::string dataPathPS = "4DSTEM_simulation/data/realslices/ppotential_fp0000/realslice";

    Array2D<PRISMATIC_FLOAT_PRECISION> refAnnular = readDataSet2D(importFile, dataPath2D);
    Array3D<PRISMATIC_FLOAT_PRECISION> refPS = readDataSet3D(importFile, dataPathPS);
    Array3D<PRISMATIC_FLOAT_PRECISION> refDPC = readDataSet3D(importFile, dataPathDPC);
    Array3D<PRISMATIC_FLOAT_PRECISION> refVD = readDataSet3D(importFile, dataPath3D);
    Array4D<PRISMATIC_FLOAT_PRECISION> refCBED = readDataSet4D(importFile, dataPath4D);

    Array2D<PRISMATIC_FLOAT_PRECISION> testAnnular = readDataSet2D(meta.filenameOutput, dataPath2D);
    Array3D<PRISMATIC_FLOAT_PRECISION> testPS = readDataSet3D(meta.filenameOutput, dataPathPS);
    Array3D<PRISMATIC_FLOAT_PRECISION> testDPC = readDataSet3D(meta.filenameOutput, dataPathDPC);
    Array3D<PRISMATIC_FLOAT_PRECISION> testVD = readDataSet3D(meta.filenameOutput, dataPath3D);
    Array4D<PRISMATIC_FLOAT_PRECISION> testCBED = readDataSet4D(meta.filenameOutput, dataPath4D);

    PRISMATIC_FLOAT_PRECISION tol = 0.001;
    PRISMATIC_FLOAT_PRECISION errorSum = 0.0;

    BOOST_TEST(compareSize(refAnnular, testAnnular));
    BOOST_TEST(compareSize(refPS, testPS));
    BOOST_TEST(compareSize(refVD, testVD));
    BOOST_TEST(compareSize(refDPC, testDPC));
    BOOST_TEST(compareSize(refCBED, testCBED));

    BOOST_TEST(compareValues(refAnnular, testAnnular) < tol);
    BOOST_TEST(compareValues(refPS, testPS) < tol);
    BOOST_TEST(compareValues(refDPC, testDPC) < tol);
    BOOST_TEST(compareValues(refVD, testVD) < tol);
    BOOST_TEST(compareValues(refCBED, testCBED) < tol);

    removeFile(importFile);
    removeFile(meta.filenameOutput);
};

BOOST_FIXTURE_TEST_CASE(fourierResampling, basicSim)
{   
    meta.potential3D = false;
    //make larger to have better test for resampling
    meta.realspacePixelSize[0] = 0.06; 
    meta.realspacePixelSize[1] = 0.06;
    meta.numGPUs = 0;
    meta.algorithm = Algorithm::Multislice;
    meta.E0 = 200e3;
    
    divertOutput(pos, fd, logPath);
    std::cout << "\n##### BEGIN TEST CASE: fourierResampling ######\n";

    std::string importFile = "../test/potentialImport.h5";
    meta.filenameOutput = "../test/potentialImport.h5";
    go(meta);

    std::cout << "\n--------------------------------------------\n";

    meta.algorithm = Algorithm::PRISM;
    meta.interpolationFactorX = 5;
    meta.interpolationFactorY = 7;
    meta.filenameOutput = "../test/potentialRerun.h5";
    meta.importFile     = "../test/potentialImport.h5";
    meta.importPath     = "4DSTEM_simulation/data/realslices/ppotential_fp0000/realslice";
    meta.importPotential = true;
    go(meta);
    std::cout << "####### END TEST CASE: fourierResampling ######\n";

    revertOutput(fd, pos);

    //read in output arrays and compare
    std::string dataPathPS = "4DSTEM_simulation/data/realslices/ppotential_fp0000/realslice";

    Array3D<PRISMATIC_FLOAT_PRECISION> refPS = readDataSet3D(importFile, dataPathPS);
    Array3D<PRISMATIC_FLOAT_PRECISION> testPS = readDataSet3D(meta.filenameOutput, dataPathPS);

    //double because mean is not precise enough within single float
    double tol = 0.00001;
    double ref_mean = 0.0;
    double test_mean = 0.0;

    //only valid numerical test is to compare means of potential arrays, since change in qx, qy alter binning in other outputs
    for(auto i = 0; i < refPS.size(); i++) ref_mean += refPS[i];
    for(auto i = 0; i < testPS.size(); i++) test_mean += testPS[i];

    std::cout << ref_mean << " " << test_mean << std::endl;    
    ref_mean /= refPS.size();
    test_mean /= testPS.size();
    std::cout << ref_mean << " " << test_mean << std::endl;
    BOOST_TEST(std::abs(ref_mean-test_mean) < tol);

    //resampled potential slices should have 80 x 84 x 3 dims do align with fx = 5, fy = 7
    //rewrite to auto check multiples of 4x(fx, fy)
    bool sizeCheck = testPS.get_dimi() == 80 && testPS.get_dimj() == 84 && testPS.get_dimk() == 3;
    BOOST_TEST(sizeCheck);

    removeFile(importFile);
    removeFile(meta.filenameOutput);

};

BOOST_FIXTURE_TEST_CASE(importSMatrix, basicSim)
{
    //run simulations

    meta.potential3D = false;
    meta.saveSMatrix = true;
    meta.savePotentialSlices = false;
    meta.interpolationFactorX = 1;
    meta.interpolationFactorY = 1;

    divertOutput(pos, fd, logPath);
    std::cout << "\n####### BEGIN TEST CASE: importSMatrix ########\n";

    std::string importFile = "../test/smatrixImport.h5";
    meta.filenameOutput = "../test/smatrixImport.h5";
    go(meta);

    std::cout << "\n--------------------------------------------\n";

    meta.filenameOutput = "../test/smatrixRerun.h5";
    meta.importFile     = "../test/smatrixImport.h5";
    meta.importSMatrix = true;
    go(meta);
    std::cout << "######### END TEST CASE: importSMatrix ########\n";

    revertOutput(fd, pos);

    //read in output arrays and compare
    std::string dataPath2D = "4DSTEM_simulation/data/realslices/annular_detector_depth0000/realslice";
    std::string dataPathDPC = "4DSTEM_simulation/data/realslices/DPC_CoM_depth0000/realslice";
    std::string dataPath3D = "4DSTEM_simulation/data/realslices/virtual_detector_depth0000/realslice";
    std::string dataPath4D = "4DSTEM_simulation/data/datacubes/CBED_array_depth0000/datacube";
    std::string dataPathSM = "4DSTEM_simulation/data/realslices/smatrix_fp0000/realslice";
    std::vector<size_t> order_sm = {2,1,0};

    Array2D<PRISMATIC_FLOAT_PRECISION> refAnnular = readDataSet2D(importFile, dataPath2D);
    Array3D<PRISMATIC_FLOAT_PRECISION> refDPC = readDataSet3D(importFile, dataPathDPC);
    Array3D<PRISMATIC_FLOAT_PRECISION> refVD = readDataSet3D(importFile, dataPath3D);
    Array4D<PRISMATIC_FLOAT_PRECISION> refCBED = readDataSet4D(importFile, dataPath4D);
    Array3D<std::complex<PRISMATIC_FLOAT_PRECISION>> refSMatrix;
    readComplexDataSet(refSMatrix, importFile, dataPathSM, order_sm);

    Array2D<PRISMATIC_FLOAT_PRECISION> testAnnular = readDataSet2D(meta.filenameOutput, dataPath2D);
    Array3D<PRISMATIC_FLOAT_PRECISION> testDPC = readDataSet3D(meta.filenameOutput, dataPathDPC);
    Array3D<PRISMATIC_FLOAT_PRECISION> testVD = readDataSet3D(meta.filenameOutput, dataPath3D);
    Array4D<PRISMATIC_FLOAT_PRECISION> testCBED = readDataSet4D(meta.filenameOutput, dataPath4D);
    Array3D<std::complex<PRISMATIC_FLOAT_PRECISION>> testSMatrix;
    readComplexDataSet(testSMatrix, importFile, dataPathSM, order_sm);

    PRISMATIC_FLOAT_PRECISION tol = 0.001;

    BOOST_TEST(compareSize(refAnnular, testAnnular));
    BOOST_TEST(compareSize(refVD, testVD));
    BOOST_TEST(compareSize(refDPC, testDPC));
    BOOST_TEST(compareSize(refCBED, testCBED));
    BOOST_TEST(compareSize(refSMatrix, testSMatrix));

    BOOST_TEST(compareValues(refAnnular, testAnnular) < tol);
    BOOST_TEST(compareValues(refDPC, testDPC) < tol);
    BOOST_TEST(compareValues(refVD, testVD) < tol);
    BOOST_TEST(compareValues(refCBED, testCBED) < tol);
    BOOST_TEST(compareValues(refSMatrix, testSMatrix) < tol);

    removeFile(importFile);
    removeFile(meta.filenameOutput);
};

BOOST_FIXTURE_TEST_CASE(importSM_multFP, basicSim)
{
    //run simulations

    meta.potential3D = false;
    meta.saveSMatrix = true;
    meta.savePotentialSlices = false;
    meta.numFP = 4;
    meta.numGPUs = 1;

    divertOutput(pos, fd, logPath);
    std::cout << "\n###### BEGIN TEST CASE: importSM_multFP #######\n";

    std::string importFile = "../test/smatrixImport.h5";
    meta.filenameOutput = "../test/smatrixImport.h5";
    go(meta);

    std::cout << "\n--------------------------------------------\n";

    meta.filenameOutput = "../test/smatrixRerun.h5";
    meta.importFile     = "../test/smatrixImport.h5";
    meta.importSMatrix = true;
    go(meta);
    std::cout << "######## END TEST CASE: importSM_multFP #######\n";

    revertOutput(fd, pos);

    //read in output arrays and compare
    std::string dataPath2D = "4DSTEM_simulation/data/realslices/annular_detector_depth0000/realslice";
    std::string dataPathDPC = "4DSTEM_simulation/data/realslices/DPC_CoM_depth0000/realslice";
    std::string dataPath3D = "4DSTEM_simulation/data/realslices/virtual_detector_depth0000/realslice";
    std::string dataPath4D = "4DSTEM_simulation/data/datacubes/CBED_array_depth0000/datacube";

    Array2D<PRISMATIC_FLOAT_PRECISION> refAnnular = readDataSet2D(importFile, dataPath2D);
    Array3D<PRISMATIC_FLOAT_PRECISION> refDPC = readDataSet3D(importFile, dataPathDPC);
    Array3D<PRISMATIC_FLOAT_PRECISION> refVD = readDataSet3D(importFile, dataPath3D);
    Array4D<PRISMATIC_FLOAT_PRECISION> refCBED = readDataSet4D(importFile, dataPath4D);
    Array3D<std::complex<PRISMATIC_FLOAT_PRECISION>> refSMatrix;

    Array2D<PRISMATIC_FLOAT_PRECISION> testAnnular = readDataSet2D(meta.filenameOutput, dataPath2D);
    Array3D<PRISMATIC_FLOAT_PRECISION> testDPC = readDataSet3D(meta.filenameOutput, dataPathDPC);
    Array3D<PRISMATIC_FLOAT_PRECISION> testVD = readDataSet3D(meta.filenameOutput, dataPath3D);
    Array4D<PRISMATIC_FLOAT_PRECISION> testCBED = readDataSet4D(meta.filenameOutput, dataPath4D);
    Array3D<std::complex<PRISMATIC_FLOAT_PRECISION>> testSMatrix;

    PRISMATIC_FLOAT_PRECISION tol = 0.001;
    PRISMATIC_FLOAT_PRECISION errorSum = 0.0;

    BOOST_TEST(compareSize(refAnnular, testAnnular));
    BOOST_TEST(compareSize(refVD, testVD));
    BOOST_TEST(compareSize(refDPC, testDPC));
    BOOST_TEST(compareSize(refCBED, testCBED));

    BOOST_TEST(compareValues(refAnnular, testAnnular) < tol);
    BOOST_TEST(compareValues(refDPC, testDPC) < tol);
    BOOST_TEST(compareValues(refVD, testVD) < tol);
    BOOST_TEST(compareValues(refCBED, testCBED) < tol);

    for(auto i = 0; i < 4; i++)
    {
        std::string dataPathSM = "4DSTEM_simulation/data/realslices/smatrix_fp" + getDigitString(i) + "/realslice";
        std::cout << "Checking frozen phonon configuration: " << i << std::endl;
        std::vector<size_t> order_sm = {2,1,0};
        readComplexDataSet(refSMatrix, importFile, dataPathSM, order_sm);
        readComplexDataSet(testSMatrix, importFile, dataPathSM, order_sm);
        BOOST_TEST(compareSize(refSMatrix, testSMatrix));
        BOOST_TEST(compareValues(refSMatrix, testSMatrix) < tol);
    }

    removeFile(importFile);
    removeFile(meta.filenameOutput);
};

BOOST_FIXTURE_TEST_CASE(attributeTest, basicSim)
{
    meta.probeStepX = 1.0;
    meta.probeStepY = 1.0;
    meta.realspacePixelSize[0] = 0.5;
    meta.realspacePixelSize[1] = 0.5;
    PRISMATIC_FLOAT_PRECISION tol = 0.00001;

    divertOutput(pos, fd, logPath);
    std::cout << "\n####### BEGIN TEST CASE: attributeTest ########\n";
    go(meta);
    std::cout << "######### END TEST CASE: attributeTest ########\n";
    revertOutput(fd, pos);

    std::string groupPath = "4DSTEM_simulation/metadata/metadata_0/original/simulation_parameters";
    std::string attProbe = "rx";
    std::string attCellDim = "c";
    std::string attFx = "fx";
    std::string attInput = "i";

    PRISMATIC_FLOAT_PRECISION probeStepCheck;
    PRISMATIC_FLOAT_PRECISION cellDimCheck[3];
    int fxCheck;
    std::string inputCheck;

    readAttribute(meta.filenameOutput, groupPath, attProbe, probeStepCheck);
    readAttribute(meta.filenameOutput, groupPath, attCellDim, cellDimCheck);
    readAttribute(meta.filenameOutput, groupPath, attFx, fxCheck);
    readAttribute(meta.filenameOutput, groupPath, attInput, inputCheck);
    
    BOOST_TEST(probeStepCheck == meta.probeStepX);
    BOOST_TEST(fxCheck == meta.interpolationFactorX);
    BOOST_TEST(inputCheck == meta.filenameAtoms);

    //meta is not updated in place consistnetly through a full sim; compare against known value
    int errCheck = 0;
    for(auto i = 0; i < 3; i ++) errCheck += (std::abs(cellDimCheck[i]-5.43) < tol) ? 0 : 1;
    BOOST_TEST(errCheck < 1);

    removeFile(meta.filenameOutput);
};

BOOST_FIXTURE_TEST_CASE(importPot_multipleFP_P, basicSim)
{
    //run simulations

    meta.potential3D = false;
    meta.numFP = 4;
    meta.includeThermalEffects = 1;
    
    divertOutput(pos, fd, logPath);
    std::cout << "\n### BEGIN TEST CASE: importPot_multipleFP_P ###\n";

    std::string importFile = "../test/potentialImport.h5";
    meta.filenameOutput = "../test/potentialImport.h5";
    go(meta);

    std::cout << "\n--------------------------------------------\n";

    meta.filenameOutput = "../test/potentialRerun.h5";
    meta.importFile     = "../test/potentialImport.h5";
    meta.importPotential = true;
    go(meta);
    std::cout << "#### END TEST CASE: importPot_multipleFP_P ####\n";

    revertOutput(fd, pos);

    //read in output arrays and compare
    std::string dataPath2D = "4DSTEM_simulation/data/realslices/annular_detector_depth0000/realslice";
    std::string dataPathDPC = "4DSTEM_simulation/data/realslices/DPC_CoM_depth0000/realslice";
    std::string dataPath3D = "4DSTEM_simulation/data/realslices/virtual_detector_depth0000/realslice";
    std::string dataPath4D = "4DSTEM_simulation/data/datacubes/CBED_array_depth0000/datacube";

    Array2D<PRISMATIC_FLOAT_PRECISION> refAnnular = readDataSet2D(importFile, dataPath2D);
    Array3D<PRISMATIC_FLOAT_PRECISION> refDPC = readDataSet3D(importFile, dataPathDPC);
    Array3D<PRISMATIC_FLOAT_PRECISION> refVD = readDataSet3D(importFile, dataPath3D);
    Array4D<PRISMATIC_FLOAT_PRECISION> refCBED = readDataSet4D(importFile, dataPath4D);

    Array2D<PRISMATIC_FLOAT_PRECISION> testAnnular = readDataSet2D(meta.filenameOutput, dataPath2D);
    Array3D<PRISMATIC_FLOAT_PRECISION> testDPC = readDataSet3D(meta.filenameOutput, dataPathDPC);
    Array3D<PRISMATIC_FLOAT_PRECISION> testVD = readDataSet3D(meta.filenameOutput, dataPath3D);
    Array4D<PRISMATIC_FLOAT_PRECISION> testCBED = readDataSet4D(meta.filenameOutput, dataPath4D);

    PRISMATIC_FLOAT_PRECISION tol = 0.001;
    PRISMATIC_FLOAT_PRECISION errorSum = 0.0;

    BOOST_TEST(compareSize(refAnnular, testAnnular));
    BOOST_TEST(compareSize(refVD, testVD));
    BOOST_TEST(compareSize(refDPC, testDPC));
    BOOST_TEST(compareSize(refCBED, testCBED));

    BOOST_TEST(compareValues(refAnnular, testAnnular) < tol);
    BOOST_TEST(compareValues(refDPC, testDPC) < tol);
    BOOST_TEST(compareValues(refVD, testVD) < tol);
    BOOST_TEST(compareValues(refCBED, testCBED) < tol);


    for(auto i = 0; i < 4; i++)
    {
        std::string dataPathPS = "4DSTEM_simulation/data/realslices/ppotential_fp" + getDigitString(i) + "/realslice";
        std::cout << "Checking frozen phonon configuration: " << i << std::endl;
        Array3D<PRISMATIC_FLOAT_PRECISION> refPS = readDataSet3D(importFile, dataPathPS);
        Array3D<PRISMATIC_FLOAT_PRECISION> testPS = readDataSet3D(meta.filenameOutput, dataPathPS);
        BOOST_TEST(compareSize(refPS, testPS));
        BOOST_TEST(compareValues(refPS, testPS) < tol);
    }

    removeFile(importFile);
    removeFile(meta.filenameOutput);
}

BOOST_FIXTURE_TEST_CASE(importPot_multipleFP_M, basicSim)
{
    //run simulations

    meta.potential3D = false;
    meta.numFP = 4;
    meta.includeThermalEffects = 1;
    meta.algorithm = Algorithm::Multislice;
    
    divertOutput(pos, fd, logPath);
    std::cout << "\n#### BEGIN TEST CASE: importPot_multipleFP_M ####\n";

    std::string importFile = "../test/potentialImport.h5";
    meta.filenameOutput = "../test/potentialImport.h5";
    go(meta);

    std::cout << "\n--------------------------------------------\n";

    meta.filenameOutput = "../test/potentialRerun.h5";
    meta.importFile     = "../test/potentialImport.h5";
    meta.importPotential = true;
    go(meta);
    std::cout << "#### END TEST CASE: importPot_multipleFP_M ####\n";

    revertOutput(fd, pos);

    //read in output arrays and compare
    std::string dataPath2D = "4DSTEM_simulation/data/realslices/annular_detector_depth0000/realslice";
    std::string dataPathDPC = "4DSTEM_simulation/data/realslices/DPC_CoM_depth0000/realslice";
    std::string dataPath3D = "4DSTEM_simulation/data/realslices/virtual_detector_depth0000/realslice";
    std::string dataPath4D = "4DSTEM_simulation/data/datacubes/CBED_array_depth0000/datacube";

    Array2D<PRISMATIC_FLOAT_PRECISION> refAnnular = readDataSet2D(importFile, dataPath2D);
    Array3D<PRISMATIC_FLOAT_PRECISION> refDPC = readDataSet3D(importFile, dataPathDPC);
    Array3D<PRISMATIC_FLOAT_PRECISION> refVD = readDataSet3D(importFile, dataPath3D);
    Array4D<PRISMATIC_FLOAT_PRECISION> refCBED = readDataSet4D(importFile, dataPath4D);

    Array2D<PRISMATIC_FLOAT_PRECISION> testAnnular = readDataSet2D(meta.filenameOutput, dataPath2D);
    Array3D<PRISMATIC_FLOAT_PRECISION> testDPC = readDataSet3D(meta.filenameOutput, dataPathDPC);
    Array3D<PRISMATIC_FLOAT_PRECISION> testVD = readDataSet3D(meta.filenameOutput, dataPath3D);
    Array4D<PRISMATIC_FLOAT_PRECISION> testCBED = readDataSet4D(meta.filenameOutput, dataPath4D);

    PRISMATIC_FLOAT_PRECISION tol = 0.001;
    PRISMATIC_FLOAT_PRECISION errorSum = 0.0;

    BOOST_TEST(compareSize(refAnnular, testAnnular));
    BOOST_TEST(compareSize(refVD, testVD));
    BOOST_TEST(compareSize(refDPC, testDPC));
    BOOST_TEST(compareSize(refCBED, testCBED));

    BOOST_TEST(compareValues(refAnnular, testAnnular) < tol);
    BOOST_TEST(compareValues(refDPC, testDPC) < tol);
    BOOST_TEST(compareValues(refVD, testVD) < tol);
    BOOST_TEST(compareValues(refCBED, testCBED) < tol);


    for(auto i = 0; i < 4; i++)
    {
        std::string dataPathPS = "4DSTEM_simulation/data/realslices/ppotential_fp" + getDigitString(i) + "/realslice";
        std::cout << "Checking frozen phonon configuration: " << i << std::endl;
        Array3D<PRISMATIC_FLOAT_PRECISION> refPS = readDataSet3D(importFile, dataPathPS);
        Array3D<PRISMATIC_FLOAT_PRECISION> testPS = readDataSet3D(meta.filenameOutput, dataPathPS);
        BOOST_TEST(compareSize(refPS, testPS));
        BOOST_TEST(compareValues(refPS, testPS) < tol);
    }

    removeFile(importFile);
    removeFile(meta.filenameOutput);
}

BOOST_FIXTURE_TEST_CASE(importPot_fpMismatch, basicSim)
{
    //run simulations

    meta.potential3D = false;
    meta.numFP = 4;
    meta.includeThermalEffects = 1;
    // meta.algorithm = Algorithm::Multislice;
    
    divertOutput(pos, fd, logPath);
    std::cout << "\n#### BEGIN TEST CASE: importPot_fpMismatch ####\n";

    std::string importFile = "../test/potentialImport.h5";
    meta.filenameOutput = "../test/potentialImport.h5";
    go(meta);

    std::cout << "\n--------------------------------------------\n";

    meta.filenameOutput = "../test/potentialRerun.h5";
    meta.importFile     = "../test/potentialImport.h5";
    meta.numFP = 1;
    meta.importPotential = true;
    go(meta);
    std::cout << "##### END TEST CASE: importPot_fpMismatch #####\n";

    revertOutput(fd, pos);

    //read in output arrays and compare
    std::string dataPath2D = "4DSTEM_simulation/data/realslices/annular_detector_depth0000/realslice";
    std::string dataPathDPC = "4DSTEM_simulation/data/realslices/DPC_CoM_depth0000/realslice";
    std::string dataPath3D = "4DSTEM_simulation/data/realslices/virtual_detector_depth0000/realslice";
    std::string dataPath4D = "4DSTEM_simulation/data/datacubes/CBED_array_depth0000/datacube";

    Array2D<PRISMATIC_FLOAT_PRECISION> refAnnular = readDataSet2D(importFile, dataPath2D);
    Array3D<PRISMATIC_FLOAT_PRECISION> refDPC = readDataSet3D(importFile, dataPathDPC);
    Array3D<PRISMATIC_FLOAT_PRECISION> refVD = readDataSet3D(importFile, dataPath3D);
    Array4D<PRISMATIC_FLOAT_PRECISION> refCBED = readDataSet4D(importFile, dataPath4D);

    Array2D<PRISMATIC_FLOAT_PRECISION> testAnnular = readDataSet2D(meta.filenameOutput, dataPath2D);
    Array3D<PRISMATIC_FLOAT_PRECISION> testDPC = readDataSet3D(meta.filenameOutput, dataPathDPC);
    Array3D<PRISMATIC_FLOAT_PRECISION> testVD = readDataSet3D(meta.filenameOutput, dataPath3D);
    Array4D<PRISMATIC_FLOAT_PRECISION> testCBED = readDataSet4D(meta.filenameOutput, dataPath4D);

    PRISMATIC_FLOAT_PRECISION tol = 0.001;
    PRISMATIC_FLOAT_PRECISION errorSum;
    PRISMATIC_FLOAT_PRECISION maxCBED = 0.0;

    BOOST_TEST(compareSize(refAnnular, testAnnular));
    BOOST_TEST(compareSize(refVD, testVD));
    BOOST_TEST(compareSize(refDPC, testDPC));
    BOOST_TEST(compareSize(refCBED, testCBED));

    BOOST_TEST(compareValues(refAnnular, testAnnular) < tol);
    BOOST_TEST(compareValues(refDPC, testDPC) < tol);
    BOOST_TEST(compareValues(refVD, testVD) < tol);
    BOOST_TEST(compareValues(refCBED, testCBED) < tol);

    int slices = (meta.numFP < 4) ? meta.numFP : 4;
    for(auto i = 0; i < slices; i++)
    {
        std::string dataPathPS = "4DSTEM_simulation/data/realslices/ppotential_fp" + getDigitString(i) + "/realslice";
        std::cout << "Checking frozen phonon configuration: " << i << std::endl;
        Array3D<PRISMATIC_FLOAT_PRECISION> refPS = readDataSet3D(importFile, dataPathPS);
        Array3D<PRISMATIC_FLOAT_PRECISION> testPS = readDataSet3D(meta.filenameOutput, dataPathPS);
        BOOST_TEST(compareSize(refPS, testPS));
        BOOST_TEST(compareValues(refPS, testPS) < tol);
    }

    removeFile(importFile);
    removeFile(meta.filenameOutput);
}

BOOST_AUTO_TEST_CASE(complexIO)
{
    //testing IO operations on complex datasets
    int seed = 10101;
    srand(seed);
    std::default_random_engine de(seed);

    Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> refArr2D = zeros_ND<2,std::complex<PRISMATIC_FLOAT_PRECISION>>({{2,7}});
    Array3D<std::complex<PRISMATIC_FLOAT_PRECISION>> refArr3D = zeros_ND<3,std::complex<PRISMATIC_FLOAT_PRECISION>>({{2,7,5}});
    Array4D<std::complex<PRISMATIC_FLOAT_PRECISION>> refArr4D = zeros_ND<4,std::complex<PRISMATIC_FLOAT_PRECISION>>({{2,7,5,3}});
    assignRandomValues(refArr2D, de);
    assignRandomValues(refArr3D, de);
    assignRandomValues(refArr4D, de);

    //create a test file
    std::string fname = "../test/testFile.h5";
    H5::H5File testFile = H5::H5File(fname.c_str(), H5F_ACC_TRUNC);
    H5::Group testGroup(testFile.createGroup("/complex_data"));

    hsize_t mdims_2D[2] = {refArr2D.get_dimi(), refArr2D.get_dimj()};
    hsize_t mdims_3D[3] = {refArr3D.get_dimi(), refArr3D.get_dimj(), refArr3D.get_dimk()};
    hsize_t mdims_4D[4] = {refArr4D.get_dimi(), refArr4D.get_dimj(), refArr4D.get_dimk(), refArr4D.get_diml()};

    std::vector<size_t> order_2D = {0,1}; 
    std::vector<size_t> order_3D = {0,1,2}; 
    std::vector<size_t> order_4D = {0,1,2,3}; 
    writeComplexDataSet(testGroup, "complex2D", &refArr2D[0], mdims_2D, 2, order_2D);
    writeComplexDataSet(testGroup, "complex3D", &refArr3D[0], mdims_3D, 3, order_3D);
    writeComplexDataSet(testGroup, "complex4D", &refArr4D[0], mdims_4D, 4, order_4D);

    testGroup.close();
    testFile.close();

    std::string datapath2D = "complex_data/complex2D";
    std::string datapath3D = "complex_data/complex3D";
    std::string datapath4D = "complex_data/complex4D";

    Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> testArr2D; 
    Array3D<std::complex<PRISMATIC_FLOAT_PRECISION>> testArr3D; 
    Array4D<std::complex<PRISMATIC_FLOAT_PRECISION>> testArr4D;
    readComplexDataSet(testArr2D, fname, datapath2D, order_2D);
    readComplexDataSet(testArr3D, fname, datapath3D, order_3D);
    readComplexDataSet(testArr4D, fname, datapath4D, order_4D);

    PRISMATIC_FLOAT_PRECISION tol = 0.0001;
    BOOST_TEST(compareSize(refArr2D, testArr2D));
    BOOST_TEST(compareSize(refArr3D, testArr3D));
    BOOST_TEST(compareSize(refArr4D, testArr4D));
    BOOST_TEST(compareValues(refArr2D, testArr2D) < tol);
    BOOST_TEST(compareValues(refArr3D, testArr3D) < tol);
    BOOST_TEST(compareValues(refArr4D, testArr4D) < tol);

    removeFile(fname);
}

BOOST_AUTO_TEST_CASE(dataGroupCount)
{
    //create a test file
    std::string fname = "../test/testFile.h5";
    H5::H5File testFile = H5::H5File(fname.c_str(), H5F_ACC_TRUNC);
    H5::Group testGroup(testFile.createGroup("/complex_data"));

    int numgroups = 10;
    std::string basename = "testName_";
    for(auto i = 0; i < numgroups; i++)
    {
        H5::Group newGroup(testGroup.createGroup(basename+getDigitString(i)));
        newGroup.close();
    }

    int counted = countDataGroups(testGroup, basename);

    testGroup.close();
    testFile.close();

    BOOST_TEST(numgroups == counted);

    removeFile(fname);
}

BOOST_AUTO_TEST_CASE(virtualDataSet)
{
    //create test file
    std::string fname = "../test/testFile.h5";
    Parameters<PRISMATIC_FLOAT_PRECISION> pars;
    pars.outputFile = H5::H5File(fname.c_str(), H5F_ACC_TRUNC);;
    setupOutputFile(pars);

    //create component datasets
    size_t ydim = 10;
    size_t xdim = 10;
    Array2D<PRISMATIC_FLOAT_PRECISION> one = ones_ND<2,PRISMATIC_FLOAT_PRECISION>({{ydim, xdim}});
    Array2D<PRISMATIC_FLOAT_PRECISION> two = one*2;
    Array2D<PRISMATIC_FLOAT_PRECISION> three = one*3;
    Array2D<PRISMATIC_FLOAT_PRECISION> four = one*4;

    //write to file
    H5::Group dataHold(pars.outputFile.createGroup("/4DSTEM_simulation/dataHold"));
    hsize_t msize[2] = {ydim, xdim};
    H5::DataSpace mspace(2, msize);
    H5::DataSpace fspace;

    H5::DataSet one_dset = dataHold.createDataSet("one", PFP_TYPE, mspace);
    H5::DataSet two_dset = dataHold.createDataSet("two", PFP_TYPE, mspace);
    H5::DataSet three_dset = dataHold.createDataSet("three", PFP_TYPE, mspace);
    H5::DataSet four_dset = dataHold.createDataSet("four", PFP_TYPE, mspace);

    fspace = one_dset.getSpace();
    one_dset.write(&one[0], PFP_TYPE, mspace, fspace);    
    fspace = two_dset.getSpace();
    two_dset.write(&two[0], PFP_TYPE, mspace, fspace);
    fspace = three_dset.getSpace();
    three_dset.write(&three[0], PFP_TYPE, mspace, fspace);
    fspace = four_dset.getSpace();
    four_dset.write(&four[0], PFP_TYPE, mspace, fspace);


    std::string sgName = "testSG";
    std::vector<H5::DataSet> datasets{one_dset, two_dset, three_dset, four_dset};

    //wrtie virtual dataset as 4x10x10
    std::vector<std::vector<size_t>> indices_3D;
    indices_3D.push_back(std::vector<size_t>{0});
    indices_3D.push_back(std::vector<size_t>{1});
    indices_3D.push_back(std::vector<size_t>{2});
    indices_3D.push_back(std::vector<size_t>{3});

    writeVirtualDataSet(dataHold, "testVDS", datasets, indices_3D);


    //write virtual dataset as 2x2x10x10
    std::vector<std::vector<size_t>> indices_4D;
    indices_4D.push_back(std::vector<size_t>{0,0});
    indices_4D.push_back(std::vector<size_t>{0,1});
    indices_4D.push_back(std::vector<size_t>{1,0});
    indices_4D.push_back(std::vector<size_t>{1,1});

    writeVirtualDataSet(dataHold, "testVD", datasets, indices_4D);
    pars.outputFile.close();

    PRISMATIC_FLOAT_PRECISION tol = 0.00001;

    Array3D<PRISMATIC_FLOAT_PRECISION> vds_read = readDataSet3D(fname, "/4DSTEM_simulation/dataHold/testVDS");
    { //restride
        Array3D<PRISMATIC_FLOAT_PRECISION> vds_tmp(vds_read);
        for(auto i = 0; i < vds_read.get_dimi(); i++)
		{
			for(auto j = 0; j < vds_read.get_dimj(); j++)
			{
				for(auto k = 0; k < vds_read.get_dimk(); k++)
				{
					vds_read[k*vds_read.get_dimi()*vds_read.get_dimj()+j*vds_read.get_dimi()+i] = vds_tmp[i*vds_read.get_dimk()*vds_read.get_dimj()+j*vds_read.get_dimk()+k];
				}
			}
		}
    }

    Array2D<PRISMATIC_FLOAT_PRECISION> one_read = subspace(vds_read, 0);
    Array2D<PRISMATIC_FLOAT_PRECISION> two_read = subspace(vds_read, 1);
    Array2D<PRISMATIC_FLOAT_PRECISION> three_read = subspace(vds_read, 2);
    Array2D<PRISMATIC_FLOAT_PRECISION> four_read = subspace(vds_read, 3);

    BOOST_TEST(compareSize(one, one_read));
    BOOST_TEST(compareValues(one, one_read) < tol);
    BOOST_TEST(compareSize(two, two_read));
    BOOST_TEST(compareValues(two, two_read) < tol);
    BOOST_TEST(compareSize(three, three_read));
    BOOST_TEST(compareValues(three, three_read) < tol);
    BOOST_TEST(compareSize(four, four_read));
    BOOST_TEST(compareValues(four, four_read) < tol);

    removeFile(fname);

}

BOOST_AUTO_TEST_CASE(datasetCopy)
{
    //make a random dataset, multidimensonal to test striding
    int seed = 23;
    srand(seed);
    std::default_random_engine de(seed);
    Array4D<PRISMATIC_FLOAT_PRECISION> refArr = zeros_ND<4, PRISMATIC_FLOAT_PRECISION>({{2,7,5,3}});
    Array4D<std::complex<PRISMATIC_FLOAT_PRECISION>> refArr_complex = zeros_ND<4, std::complex<PRISMATIC_FLOAT_PRECISION>>({{2,7,5,3}});
    assignRandomValues(refArr, de);
    assignRandomValues(refArr_complex, de);

    //create a test file
    std::string fname = "../test/testFile.h5";
    H5::H5File testFile = H5::H5File(fname.c_str(), H5F_ACC_TRUNC);
    H5::Group sourceGroup(testFile.createGroup("/source"));
    H5::Group targetGroup(testFile.createGroup("/target"));

    //create datasets with attributes in source group
    hsize_t mdims_4D[4] = {refArr.get_dimi(), refArr.get_dimj(), refArr.get_dimk(), refArr.get_diml()};
    H5::DataSpace mspace(4, mdims_4D);
    
    H5::DataSet refDS = sourceGroup.createDataSet("ref_ds", PFP_TYPE, mspace);
    H5::DataSpace fspace = refDS.getSpace();

    refDS.write(&refArr[0], PFP_TYPE, mspace, fspace);

    H5::DataSpace str_name_ds(H5S_SCALAR);
    H5::StrType strdatatype(H5::PredType::C_S1, 256);
    H5::Attribute ref_ds_name = refDS.createAttribute("name", strdatatype, str_name_ds);
    H5::Attribute ref_ds_unit = refDS.createAttribute("units", strdatatype, str_name_ds);

    const H5std_string ref_ds_name_str("i'm_the_name");
    const H5std_string ref_ds_unit_str("i'm_the_units");
    ref_ds_name.write(strdatatype, ref_ds_name_str);
    ref_ds_unit.write(strdatatype, ref_ds_unit_str);

    std::vector<size_t> order = {0,1,2,3};
    writeComplexDataSet(sourceGroup, "ref_ds_complex", &refArr_complex[0], mdims_4D, 4, order);
    H5::DataSet refDS_complex = sourceGroup.openDataSet("ref_ds_complex");

    //copy datasets then close all open objects
    copyDataSet(targetGroup, refDS);
    copyDataSet(targetGroup, refDS_complex);

    refDS.close();
    mspace.close();
    sourceGroup.close();
    targetGroup.close();
    testFile.close();

    //read dataset values and compare
    Array4D<PRISMATIC_FLOAT_PRECISION> refArr_read = readDataSet4D(fname, "/source/ref_ds");
    Array4D<PRISMATIC_FLOAT_PRECISION> testArr_read = readDataSet4D(fname, "/target/ref_ds");
    Array4D<std::complex<PRISMATIC_FLOAT_PRECISION>> refArr_complex_read;
    Array4D<std::complex<PRISMATIC_FLOAT_PRECISION>> testArr_complex_read;

    readComplexDataSet(refArr_complex_read, fname, "/source/ref_ds_complex", order);
    readComplexDataSet(testArr_complex_read, fname, "/target/ref_ds_complex", order);

    PRISMATIC_FLOAT_PRECISION tol = 0.00001;
    BOOST_TEST(compareSize(refArr_read, testArr_read));
    BOOST_TEST(compareSize(refArr_complex_read, testArr_complex_read));
    BOOST_TEST(compareValues(refArr_read, testArr_read) < tol);
    BOOST_TEST(compareValues(refArr_complex_read, testArr_complex_read) < tol);

    //reopen real datasets and compare attributes
    H5::H5File reOpen = H5::H5File(fname.c_str(), H5F_ACC_RDONLY);
    H5::DataSet src_read = reOpen.openDataSet("/source/ref_ds");
    H5::DataSet tar_read = reOpen.openDataSet("/target/ref_ds");

    for(auto i = 0; i < src_read.getNumAttrs(); i++)
    {
        H5::Attribute src_attr = src_read.openAttribute(i);
        BOOST_TEST(tar_read.attrExists(src_attr.getName()));
        if(tar_read.attrExists(src_attr.getName()))
        {
            //if exists, check data equivalence
            H5::Attribute tar_attr = tar_read.openAttribute(src_attr.getName());

            char src_buffer[src_attr.getInMemDataSize()];
            src_attr.read(src_attr.getDataType(), src_buffer);

            char tar_buffer[tar_attr.getInMemDataSize()];
            tar_attr.read(tar_attr.getDataType(), tar_buffer);
            std::string src_str = src_buffer;
            std::string tar_str = tar_buffer;
            bool val_check = src_str == tar_str;
            BOOST_TEST(val_check);
        }

    }

    removeFile(fname);
}

BOOST_FIXTURE_TEST_CASE(supergroup, basicSim)
{
    //generate virtual detector depth series
    std::string fname = "../test/supergroup.h5";
    meta.potential3D = false;
    meta.save2DOutput = false;
    meta.save4DOutput = false;
    meta.saveDPC_CoM = false;
    meta.numSlices = 1;
    meta.algorithm = Algorithm::Multislice;
    meta.filenameOutput = fname;
    meta.probeStepX = 0.4;
    meta.probeStepY = 0.3;
    meta.realspacePixelSize[0] = 0.25;
    meta.realspacePixelSize[1] = 0.1;

    divertOutput(pos, fd, logPath);
    std::cout << "\n######### BEGIN TEST CASE: supergroup #########\n";
    go(meta);
    std::cout << "########### END TEST CASE: supergroup #########\n";
    revertOutput(fd, pos);

    //create supergroup of depth series
    H5::H5File output = H5::H5File(fname.c_str(), H5F_ACC_RDWR);
    depthSeriesSG(output);

    //check equivalance of indvidual datasets to component datasets in VDS
    Array4D<PRISMATIC_FLOAT_PRECISION> vds_read = readDataSet4D_keepOrder(fname, "/4DSTEM_simulation/data/supergroups/vd_depth_series/supergroup");
    Array3D<PRISMATIC_FLOAT_PRECISION> vds_tmp = readDataSet3D(fname, "/4DSTEM_simulation/data/realslices/virtual_detector_depth0000/realslice");
    Array4D<PRISMATIC_FLOAT_PRECISION> ref_array = zeros_ND<4, PRISMATIC_FLOAT_PRECISION>({{3, vds_tmp.get_dimk(), vds_tmp.get_dimj(), vds_tmp.get_dimi()}});
    
    size_t strides = vds_tmp.get_dimk()*vds_tmp.get_dimj()*vds_tmp.get_dimi();
    for(auto i = 0; i < 3; i++)
    {   //read each vd and add to 4D array
        std::string path = "/4DSTEM_simulation/data/realslices/virtual_detector_depth" + getDigitString(i) + "/realslice";
        Array3D<PRISMATIC_FLOAT_PRECISION> vds_depth_i = readDataSet3D(fname, path);
        std::copy(vds_depth_i.begin(), vds_depth_i.end(), &ref_array[i*strides]);
    }

    PRISMATIC_FLOAT_PRECISION tol = 0.00001;
    std::array<size_t, 4> dims_in = {vds_read.get_dimarr()}; //as long as last dim is first, this should be good
    std::array<size_t, 4> order = {1, 2, 3, 0};
    vds_read = restride(vds_read, dims_in, order);

    BOOST_TEST(compareSize(vds_read, ref_array));
    BOOST_TEST(compareValues(vds_read, ref_array) < tol);
    
    //check source dims vs dims in supergroup
    //check total number of dims in group correspond to comparable ranks in datasets
    H5::Group testSG = output.openGroup("/4DSTEM_simulation/data/supergroups/vd_depth_series/");
    size_t dim_rank = 0;
    dim_rank += countDimensions(testSG, "dim");
    dim_rank += countDimensions(testSG, "sgdim");
    BOOST_TEST(dim_rank == vds_read.get_rank());

    std::vector<PRISMATIC_FLOAT_PRECISION> depth_check({2.0, 4.0, 6.0});
    std::vector<PRISMATIC_FLOAT_PRECISION> depth_read(3);
    H5::DataSet sgdim = testSG.openDataSet("sgdim1");

    hsize_t msize[1] = {3};
    H5::DataSpace depth_mspace(1, msize);
    H5::DataSpace depth_fspace = sgdim.getSpace();

	int rank = depth_fspace.getSimpleExtentNdims();
	hsize_t dims_out[rank];
	int ndims = depth_fspace.getSimpleExtentDims(dims_out, NULL); //nidms and rank are

    sgdim.read(&depth_read[0], PFP_TYPE, depth_mspace, depth_fspace);
    PRISMATIC_FLOAT_PRECISION errSum = std::abs(depth_check[0]-depth_read[0]);
    errSum += std::abs(depth_check[1]-depth_read[1]);
    errSum += std::abs(depth_check[2]-depth_read[2]);
    BOOST_TEST(errSum < tol);

    removeFile(fname);

}

BOOST_FIXTURE_TEST_CASE(complexOutputWave_M, basicSim)
{
    
    //run a complex output sim

    meta.potential3D = false;
    meta.algorithm = Algorithm::Multislice;
    meta.numGPUs = 1; //change later, test CPU first
    meta.filenameOutput = "../test/complexOutputWave_amplitude.h5";
    meta.savePotentialSlices = false;
    meta.alsoDoCPUWork = 0;
    meta.numGPUs = 1;
    meta.numStreamsPerGPU = 4;
    meta.probeStepX = 0.25;
    meta.probeStepY = 0.2;

    divertOutput(pos, fd, logPath);
    std::cout << "\n#### BEGIN TEST CASE: complexOutputWave_M #####\n";

    go(meta);

    std::cout << "\n--------------------------------------------\n";

    meta.filenameOutput = "../test/complexOutputWave.h5";
    meta.saveComplexOutputWave = true;
    go(meta);
    std::cout << "###### END TEST CASE: complexOutputWave_M #####\n";

    revertOutput(fd, pos);
    
    std::string complexFile = "../test/complexOutputWave.h5";
    std::string amplitudeFile = "../test/complexOutputWave_amplitude.h5";

    //read in output arrays
    std::string dataPath2D = "4DSTEM_simulation/data/realslices/annular_detector_depth0000/realslice";
    std::string dataPath3D = "4DSTEM_simulation/data/realslices/virtual_detector_depth0000/realslice";
    std::string dataPath4D = "4DSTEM_simulation/data/datacubes/CBED_array_depth0000/datacube";

    std::vector<size_t> order_2D = {0,1}; 
    std::vector<size_t> order_3D = {0,1,2}; 
    std::vector<size_t> order_4D = {2,3,0,1}; 

    Array2D<PRISMATIC_FLOAT_PRECISION> refAnnular;
    Array3D<PRISMATIC_FLOAT_PRECISION> refVD;
    Array4D<PRISMATIC_FLOAT_PRECISION> refCBED;
    readRealDataSet(refAnnular, amplitudeFile, dataPath2D, order_2D);
    readRealDataSet(refVD, amplitudeFile, dataPath3D, order_3D);
    readRealDataSet(refCBED, amplitudeFile, dataPath4D, order_4D);

    Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> testAnnular;
    Array3D<std::complex<PRISMATIC_FLOAT_PRECISION>> testVD;
    Array4D<std::complex<PRISMATIC_FLOAT_PRECISION>> testCBED;

    readComplexDataSet(testAnnular, complexFile, dataPath2D, order_2D);
    readComplexDataSet(testVD, complexFile, dataPath3D, order_3D);
    readComplexDataSet(testCBED, complexFile, dataPath4D, order_4D);

    //test that magnitude of values is equivalent
    Array2D<PRISMATIC_FLOAT_PRECISION> testAnnularAmp = getAmplitude(testAnnular);
    Array3D<PRISMATIC_FLOAT_PRECISION> testVDAmp = getAmplitude(testVD);
    Array4D<PRISMATIC_FLOAT_PRECISION> testCBEDAmp = getAmplitude(testCBED);

    BOOST_TEST(compareSize(refAnnular, testAnnularAmp));
    BOOST_TEST(compareSize(refVD, testVDAmp));
    BOOST_TEST(compareSize(refCBED, testCBEDAmp));

    PRISMATIC_FLOAT_PRECISION tol = 0.0001; //high because we are using total error
    BOOST_TEST(compareValues(refCBED, testCBEDAmp) < tol); //phase information is lost in 2D/3D, so can't compare directly w/o accounting for off diagonal errors

    removeFile(complexFile);
    removeFile(amplitudeFile);
}

BOOST_FIXTURE_TEST_CASE(complexOutputWave_P, basicSim)
{
    
    //run a complex output sim

    meta.potential3D = false;
    meta.algorithm = Algorithm::PRISM;
    meta.filenameOutput = "../test/complexOutputWave_amplitude.h5";
    meta.savePotentialSlices = false;
    meta.transferMode = StreamingMode::Stream;
    divertOutput(pos, fd, logPath);
    std::cout << "\n#### BEGIN TEST CASE: complexOutputWave_P #####\n";

    go(meta);

    std::cout << "\n--------------------------------------------\n";

    meta.filenameOutput = "../test/complexOutputWave.h5";
    meta.saveComplexOutputWave = true;
    go(meta);
    std::cout << "###### END TEST CASE: complexOutputWave_P #####\n";

    revertOutput(fd, pos);
    
    std::string complexFile = "../test/complexOutputWave.h5";
    std::string amplitudeFile = "../test/complexOutputWave_amplitude.h5";

    //read in output arrays
    std::string dataPath2D = "4DSTEM_simulation/data/realslices/annular_detector_depth0000/realslice";
    std::string dataPath3D = "4DSTEM_simulation/data/realslices/virtual_detector_depth0000/realslice";
    std::string dataPath4D = "4DSTEM_simulation/data/datacubes/CBED_array_depth0000/datacube";

    std::vector<size_t> order_2D = {0,1}; 
    std::vector<size_t> order_3D = {0,1,2}; 
    std::vector<size_t> order_4D = {2,3,0,1}; 

    Array2D<PRISMATIC_FLOAT_PRECISION> refAnnular;
    Array3D<PRISMATIC_FLOAT_PRECISION> refVD;
    Array4D<PRISMATIC_FLOAT_PRECISION> refCBED;
    readRealDataSet(refAnnular, amplitudeFile, dataPath2D, order_2D);
    readRealDataSet(refVD, amplitudeFile, dataPath3D, order_3D);
    readRealDataSet(refCBED, amplitudeFile, dataPath4D, order_4D);

    Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> testAnnular;
    Array3D<std::complex<PRISMATIC_FLOAT_PRECISION>> testVD;
    Array4D<std::complex<PRISMATIC_FLOAT_PRECISION>> testCBED;

    readComplexDataSet(testAnnular, complexFile, dataPath2D, order_2D);
    readComplexDataSet(testVD, complexFile, dataPath3D, order_3D);
    readComplexDataSet(testCBED, complexFile, dataPath4D, order_4D);

    //test that magnitude of values is equivalent
    Array2D<PRISMATIC_FLOAT_PRECISION> testAnnularAmp = getAmplitude(testAnnular);
    Array3D<PRISMATIC_FLOAT_PRECISION> testVDAmp = getAmplitude(testVD);
    Array4D<PRISMATIC_FLOAT_PRECISION> testCBEDAmp = getAmplitude(testCBED);

    BOOST_TEST(compareSize(refAnnular, testAnnularAmp));
    BOOST_TEST(compareSize(refVD, testVDAmp));
    BOOST_TEST(compareSize(refCBED, testCBEDAmp));

    PRISMATIC_FLOAT_PRECISION tol = 0.0001; //high because we are using total error; error seems to propagate with num streams 
    BOOST_TEST(compareValues(refCBED, testCBEDAmp) < tol); //phase information is lost in 2D/3D, so can't compare directly w/o accounting for off diagonal errors

    removeFile(complexFile);
    removeFile(amplitudeFile);
}

BOOST_AUTO_TEST_CASE(hdfStride)
{

    //set up a prime dimensioned dataset to test striding
    const size_t Nz = 2; const size_t Ny = 3; const size_t Nx = 5;
    Array3D<PRISMATIC_FLOAT_PRECISION> refData = zeros_ND<3, PRISMATIC_FLOAT_PRECISION>({{Nz, Ny, Nx}});
    for(auto i = 0; i < refData.size(); i++) refData[i] = i;

    //create a test file and write to test dataset
    std::string fname = "../test/hdfStride.h5";
    H5::H5File testFile = H5::H5File(fname.c_str(), H5F_ACC_TRUNC);
    hsize_t mdims[3] = {Nx, Ny, Nz};
    H5::DataSpace mspace(3,mdims);
    H5::DataSet testds = testFile.createDataSet("testds", PFP_TYPE, mspace);
    H5::DataSpace fspace = testds.getSpace();

    std::vector<size_t> vdims = {Nx,Ny,Nz};
    std::vector<size_t> vorder = {0,1,2};
    restrideElements(fspace, vdims, vorder);
    testds.write(&refData[0], PFP_TYPE, mspace, fspace);

    testds.close();
    testFile.close();

    //read from test dataset
    Array3D<PRISMATIC_FLOAT_PRECISION> testData;
    readRealDataSet(testData, fname, "testds", vorder);

    //compare values
    BOOST_TEST(compareValues(refData,testData) < 0.001);

    removeFile(fname);
}

BOOST_AUTO_TEST_CASE(CBEDoperator)
{
    //testing an operator to average the CBED arrays in place rather than stride, read-stride, add, write process
    
    //create function pointer for operator
    herr_t (*foo)(void*, hid_t, unsigned, const hsize_t*, void*);
    foo = &CBED_process;

    //set up test data
    int seed = 10101;
    srand(seed);
    std::default_random_engine de(seed);
    
    size_t Ny = 2; size_t Nx = 5;
    Array2D<PRISMATIC_FLOAT_PRECISION> testArr = ones_ND<2,PRISMATIC_FLOAT_PRECISION>({{Ny,Nx}});
    for(auto i = 0; i < Nx*Ny; i++) testArr[i] = i;

    for(auto i = 0; i < Ny; i++)
    {
        for(auto j =0; j < Nx; j++)
        {
            std::cout << testArr.at(i,j) << " ";
        }
        std::cout << std::endl;
    }

    Array2D<PRISMATIC_FLOAT_PRECISION> opData = zeros_ND<2,PRISMATIC_FLOAT_PRECISION>({{Ny,Nx}});
    assignRandomValues(opData, de);
    opData/=10;

    //create output file and store 'initial' configuration
    std::string fname = "../test/CBEDoperator.h5";
    H5::H5File testFile = H5::H5File(fname.c_str(), H5F_ACC_TRUNC);

    hsize_t mdims[2] = {Nx,Ny};
    H5::DataSpace mspace(2, mdims);
    H5::DataSet testds = testFile.createDataSet("testds", PFP_TYPE, mspace);
    
    H5::DataSpace fspace = testds.getSpace();
    std::vector<size_t> vdims = {Nx, Ny};
    std::vector<size_t> vorder = {0,1};
    restrideElements(fspace, vdims, vorder);
    testds.write(&testArr[0], PFP_TYPE, mspace, fspace);

    // actual procedure with averaging FP
    Array2D<PRISMATIC_FLOAT_PRECISION> readArr = zeros_ND<2,PRISMATIC_FLOAT_PRECISION>({{Ny,Nx}});
    testds.read(&readArr[0], PFP_TYPE, mspace, fspace);
    for(auto i =0; i < Nx*Ny; i++) readArr[i] += opData[i];
    testds.write(&readArr[0], PFP_TYPE, mspace, fspace);

    //read from disk for final output check
    testds.read(&readArr[0], PFP_TYPE, mspace, fspace);
    PRISMATIC_FLOAT_PRECISION errSum = 0;
    PRISMATIC_FLOAT_PRECISION tol = 0.00001;
    for(auto i =0; i < Nx*Ny; i++) errSum += std::abs(opData[i] - (readArr[i]-testArr[i]));
    BOOST_TEST((errSum/(Ny*Nx)) < tol);
    for(auto i = 0; i < Nx*Ny; i++)
    {
        std::cout << testArr[i] << " " << opData[i] << " " << readArr[i] << std::endl;
    }

    for(auto i = 0; i < Ny; i++)
    {
        for(auto j =0; j < Nx; j++)
        {
            std::cout << readArr.at(i,j) << " ";
        }
        std::cout << std::endl;
    }

}

BOOST_FIXTURE_TEST_CASE(fileSizeCheck, basicSim)
{
    meta.filenameOutput = "../test/fileSizeCheck.h5";
    // meta.filenameAtoms = "../test/au_np.xyz";
    meta.potential3D = false;
    meta.save2DOutput = false;
    meta.save3DOutput = true;
    meta.save4DOutput = false;
    meta.savePotentialSlices = false;
    meta.saveSMatrix = false;
    meta.saveDPC_CoM = false;
    meta.maxFileSize = 100;
    meta.algorithm = Algorithm::PRISM;

    divertOutput(pos, fd, logPath);
    bool errorCheck = false;
    std::cout << "\n####### BEGIN TEST CASE: fileSizeCheck ########\n";
    try
    {
        Parameters<PRISMATIC_FLOAT_PRECISION> testParams(meta);
    }
    catch (...)
    {
        errorCheck = true;
    }
    std::cout << "######### END TEST CASE: fileSizeCheck ########\n";
    revertOutput(fd, pos);

    BOOST_TEST(errorCheck);
}

BOOST_AUTO_TEST_SUITE_END();
}