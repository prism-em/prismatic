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
    PRISMATIC_FLOAT_PRECISION dummy = 1.0; //dummy float for IO overlaoding
    int seed = 10101;
    srand(seed);
    std::default_random_engine de(seed);

    Array2D<PRISMATIC_FLOAT_PRECISION> testArr2D = zeros_ND<2,PRISMATIC_FLOAT_PRECISION>({{2,7}});
    Array3D<PRISMATIC_FLOAT_PRECISION> testArr3D = zeros_ND<3,PRISMATIC_FLOAT_PRECISION>({{2,7,5}});
    Array4D<PRISMATIC_FLOAT_PRECISION> testArr4D = zeros_ND<4,PRISMATIC_FLOAT_PRECISION>({{2,7,5,3}});
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
        H5::DataSet test4D_data = test4D_group.createDataSet("datacube", H5::PredType::NATIVE_FLOAT, mspace);
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
                writeDatacube4D(pars,&testArr4D[arrayStart],mdims,offset,numFP,nameString.c_str());

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
        H5::DataSet test3D_data = test3D_group.createDataSet("realslice", H5::PredType::NATIVE_FLOAT, mspace);
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
        H5::DataSet test2D_data = test2D_group.createDataSet("realslice", H5::PredType::NATIVE_FLOAT, mspace);
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

    Array2D<PRISMATIC_FLOAT_PRECISION> read2D = readDataset2D(pars.meta.filenameOutput, dataPath2D);
    Array3D<PRISMATIC_FLOAT_PRECISION> read3D = readDataset3D(pars.meta.filenameOutput, dataPath3D);
    Array4D<PRISMATIC_FLOAT_PRECISION> read4D = readDataset4D(pars.meta.filenameOutput, dataPath4D);

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

BOOST_FIXTURE_TEST_CASE(importPotential2D, basicSim)
{
    //run simulations

    meta.potential3D = false;

    divertOutput(pos, fd, logPath);
    std::cout << "\n##### BEGIN TEST CASE: importPotential2D ######\n";

    std::string importFile = "../test/potentialImport.h5";
    meta.filenameOutput = "../test/potentialImport.h5";
    go(meta);

    std::cout << "\n--------------------------------------------\n";

    meta.filenameOutput = "../test/potentialRerun.h5";
    meta.importFile     = "../test/potentialImport.h5";
    meta.importPath     = "4DSTEM_simulation/data/realslices/ppotential/realslice";
    meta.importPotential = true;
    go(meta);
    std::cout << "####### END TEST CASE: importPotential2D ######\n";

    revertOutput(fd, pos);

    //read in output arrays and compare
    std::string dataPath2D = "4DSTEM_simulation/data/realslices/annular_detector_depth0000/realslice";
    std::string dataPathDPC = "4DSTEM_simulation/data/realslices/DPC_CoM_depth0000/realslice";
    std::string dataPath3D = "4DSTEM_simulation/data/realslices/virtual_detector_depth0000/realslice";
    std::string dataPath4D = "4DSTEM_simulation/data/datacubes/CBED_array_depth0000/datacube";
    std::string dataPathPS = "4DSTEM_simulation/data/realslices/ppotential/realslice";

    Array2D<PRISMATIC_FLOAT_PRECISION> refAnnular = readDataset2D(importFile, dataPath2D);
    Array3D<PRISMATIC_FLOAT_PRECISION> refPS = readDataset3D(importFile, dataPathPS);
    Array3D<PRISMATIC_FLOAT_PRECISION> refDPC = readDataset3D(importFile, dataPathDPC);
    Array3D<PRISMATIC_FLOAT_PRECISION> refVD = readDataset3D(importFile, dataPath3D);
    Array4D<PRISMATIC_FLOAT_PRECISION> refCBED = readDataset4D(importFile, dataPath4D);

    Array2D<PRISMATIC_FLOAT_PRECISION> testAnnular = readDataset2D(meta.filenameOutput, dataPath2D);
    Array3D<PRISMATIC_FLOAT_PRECISION> testPS = readDataset3D(meta.filenameOutput, dataPathPS);
    Array3D<PRISMATIC_FLOAT_PRECISION> testDPC = readDataset3D(meta.filenameOutput, dataPathDPC);
    Array3D<PRISMATIC_FLOAT_PRECISION> testVD = readDataset3D(meta.filenameOutput, dataPath3D);
    Array4D<PRISMATIC_FLOAT_PRECISION> testCBED = readDataset4D(meta.filenameOutput, dataPath4D);

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

BOOST_FIXTURE_TEST_CASE(importPotential3D, basicSim)
{
    //run simulations

    meta.potential3D = true;

    divertOutput(pos, fd, logPath);
    std::cout << "\n##### BEGIN TEST CASE: importPotential3D ######\n";

    std::string importFile = "../test/potentialImport.h5";
    meta.filenameOutput = "../test/potentialImport.h5";
    go(meta);

    std::cout << "\n--------------------------------------------\n";

    meta.filenameOutput = "../test/potentialRerun.h5";
    meta.importFile     = "../test/potentialImport.h5";
    meta.importPath     = "4DSTEM_simulation/data/realslices/ppotential/realslice";
    meta.importPotential = true;
    go(meta);
    std::cout << "####### END TEST CASE: importPotential3D ######\n";

    revertOutput(fd, pos);

    //read in output arrays and compare
    std::string dataPath2D = "4DSTEM_simulation/data/realslices/annular_detector_depth0000/realslice";
    std::string dataPathDPC = "4DSTEM_simulation/data/realslices/DPC_CoM_depth0000/realslice";
    std::string dataPath3D = "4DSTEM_simulation/data/realslices/virtual_detector_depth0000/realslice";
    std::string dataPath4D = "4DSTEM_simulation/data/datacubes/CBED_array_depth0000/datacube";
    std::string dataPathPS = "4DSTEM_simulation/data/realslices/ppotential/realslice";

    Array2D<PRISMATIC_FLOAT_PRECISION> refAnnular = readDataset2D(importFile, dataPath2D);
    Array3D<PRISMATIC_FLOAT_PRECISION> refPS = readDataset3D(importFile, dataPathPS);
    Array3D<PRISMATIC_FLOAT_PRECISION> refDPC = readDataset3D(importFile, dataPathDPC);
    Array3D<PRISMATIC_FLOAT_PRECISION> refVD = readDataset3D(importFile, dataPath3D);
    Array4D<PRISMATIC_FLOAT_PRECISION> refCBED = readDataset4D(importFile, dataPath4D);

    Array2D<PRISMATIC_FLOAT_PRECISION> testAnnular = readDataset2D(meta.filenameOutput, dataPath2D);
    Array3D<PRISMATIC_FLOAT_PRECISION> testPS = readDataset3D(meta.filenameOutput, dataPathPS);
    Array3D<PRISMATIC_FLOAT_PRECISION> testDPC = readDataset3D(meta.filenameOutput, dataPathDPC);
    Array3D<PRISMATIC_FLOAT_PRECISION> testVD = readDataset3D(meta.filenameOutput, dataPath3D);
    Array4D<PRISMATIC_FLOAT_PRECISION> testCBED = readDataset4D(meta.filenameOutput, dataPath4D);

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
    
    divertOutput(pos, fd, logPath);
    std::cout << "\n##### BEGIN TEST CASE: fourierResampling ######\n";

    std::string importFile = "../test/potentialImport.h5";
    meta.filenameOutput = "../test/potentialImport.h5";
    go(meta);

    std::cout << "\n--------------------------------------------\n";

    //force 
    meta.interpolationFactorX = 5;
    meta.interpolationFactorY = 7;
    meta.filenameOutput = "../test/potentialRerun.h5";
    meta.importFile     = "../test/potentialImport.h5";
    meta.importPath     = "4DSTEM_simulation/data/realslices/ppotential/realslice";
    meta.importPotential = true;
    go(meta);
    std::cout << "####### END TEST CASE: fourierResampling ######\n";

    revertOutput(fd, pos);

    //read in output arrays and compare
    std::string dataPath2D = "4DSTEM_simulation/data/realslices/annular_detector_depth0000/realslice";
    std::string dataPathDPC = "4DSTEM_simulation/data/realslices/DPC_CoM_depth0000/realslice";
    std::string dataPath3D = "4DSTEM_simulation/data/realslices/virtual_detector_depth0000/realslice";
    std::string dataPathPS = "4DSTEM_simulation/data/realslices/ppotential/realslice";

    Array2D<PRISMATIC_FLOAT_PRECISION> refAnnular = readDataset2D(importFile, dataPath2D);
    Array3D<PRISMATIC_FLOAT_PRECISION> refPS = readDataset3D(importFile, dataPathPS);
    Array3D<PRISMATIC_FLOAT_PRECISION> refDPC = readDataset3D(importFile, dataPathDPC);
    Array3D<PRISMATIC_FLOAT_PRECISION> refVD = readDataset3D(importFile, dataPath3D);

    std::cout << refPS.get_dimi() << " " << refPS.get_dimj() << " " << refPS.get_dimk() << std::endl; 

    Array2D<PRISMATIC_FLOAT_PRECISION> testAnnular = readDataset2D(meta.filenameOutput, dataPath2D);
    Array3D<PRISMATIC_FLOAT_PRECISION> testPS = readDataset3D(meta.filenameOutput, dataPathPS);
    Array3D<PRISMATIC_FLOAT_PRECISION> testDPC = readDataset3D(meta.filenameOutput, dataPathDPC);
    Array3D<PRISMATIC_FLOAT_PRECISION> testVD = readDataset3D(meta.filenameOutput, dataPath3D);

    PRISMATIC_FLOAT_PRECISION tol = 0.001;
    PRISMATIC_FLOAT_PRECISION errorSum = 0.0;

    //only compare subset of outputs
    BOOST_TEST(compareSize(refAnnular, testAnnular));
    BOOST_TEST(compareSize(refVD, testVD));
    BOOST_TEST(compareSize(refDPC, testDPC));

    BOOST_TEST(compareValues(refAnnular, testAnnular) < tol);
    BOOST_TEST(compareValues(refDPC, testDPC) < tol);
    BOOST_TEST(compareValues(refVD, testVD) < tol);

    //resampled potential slices should have 80 x 84 x 3 dims do align with fx = 5, fy = 7
    //rewrite to auto check multiples of 4x(fx, fy)
    std::cout << testPS.get_dimi() << " " << testPS.get_dimj() << " " << std::endl;
    bool sizeCheck = testPS.get_dimi() == 80 && testPS.get_dimj() == 84 && testPS.get_dimk() == 3;
    BOOST_TEST(sizeCheck);

    // removeFile(importFile);
    // removeFile(meta.filenameOutput);

    unsigned long  seven = 7;
    unsigned long eightfour = 84;
    std::cout << -1 % seven << std::endl;
    std::cout << -1 % eightfour << std::endl;
};

BOOST_FIXTURE_TEST_CASE(importSMatrix, basicSim)
{

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

BOOST_AUTO_TEST_CASE(complexIO)
{
    //testing IO operations on complex datasets
    PRISMATIC_FLOAT_PRECISION dummy = 1.0; //dummy float for IO overlaoding
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

    writeComplexDataset(testGroup, "complex2D", &refArr2D[0], mdims_2D, 2);
    writeComplexDataset(testGroup, "complex3D", &refArr3D[0], mdims_3D, 3);
    writeComplexDataset(testGroup, "complex4D", &refArr4D[0], mdims_4D, 4);

    testGroup.close();
    testFile.close();

    std::string datapath2D = "complex_data/complex2D";
    std::string datapath3D = "complex_data/complex3D";
    std::string datapath4D = "complex_data/complex4D";

    Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> testArr2D; 
    Array3D<std::complex<PRISMATIC_FLOAT_PRECISION>> testArr3D; 
    Array4D<std::complex<PRISMATIC_FLOAT_PRECISION>> testArr4D; 
    readComplexDataset(testArr2D, fname, datapath2D);
    readComplexDataset(testArr3D, fname, datapath3D);
    readComplexDataset(testArr4D, fname, datapath4D);

    PRISMATIC_FLOAT_PRECISION tol = 0.0001;
    BOOST_TEST(compareSize(refArr2D, testArr2D));
    BOOST_TEST(compareSize(refArr3D, testArr3D));
    BOOST_TEST(compareSize(refArr4D, testArr4D));
    BOOST_TEST(compareValues(refArr2D, testArr2D) < tol);
    BOOST_TEST(compareValues(refArr3D, testArr3D) < tol);
    BOOST_TEST(compareValues(refArr4D, testArr4D) < tol);

}

BOOST_AUTO_TEST_SUITE_END();


}