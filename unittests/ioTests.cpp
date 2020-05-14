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

namespace Prismatic{

class basicSim{

    public:
    basicSim()     {setupSim(),BOOST_TEST_MESSAGE( "Setting up fixture");}
    ~basicSim()    {BOOST_TEST_MESSAGE( "Tearing down fixture");}
    Metadata<PRISMATIC_FLOAT_PRECISION> meta;
    Parameters<PRISMATIC_FLOAT_PRECISION> pars;

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

BOOST_AUTO_TEST_SUITE(ioTests);

/*
BOOST_FIXTURE_TEST_CASE(operationReorganization, basicSim)
{
    //make sure nothing is broken by moving all file IO operations to their own source
    std::cout << "Opening log file ioTests.log to capture simulation output." << std::endl;
    std::string logFile = "ioTests.log";
    fpos_t pos;
    fgetpos(stdout, &pos);
    int fd = dup(fileno(stdout));
    freopen(logFile.c_str(),"w",stdout);
    go(meta);

    std::cout << "\n--------------------------------------------\n";
    meta.algorithm = Algorithm::Multislice;
    go(meta);
    //clean up files
    fflush(stdout);
    dup2(fd,fileno(stdout));
    close(fd);
    clearerr(stdout);
    fsetpos(stdout, &pos);

    std::cout << "Log file closed. Output returning to terminal." << std::endl;
    if( remove( meta.filenameOutput.c_str() ) != 0 )
        perror( "Error deleting file" );
    else
        puts( "Test file successfully deleted" );

};
*/

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
    bool sizeCheck2D = read2D.size()==testArr2D.size() 
                        && read2D.get_dimi() == testArr2D.get_dimi() 
                        && read2D.get_dimj() == testArr2D.get_dimj();

    BOOST_TEST(sizeCheck2D);
            
    bool sizeCheck3D = read3D.size()==testArr3D.size() 
                        && read3D.get_dimi() == testArr3D.get_dimi() 
                        && read3D.get_dimj() == testArr3D.get_dimj()
                        && read3D.get_dimk() == testArr3D.get_dimk();

    BOOST_TEST(sizeCheck3D);

    bool sizeCheck4D = read4D.size()==testArr4D.size() 
                        && read4D.get_dimi() == testArr4D.get_dimi() 
                        && read4D.get_dimj() == testArr4D.get_dimj()
                        && read4D.get_dimk() == testArr4D.get_dimk()
                        && read4D.get_diml() == testArr4D.get_diml();

    BOOST_TEST(sizeCheck4D);

    //check value equivalence
    PRISMATIC_FLOAT_PRECISION tol = 0.00001;
    PRISMATIC_FLOAT_PRECISION errorSum = 0.0;
    for(auto i = 0; i < read2D.size(); i++) errorSum += std::abs(read2D[i]-testArr2D[i]);
    BOOST_TEST(errorSum < tol);

    errorSum = 0.0;
    for(auto i = 0; i < read3D.size(); i++) errorSum += std::abs(read3D[i]-testArr3D[i]);
    BOOST_TEST(errorSum < tol);

    errorSum = 0.0;
    for(auto i = 0; i < read4D.size(); i++) errorSum += std::abs(read4D[i]-testArr4D[i]);
    BOOST_TEST(errorSum < tol);

    if( remove( pars.meta.filenameOutput.c_str() ) != 0 )
        perror( "Error deleting file" );
    else
        puts( "Test file successfully deleted" );

};

BOOST_AUTO_TEST_SUITE_END();


}