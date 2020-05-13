#include "ioTests.h"
#include <boost/test/unit_test.hpp>
#include "ArrayND.h"
#include <iostream>
#include <vector>
#include "go.h"
#include "meta.h"
#include "params.h"
#include <stdio.h>

namespace Prismatic{

class basicSim{

    public:
    basicSim()     {setupSim(),BOOST_TEST_MESSAGE( "Setting up fixture");}
    ~basicSim()    {BOOST_TEST_MESSAGE( "Tearing down fixture");}
    Metadata<PRISMATIC_FLOAT_PRECISION> meta;

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
    }
    

};

BOOST_AUTO_TEST_SUITE(ioTests);

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

BOOST_AUTO_TEST_CASE(readH5)
{
    BOOST_TEST(1==0);
};

BOOST_AUTO_TEST_SUITE_END();


}