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
#include <chrono>

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

void removeFile(const std::string &filepath);

BOOST_GLOBAL_FIXTURE(logFile);

BOOST_AUTO_TEST_SUITE(ioSpeedTests);

BOOST_AUTO_TEST_CASE(arrayWrites)
{
    //prepare file and datasets for IO operations
    std::string fname = "../test/arrayWrites.h5";
    H5::H5File file = H5::H5File(fname.c_str(), H5F_ACC_TRUNC);
    H5::Group dsgroup = file.createGroup("/datasets");

    int seed = 10101;
    srand(seed);
    std::default_random_engine de(seed);

    size_t N_runs = 100;
    std::vector<size_t> dim_sizes = {10, 50, 100, 250, 1000};
    size_t num_dims = 4;
    for(auto i = 0; i < num_dims; i++)
    {
        for(auto j = 0; j < num_dims; j++)
        {
            size_t Nx = dim_sizes[i];
            size_t Ny = dim_sizes[j];
            std::cout << "2D array: " << Nx << " x " << Ny << std::endl;
            std::string dsname = "testds" + std::to_string(i) + "_" + std::to_string(j);
            Array2D<PRISMATIC_FLOAT_PRECISION> testArr = zeros_ND<2, PRISMATIC_FLOAT_PRECISION>({{Ny, Nx}});
            assignRandomValues(testArr, de);

            hsize_t dims2D[2] = {Nx, Ny};
            H5::DataSpace mspace(2, dims2D);
            H5::DataSet testds = dsgroup.createDataSet(dsname.c_str(), PFP_TYPE, mspace);
            mspace.close();
            testds.close();

            std::vector<size_t> order = {0,1};

            //baseline fastest speed
            auto t1 = std::chrono::high_resolution_clock::now();
            for(auto n = 0; n < N_runs; n++)
            {
                writeRealDataSet_inOrder(dsgroup, dsname, &testArr[0], dims2D, 2);
            }
            auto t2 = std::chrono::high_resolution_clock::now();

            double baseline = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();


            // //with element by element selection
            // t1 = std::chrono::high_resolution_clock::now();
            // for(auto n = 0; n < N_runs; n++)
            // {
            //     writeRealDataSet(dsgroup, dsname, &testArr[0], dims2D, 2, order);
            // }
            // t2 = std::chrono::high_resolution_clock::now();

            // double strat1 = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

            //restride array in memory, then write in order
            t1 = std::chrono::high_resolution_clock::now();
            for(auto n = 0; n < N_runs; n++)
            {
                std::array<size_t, 2> dims_in = {Ny, Nx};
                std::array<size_t, 2> dims_order = {0,1};
                Array2D<PRISMATIC_FLOAT_PRECISION> writeArr = restride(testArr, dims_in, dims_order);
                writeRealDataSet_inOrder(dsgroup, dsname, &writeArr[0], dims2D, 2);
            }
            t2 = std::chrono::high_resolution_clock::now();

            double strat2 = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
  
            //more optimal restride
            t1 = std::chrono::high_resolution_clock::now();
            for(auto n = 0; n < N_runs; n++)
            {
                Array2D<PRISMATIC_FLOAT_PRECISION> writeArr = zeros_ND<2, PRISMATIC_FLOAT_PRECISION>({{testArr.get_dimi(), testArr.get_dimj()}});
                for(auto jj = 0; jj < testArr.get_dimj(); jj++)
                {
                    for(auto ii = 0; ii < testArr.get_dimi(); ii++)
                    {
                        writeArr.at(ii,jj) = testArr.at(jj,ii);
                    }
                }
                writeRealDataSet_inOrder(dsgroup, dsname, &writeArr[0], dims2D, 2);
            }
            t2 = std::chrono::high_resolution_clock::now();

            double strat3 = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

            std::cout << "average baseline write time: " << baseline/(1000*(double) N_runs) << " milliseconds" << std::endl;
            // std::cout << "  average strat1 write time: " << strat1/(1000*(double) N_runs) << " milliseconds" << std::endl;
            std::cout << "  average strat2 write time: " << strat2/(1000*(double) N_runs) << " milliseconds" << std::endl;
            std::cout << "  average strat3 write time: " << strat3/(1000*(double) N_runs) << " milliseconds" << std::endl;
        }
    }

    for(auto i = 0; i < num_dims; i++)
    {
        for(auto j = 0; j < num_dims; j++)
        {
            for(auto k = 0; k < num_dims; k++)
            {

                size_t Nx = dim_sizes[i];
                size_t Ny = dim_sizes[j];
                size_t Nz = dim_sizes[k];
                std::cout << "3D array: " << Nx << " x " << Ny << " x " << Nz << std::endl;
                std::string dsname = "testds" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k);
                Array3D<PRISMATIC_FLOAT_PRECISION> testArr = zeros_ND<3, PRISMATIC_FLOAT_PRECISION>({{Nz, Ny, Nx}});
                assignRandomValues(testArr, de);

                hsize_t dims3D[3] = {Nx, Ny, Nz};
                H5::DataSpace mspace(3, dims3D);
                H5::DataSet testds = dsgroup.createDataSet(dsname.c_str(), PFP_TYPE, mspace);
                mspace.close();
                testds.close();

                std::vector<size_t> order = {0,1, 2};

                //baseline fastest speed
                auto t1 = std::chrono::high_resolution_clock::now();
                for(auto n = 0; n < N_runs; n++)
                {
                    writeRealDataSet_inOrder(dsgroup, dsname, &testArr[0], dims3D, 3);
                }
                auto t2 = std::chrono::high_resolution_clock::now();

                double baseline = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();


                // //with element by element selection
                // t1 = std::chrono::high_resolution_clock::now();
                // for(auto n = 0; n < N_runs; n++)
                // {
                //     writeRealDataSet(dsgroup, dsname, &testArr[0], dims3D, 3, order);
                // }
                // t2 = std::chrono::high_resolution_clock::now();

                // double strat1 = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

                //restride array in memory, then write in order
                t1 = std::chrono::high_resolution_clock::now();
                for(auto n = 0; n < N_runs; n++)
                {
                    std::array<size_t, 3> dims_in = {Nz, Ny, Nx};
                    std::array<size_t, 3> dims_order = {0,1,2};
                    Array3D<PRISMATIC_FLOAT_PRECISION> writeArr = restride(testArr, dims_in, dims_order);
                    writeRealDataSet_inOrder(dsgroup, dsname, &writeArr[0], dims3D, 3);
                }
                t2 = std::chrono::high_resolution_clock::now();

                double strat2 = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
                
                //more optimal restride
                t1 = std::chrono::high_resolution_clock::now();
                for(auto n = 0; n < N_runs; n++)
                {
                    Array3D<PRISMATIC_FLOAT_PRECISION> writeArr = zeros_ND<3, PRISMATIC_FLOAT_PRECISION>({{testArr.get_dimi(), testArr.get_dimj(), testArr.get_dimk()}});
                    for(auto kk = 0; kk < testArr.get_dimk(); kk++)
                    {
                        for(auto jj = 0; jj < testArr.get_dimj(); jj++)
                        {
                            for(auto ii = 0; ii < testArr.get_dimi(); ii++)
                            {
                                writeArr.at(ii,jj,kk) = testArr.at(kk,jj,ii);
                            }
                        }
                    }
                    writeRealDataSet_inOrder(dsgroup, dsname, &writeArr[0], dims3D, 3);
                }
                t2 = std::chrono::high_resolution_clock::now();

                double strat3 = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();


                //keep array in original order

                std::cout << "average baseline write time: " << baseline/(1000*(double) N_runs) << " milliseconds" << std::endl;
                // std::cout << "  average strat1 write time: " << strat1/(1000*(double) N_runs) << " milliseconds" << std::endl;
                std::cout << "  average strat2 write time: " << strat2/(1000*(double) N_runs) << " milliseconds" << std::endl;
                std::cout << "  average strat3 write time: " << strat3/(1000*(double) N_runs) << " milliseconds" << std::endl;
            }
        }
    }

    // removeFile(fname);
}

BOOST_AUTO_TEST_CASE(arrayReads)
{
        //prepare file and datasets for IO operations
    std::string fname = "../test/arrayWrites.h5";

    size_t N_runs = 100;
    std::vector<size_t> dim_sizes = {10, 50, 100, 250, 1000};
    size_t num_dims = 4;
    for(auto i = 0; i < num_dims; i++)
    {
        for(auto j = 0; j < num_dims; j++)
        {
            size_t Nx = dim_sizes[i];
            size_t Ny = dim_sizes[j];
            std::cout << "2D array: " << Nx << " x " << Ny << std::endl;
            std::string dsname = "/datasets/testds" + std::to_string(i) + "_" + std::to_string(j);
            Array2D<PRISMATIC_FLOAT_PRECISION> testArr;
            std::vector<size_t> order = {0,1};

            //baseline fastest speed
            auto t1 = std::chrono::high_resolution_clock::now();
            for(auto n = 0; n < N_runs; n++)
            {
                readRealDataSet_inOrder(testArr, fname, dsname);
            }
            auto t2 = std::chrono::high_resolution_clock::now();

            double baseline = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();


            //with element by element selection
            t1 = std::chrono::high_resolution_clock::now();
            for(auto n = 0; n < N_runs; n++)
            {
                readRealDataSet(testArr, fname, dsname, order);
            }
            t2 = std::chrono::high_resolution_clock::now();

            double strat1 = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

            //more optimal restride
            t1 = std::chrono::high_resolution_clock::now();
            for(auto n = 0; n < N_runs; n++)
            {
                readRealDataSet_inOrder(testArr, fname, dsname);
                Array2D<PRISMATIC_FLOAT_PRECISION> finalArr = zeros_ND<2, PRISMATIC_FLOAT_PRECISION>({{testArr.get_dimi(), testArr.get_dimj()}});
                for(auto jj = 0; jj < testArr.get_dimj(); jj++)
                {
                    for(auto ii = 0; ii < testArr.get_dimi(); ii++)
                    {
                        finalArr.at(ii,jj) = testArr.at(jj,ii);
                    }
                }

            }
            t2 = std::chrono::high_resolution_clock::now();

            double strat3 = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

            std::cout << "average baseline read time: " << baseline/(1000*(double) N_runs) << " milliseconds" << std::endl;
            std::cout << "  average strat1 read time: " << strat1/(1000*(double) N_runs) << " milliseconds" << std::endl;
            std::cout << "  average strat3 read time: " << strat3/(1000*(double) N_runs) << " milliseconds" << std::endl;
        }
    }

    for(auto i = 0; i < num_dims; i++)
    {
        for(auto j = 0; j < num_dims; j++)
        {
            for(auto k = 0; k < num_dims; k++)
            {
                size_t Nx = dim_sizes[i];
                size_t Ny = dim_sizes[j];
                size_t Nz = dim_sizes[k];
                std::cout << "3D array: " << Nx << " x " << Ny << " x " << Nz << std::endl;
                std::string dsname = "/datasets/testds" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k);
                Array3D<PRISMATIC_FLOAT_PRECISION> testArr;


                std::vector<size_t> order = {0,1,2};

                //baseline fastest speed
                auto t1 = std::chrono::high_resolution_clock::now();
                for(auto n = 0; n < N_runs; n++)
                {
                    readRealDataSet_inOrder(testArr, fname, dsname);
                }
                auto t2 = std::chrono::high_resolution_clock::now();

                double baseline = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();


                // //with element by element selection
                // t1 = std::chrono::high_resolution_clock::now();
                // for(auto n = 0; n < N_runs; n++)
                // {
                //     readRealDataSet(testArr, fname, dsname, order);
                // }
                // t2 = std::chrono::high_resolution_clock::now();

                // double strat1 = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

                //more optimal restride
                t1 = std::chrono::high_resolution_clock::now();
                for(auto n = 0; n < N_runs; n++)
                {
                    readRealDataSet_inOrder(testArr, fname, dsname);
                    Array3D<PRISMATIC_FLOAT_PRECISION> finalArr = zeros_ND<3, PRISMATIC_FLOAT_PRECISION>({{testArr.get_dimi(), testArr.get_dimj(), testArr.get_dimk()}});
                    for(auto kk = 0; kk < testArr.get_dimk(); kk++)
                    {
                        for(auto jj = 0; jj < testArr.get_dimj(); jj++)
                        {
                            for(auto ii = 0; ii < testArr.get_dimi(); ii++)
                            {
                                finalArr.at(ii,jj,kk) = testArr.at(kk,jj,ii);
                            }
                        }
                    }
                }
                t2 = std::chrono::high_resolution_clock::now();

                double strat3 = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

                std::cout << "average baseline read time: " << baseline/(1000*(double) N_runs) << " milliseconds" << std::endl;
                // std::cout << "  average strat1 read time: " << strat1/(1000*(double) N_runs) << " milliseconds" << std::endl;
                std::cout << "  average strat3 read time: " << strat3/(1000*(double) N_runs) << " milliseconds" << std::endl;
            }
        }
    }

}

BOOST_AUTO_TEST_SUITE_END();
}