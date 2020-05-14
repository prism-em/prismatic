#ifndef IO_TESTS_H
#define IO_TESTS_H
#include <random>
#include "ArrayND.h"
#include <vector>
#include "params.h"

namespace Prismatic{


void assignRandomValues(Array1D<PRISMATIC_FLOAT_PRECISION> &arr, std::default_random_engine &de)
{
    std::normal_distribution<PRISMATIC_FLOAT_PRECISION> randn(0,1);
    for (auto i = 0; i < arr.size(); i++) arr[i] = randn(de);
};

void assignRandomValues(Array2D<PRISMATIC_FLOAT_PRECISION> &arr, std::default_random_engine &de)
{
    std::normal_distribution<PRISMATIC_FLOAT_PRECISION> randn(0,1);
    for (auto i = 0; i < arr.size(); i++) arr[i] = randn(de);
};

void assignRandomValues(Array3D<PRISMATIC_FLOAT_PRECISION> &arr, std::default_random_engine &de)
{
    std::normal_distribution<PRISMATIC_FLOAT_PRECISION> randn(0,1);
    for (auto i = 0; i < arr.size(); i++) arr[i] = randn(de);
};

void assignRandomValues(Array4D<PRISMATIC_FLOAT_PRECISION> &arr, std::default_random_engine &de)
{
    std::normal_distribution<PRISMATIC_FLOAT_PRECISION> randn(0,1);
    for (auto i = 0; i < arr.size(); i++) arr[i] = randn(de);
};

} //namespace Prismatic
#endif