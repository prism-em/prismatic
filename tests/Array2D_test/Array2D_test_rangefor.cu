//
// Created by AJ Pryor on 1/30/17.
//
//Test range-based for loops on GPU with PRISM::Array2D

#include <iostream>
#include "Array2D.h"
#include "thrust/host_vector.h"
#include "thrust/device_vector.h"
#include "../../../../../../../usr/local/cuda/include/thrust/detail/copy.h"

using namespace std;
int main(){


    thrust::host_vector<int> test_h(6,0);
    for (auto i = 0; i < test_h.size(); ++i)test_h[i] = i + 1;
    PRISM::Array2D< thrust::host_vector<int> > arr_h(test_h,2,3);
    cout << "Printing a host-side 2x3 PRISM::Array2D with thrust::host_vector" << endl;
    for (auto& i:arr_h) std::cout << i << std::endl;


    thrust::device_vector<int> test_d(test_h);
    PRISM::Array2D< thrust::device_vector<int> > arr_d(test_d,2,3);
    thrust::copy(arr_d.begin(), arr_d.end(),arr_h.begin());
    cout << "Printing a device-side 2x3 PRISM::Array2D with thrust::device_vector" << endl;
    for (auto& i:arr_h) std::cout << i << std::endl;

    cout << "Test Successful!" << endl;

}