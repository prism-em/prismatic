//
// Created by AJ Pryor on 1/27/17.
//
//
//Test construction of a PRISM:Array2D with Thrust vectors
#include <iostream>
#include "Array2D.h"
#include "thrust/host_vector.h"
#include "thrust/device_vector.h"

using namespace std;
int main(){
    cout << "Creating a host-side 2x3 PRISM::Array2D with std::vector" << endl;

    vector<int> test{1,2,3,4,5,6};
    PRISM::Array2D< std::vector<int> > arr(test,2,3);
    cout << "nrows = " << arr.get_nrows() << endl;
    cout << "ncols = " << arr.get_ncols() << endl;

    cout << "Creating a host-side 2x3 PRISM::Array2D with thrust::host_vector" << endl;
    thrust::host_vector<int> test_h(6,0);
    PRISM::Array2D< thrust::host_vector<int> > arr_h(test_h,2,3);
    cout << "nrows = " << arr.get_nrows() << endl;
    cout << "ncols = " << arr.get_ncols() << endl;

    cout << "Creating a device-side 2x3 PRISM::Array2D with thrust::device_vector" << endl;
    thrust::device_vector<int> test_d(6,0);
    PRISM::Array2D< thrust::device_vector<int> > arr_d(test_d,2,3);
    cout << "nrows = " << arr.get_nrows() << endl;
    cout << "ncols = " << arr.get_ncols() << endl;
    cout << "Test Successful!" << endl;
}
