//
// Created by AJ Pryor on 1/27/17.
//
// This code creates a 1-D array and doubles all elements separately on the GPU and CPU using
// the Thrust library, then verifies the result

#include <iostream>
#include <algorithm>
#include "thrust/host_vector.h"
#include "thrust/device_vector.h"

#define ARR_SIZE 1000

__global__ void times_two(double* d_a, size_t N){
    auto idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N) d_a[idx] *= 2;
}


// this function checks if the elements of two vectors are equal
template <typename T>
bool all_equal(typename thrust::host_vector<T>::iterator a_first,
               typename thrust::host_vector<T>::iterator a_last,
               typename thrust::host_vector<T>::iterator b_first){
    while (a_first != a_last){
        if (*a_first != *b_first){
            return false;
        }
        ++a_first;
        ++b_first;
    };
    return true;
};

using namespace std;
int main(){
    cout << "Testing a basic kernel that doubles values in an array using Thrust" << endl;

    // allocate a vector that is equivalent to calling MATLAB 1:ARR_SIZE
    thrust::host_vector<double> h_a(ARR_SIZE,0);
    thrust::transform(h_a.begin(), h_a.end(), h_a.begin(), [](const double& a){
        static double counter = 0;
        return counter++;
    });

    // allocate memory on the GPU
    thrust::device_vector<double> d_a = h_a;

    // perform the double on the device in CUDA
    times_two<<<1, ARR_SIZE>>>(thrust::raw_pointer_cast(d_a.data()), d_a.size());

    // copy the result back to host
    thrust::host_vector<double> d_answer = d_a;

    // perform the doubling on the host side in C++
    for (auto &i:h_a)i*=2;

    // for (auto &i:h_a)cout << i << endl;

    // print the first value of both arrays
    cout << "First 3 elements of device array: " << d_answer[0] << ", "
                                                 << d_answer[1] << ", "
                                                 << d_answer[2] << ", "
                                                 << endl;

    cout << "First 3 elements of host array: "   << h_a[0] << ", "
                                                 << h_a[1] << ", "
                                                 << h_a[2] << ", "
                                                 << endl;

    // make sure the GPU is finished
    cudaDeviceSynchronize();

    // compare results
    cout << "Verifying all elements are equal..." << endl;
    if (all_equal<double>(h_a.begin(),h_a.end(),d_answer.begin()))cout << "CPU and GPU results are equal" << endl;
}