//
// Created by AJ Pryor on 1/27/17.
//
// This code creates a 1-D array and doubles all elements separately on the GPU and CPU using
// low-level CUDA API calls, then verifies the result

#include <iostream>
#include <algorithm>

#define ARR_SIZE 1000

__global__ void times_two(double* d_a, size_t N){
    auto idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N) d_a[idx] *= 2;
}


// this function checks if the elements of two vectors are equal
template <typename T>
bool all_equal(T* a_first,
               T* a_last,
               T* b_first){
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
    cout << "Testing a basic kernel that doubles values in an array using low-level CUDA API calls" << endl;

    // allocate a vector that is equivalent to calling MATLAB 1:ARR_SIZE
    double* h_a = new double[ARR_SIZE];

    std::transform(h_a, h_a + ARR_SIZE, h_a, [](const double& a){
        static double counter = 0;
        return counter++;
    });


    // allocate memory on the GPU
    double* d_a;
    cudaMalloc((void**)&d_a, ARR_SIZE * sizeof(double));
    double* d_answer = new double[ARR_SIZE];
    //copy from host to device
    cudaMemcpy(d_a, h_a, ARR_SIZE * sizeof(double), cudaMemcpyHostToDevice);

    // perform the double on the device in CUDA
    times_two<<<1, ARR_SIZE>>>(d_a, ARR_SIZE);

    // copy the result back to host
    cudaMemcpy(d_answer, d_a, ARR_SIZE * sizeof(double), cudaMemcpyDeviceToHost);

    // perform the doubling on the host side in C++
    for (double* i = h_a; i != h_a + ARR_SIZE; ++i)*i *=2;

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
    if (all_equal<double>(h_a, h_a + ARR_SIZE, d_answer)) cout << "CPU and GPU results are equal" << endl;
    cout << "Test Successful!" << endl;
    // free memory
    delete(h_a);
    delete(d_answer);
    cudaFree(d_a);

 }
