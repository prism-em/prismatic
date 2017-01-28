// This is currently a placeholder program that just computes an FFT

#include <iostream>
#include <complex>
#include "fftw3.h"

#define N 1000
int main() {
    double vec[N];
    std::complex<double> k[N];
    for (auto i=0; i<N; ++i) vec[i]=i;
    fftw_plan plan = fftw_plan_dft_r2c_1d(N,
                                          vec,
                                          reinterpret_cast<fftw_complex*>(k),
                                          FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    for (auto i = 0,j=0; j<10;i+=5,++j){
        std::cout << "k[" <<  i << "] =" << k[i] << std::endl;

    }
    return 0;
}