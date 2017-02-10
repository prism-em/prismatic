#include <iostream>
#include <vector>
#include "mex.h"

void timesTwo(double* a, const int& N);


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){
    double* a = mxGetPr(prhs[0]);
size_t N = (size_t)(mxGetN(prhs[0]) * mxGetM(prhs[0]));
mexPrintf("a[1] = %f\n",a[1]);
    timesTwo(a, N);   
    mexPrintf("a[1] = %f\n",a[1]);
    plhs[0] = mxCreateDoubleMatrix(1,N,mxREAL);
    double *o = mxGetPr(plhs[0]);
    for (int i = 0; i < N; ++i)o[i] = a[i];
}
