#include <vector>
#include "fftw3.h"
#include "mex.h"

#define _probeDefocusArray 1
#define _probeSemiangleArray 1
#define _probeXtiltArray 1
#define _probeYtiltArray 1
#define _qxaReduce 1
#define _qyaReduce 1
#define _xp 1
#define _yp 1
#define _imageSize 1
#define _imageSizeReduce 1
#define _scale 1
#define _interpolationFactor 1
#define _pixelSize 1
#define _lambda 1

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){
    double *probeDefocusArray = mxGetPr(prhs[_probeDefocusArray]);
    mexPrintf("hey there\n");
//     plhs[0] = mxDuplicateArray(prhs[0]);
    
    size_t mrows = mxGetM(prhs[0]);
    size_t ncols = mxGetN(prhs[0]);
    mexPrintf("mrows = %d\n", mrows);
    mexPrintf("mcols = %d\n", ncols);

//     plhs[0] = mxCreateDoubleMatrix((mwSize)mrows, (mwSize)ncols, mxREAL);
//     double v[mrows*ncols];
    double* v = (double*)mxMalloc(mrows*ncols*sizeof(double));
    for (int i = 0; i < mrows*ncols; ++i){
        v[i] = (double)i;
    };
    plhs[0] = mxCreateNumericMatrix(0,0,mxDOUBLE_CLASS,mxREAL);

    mxSetPr(plhs[0], v); 
    mxSetM(plhs[0],mrows);
    mxSetN(plhs[0],ncols);
    
    
//     plhs[0] = mxGetPr(prhs[0]);
    
     
      //initialize:
   /*
    q1 = zeros(imageSizeReduce);
    q2 = zeros(imageSizeReduce);
    PsiProbeInit = zeros(imageSizeReduce);
    psi = zeros(imageSizeReduce);
    intOutput = zeros(imageSizeReduce);
     */
    //dq = mean([qxaReduce(2,1) qyaReduce(1,2)]);
    //qProbeMax = emdSTEM.probeSemiangleArray(a1) / emdSTEM.lambda;
}