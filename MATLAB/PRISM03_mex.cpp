#include <vector>
#include <stdexcept>
#include <string>
#include "fftw3.h"
#include "mex.h"

// define the position of each input here for clarity, since there are many
#define POS_probeDefocusArray 0
#define POS_probeSemiangleArray 1
#define POS_probeXtiltArray 2
#define POS_probeYtiltArray 3
#define POS_qxaReduce 4
#define POS_qyaReduce 5
#define POS_xp 6
#define POS_yp 7
#define POS_imageSize 8
#define POS_imageSizeReduce 9
#define POS_scale 10
#define POS_interpolationFactor 11
#define POS_pixelSize 12
#define POS_lambda 13
#define NUM_INPUTS 13

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){
    
    // Check we for correct number of arguments
    if (nrhs != NUM_INPUTS){
        throw std::invalid_argument("Incorrect number of inputs to mex function, should be " + std::to_string(NUM_INPUTS));
    }
    
    // get points to array variables
    double *probeDefocusArray   = mxGetPr(prhs[POS_probeDefocusArray]);
    double *probeSemiangleArray = mxGetPr(prhs[POS_probeSemiangleArray]);
    double *probeXtiltArray     = mxGetPr(prhs[POS_probeXtiltArray]);
    double *probeYtiltArray     = mxGetPr(prhs[POS_probeYtiltArray]);
    double *qxaReduce           = mxGetPr(prhs[POS_qxaReduce]);
    double *qyaReduce           = mxGetPr(prhs[POS_qyaReduce]);
    double *xp                  = mxGetPr(prhs[POS_xp]);
    double *yp                  = mxGetPr(prhs[POS_yp]);
    double *imageSize           = mxGetPr(prhs[POS_imageSize]);
    double *imageSizeReduce     = mxGetPr(prhs[POS_imageSizeReduce]);
    
    // get scalars
    double scale                = mxGetScalar(prhs[POS_scale]);
    double interpolationFactor  = mxGetScalar(prhs[POS_interpolationFactor]);
    double pixelSize            = mxGetScalar(prhs[POS_pixelSize]);
    double lambda               = mxGetScalar(prhs[POS_lambda]);
      
    size_t mrows = mxGetM(prhs[0]);
    size_t ncols = mxGetN(prhs[0]);
    mexPrintf("mrows = %d\n", mrows);
    mexPrintf("mcols = %d\n", ncols);

//     plhs[0] = mxCreateDoubleMatrix((mwSize)mrows, (mwSize)ncols, mxREAL);
//     double v[mrows*ncols];
    double* v = (double*)mxMalloc(mrows*ncols*sizeof(double));
    for (int i = 0; i < mrows*ncols; ++i){
//         v[i] = i;
        v[i] = probeDefocusArray[i];
        mexPrintf("val %f\n",v[i]); 
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