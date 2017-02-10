#include <vector>
#include <stdexcept>
#include <string>
#include "fftw3.h"
#include "mex.h"
#include "Array2D.h"

// define the position of each input here for clarity, since there are so many
#define POS_Scompact            0
#define POS_stack               1
#define POS_probeDefocusArray   2
#define POS_probeSemiangleArray 3
#define POS_probeXtiltArray     4
#define POS_probeYtiltArray     5
#define POS_qxaReduce           6
#define POS_qyaReduce           7
#define POS_xp                  8
#define POS_yp                  9
#define POS_beamsIndex          10
#define POS_xyBeams             11
#define POS_xVec                12
#define POS_yVec                13
#define POS_imageSizeOutput     14
#define POS_detectorAngles      15
#define POS_cellDim             16
#define POS_scale               17
#define POS_lambda              18
#define POS_dr                  19
#define POS_dq                  20
#define POS_Ndet                21
#define POS_numFP               22
#define NUM_INPUTS              22

template <typename T>
PRISM::Array2D< std::vector<T> > mat2DtoPRISM2D(const mxArray *array){
    size_t nrows = mxGetM(array);
    size_t ncols = mxGetN(array);
    mexPrintf("nrows = %i\n",nrows);
    mexPrintf("ncols = %i\n",ncols);
    size_t N     = nrows*ncols;
    double *ptr  = mxGetPr(array);
    PRISM::Array2D< std::vector<T> > arr( std::vector<T>(N, 0), nrows, ncols);
    for (auto& i:arr)i=(T)*ptr++;
    for (auto& i:arr)mexPrintf("%f\n",i);
    return arr;
}


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){
    
    // Check we for correct number of arguments
//     if (nrhs != NUM_INPUTS)mexErrMsgTxt("Incorrect number of inputs to mex function, should be 22");
    using PRISM_FLOAT_TYPE = double;
    PRISM::Array2D< std::vector<PRISM_FLOAT_TYPE> > qxaReduce = mat2DtoPRISM2D<double>(prhs[0]); // change from 0

    for (auto &i : qxaReduce)mexPrintf("%f\n",i);
    plhs[0] = mxCreateDoubleMatrix(qxaReduce.get_nrows(), qxaReduce.get_ncols(), mxREAL);
    
    double * ptr = mxGetPr(plhs[0]);
    for (auto &i:qxaReduce)*ptr++=i;
//    // get points to array variables
//    double *probeDefocusArray   = mxGetPr(prhs[POS_probeDefocusArray]);
//    double *probeSemiangleArray = mxGetPr(prhs[POS_probeSemiangleArray]);
//    double *probeXtiltArray     = mxGetPr(prhs[POS_probeXtiltArray]);
//    double *probeYtiltArray     = mxGetPr(prhs[POS_probeYtiltArray]);
//    double *qxaReduce           = mxGetPr(prhs[POS_qxaReduce]);
//    double *qyaReduce           = mxGetPr(prhs[POS_qyaReduce]);
//    double *xp                  = mxGetPr(prhs[POS_xp]);
//    double *yp                  = mxGetPr(prhs[POS_yp]);
//    double *imageSize           = mxGetPr(prhs[POS_imageSize]);
//    double *imageSizeReduce     = mxGetPr(prhs[POS_imageSizeReduce]);
//
//    // get scalars
//    double scale                = mxGetScalar(prhs[POS_scale]);
//    double interpolationFactor  = mxGetScalar(prhs[POS_interpolationFactor]);
//    double pixelSize            = mxGetScalar(prhs[POS_pixelSize]);
//    double lambda               = mxGetScalar(prhs[POS_lambda]);
//
//    size_t mrows = mxGetM(prhs[0]);
//    size_t ncols = mxGetN(prhs[0]);
//    mexPrintf("mrows = %d\n", mrows);
//    mexPrintf("mcols = %d\n", ncols);
//
////     plhs[0] = mxCreateDoubleMatrix((mwSize)mrows, (mwSize)ncols, mxREAL);
////     double v[mrows*ncols];
//    double* v = (double*)mxMalloc(mrows*ncols*sizeof(double));
//    for (int i = 0; i < mrows*ncols; ++i){
////         v[i] = i;
//        v[i] = probeDefocusArray[i];
//        mexPrintf("val %f\n",v[i]);
//    };
//    plhs[0] = mxCreateNumericMatrix(0,0,mxDOUBLE_CLASS,mxREAL);
//
//    mxSetPr(plhs[0], v);
//    mxSetM(plhs[0],mrows);
//    mxSetN(plhs[0],ncols);
//
//
////     plhs[0] = mxGetPr(prhs[0]);
//
//
//      //initialize:
//   /*
//    q1 = zeros(imageSizeReduce);
//    q2 = zeros(imageSizeReduce);
//    PsiProbeInit = zeros(imageSizeReduce);
//    psi = zeros(imageSizeReduce);
//    intOutput = zeros(imageSizeReduce);
//     */
//    //dq = mean([qxaReduce(2,1) qyaReduce(1,2)]);
//    //qProbeMax = emdSTEM.probeSemiangleArray(a1) / emdSTEM.lambda;
}